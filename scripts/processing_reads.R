# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: 16S Read Processing ==========
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(dada2)
library(Biostrings)
library(decontam)
library(microbiome)
library(phyloseq)
library(microViz)
library(vegan)
library(openssl)
library(ape)
library(here)
library(progress)

# ------------------------------------------------------------------------------
# ----- 1. Read Filtering ------------------------------------------------------
# ------------------------------------------------------------------------------

# NOTE: The raw files used in this section are not included in the GITHUB
#       repo due to file size limits, but can be downloaded from the NCBI
#       sequence read archive (SRA).

#       The fastq.gz files output from this section are included in the GITHUB
#       repo, so users without the raw read files can start from section 2.

# ----- 1.1. Quality Profiles --------------------------------------------------

read_files <- list.files(here("reads", "raw"), pattern = "fastq.gz", full.names = TRUE)

sample_ID <- read_files %>% basename() %>% str_split("_", simplify = TRUE) %>% .[, 1]

direction <- read_files %>% basename()
direction <- sub(".*_(R[12])\\.fastq", "\\1", direction)

read_files <- data.frame(filename = read_files,
                         sample = sample_ID,
                         direction = direction,
                         filesize = file.size(read_files),
                         stringsAsFactors = FALSE)

read_paths_F <- sort(list.files(here("reads", "raw"), pattern = "R1.fastq.gz", full.names = TRUE))
read_paths_R <- sort(list.files(here("reads", "raw"), pattern = "R2.fastq.gz", full.names = TRUE))

pb <- progress_bar$new(format = "  Generating quality profile plots :current/:total [:bar] :percent | ETA: :eta",
                       total = length(read_paths_F), clear = FALSE, width = 80)

for (i in seq_along(read_paths_F)) {
  pb$tick()
  out_file <- here("reads", "QC", sprintf("%s.png", i))
  if (!file.exists(out_file)) {
    plot <- plotQualityProfile(c(read_paths_F[i], read_paths_R[i]))
    ggsave(out_file, plot, width = 10, height = 8, dpi = 300)
  }
}

# ----- 1.2. Filter Reads ------------------------------------------------------

read_paths_F_filtered <- here("reads", "filtered", basename(read_paths_F))
read_paths_R_filtered <- here("reads", "filtered", basename(read_paths_R))

# ETA: 30 minutes
filtered_read_info <- filterAndTrim(read_paths_F, read_paths_F_filtered, read_paths_R, read_paths_R_filtered, truncLen=c(240,225),
                                    maxEE=c(2,2))

# ------------------------------------------------------------------------------
# ----- 2. Error Correction ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Learn Error Rates -------------------------------------------------

read_paths_F_filtered <- sort(list.files(here("reads", "filtered"), pattern = "R1.fastq.gz", full.names = TRUE))
read_paths_R_filtered <- sort(list.files(here("reads", "filtered"), pattern = "R2.fastq.gz", full.names = TRUE))

# ETA: 20 minutes
errors_F <- learnErrors(read_paths_F_filtered, nbases = 1e+09)
errors_R <- learnErrors(read_paths_R_filtered, nbases = 1e+09)

plotErrors(errors_F, nominalQ=TRUE)
plotErrors(errors_R, nominalQ=TRUE)

# ----- 2.2. Denoise Reads -----------------------------------------------------

# ETA: 30 minutes
reads_F_denoised <- dada(read_paths_F_filtered, err=errors_F)
reads_R_denoised <- dada(read_paths_R_filtered, err=errors_R)

# ----- 2.3. Merge Reads -----------------------------------------------------

reads_merged <- mergePairs(reads_F_denoised, read_paths_F_filtered, reads_R_denoised, read_paths_R_filtered, verbose=TRUE)

# ----- 2.4. Filter Merged Reads by Length -------------------------------------

sequence_table <- makeSequenceTable(reads_merged)
dim(sequence_table)

table(nchar(getSequences(sequence_table)))

sequence_table_filtered <- sequence_table[,nchar(colnames(sequence_table)) %in% 400:450]

# ----- 2.5. Remove Chimeras ---------------------------------------------------

runif(1, 0, 100000) # returned 74963.76
set.seed(74963.76)

# ETA: 2 minutes
sequence_table_no_chimera <- removeBimeraDenovo(sequence_table_filtered, method="consensus", verbose=TRUE)

# ----- 2.6. Read N Check ------------------------------------------------------

getN <- function(x) sum(getUniques(x))

if (exists("filtered_read_info")){
  read_deltas <- cbind(filtered_read_info, sapply(reads_F_denoised, getN), sapply(reads_R_denoised, getN), sapply(reads_merged, getN), rowSums(sequence_table_no_chimera))
  colnames(read_deltas) <- c("Input", "Filtered", "Denoised_F", "Denoised_R", "Merged", "Non_Chimeric")
} else {
  read_deltas <- cbind(sapply(reads_F_denoised, getN), sapply(reads_R_denoised, getN), sapply(reads_merged, getN), rowSums(sequence_table_no_chimera))
  colnames(read_deltas) <- c("Denoised_F", "Denoised_R", "Merged", "Non_Chimeric")
}

rownames(read_deltas) <- basename(read_paths_F_filtered) %>% str_remove_all("_R1.fastq.gz")
head(read_deltas)

# ----- 2.7. Assign Taxonomy (Silva v138.1) ------------------------------------

runif(1, 0, 100000) # returned 47454.27
set.seed(47454.27)

# Reconstruct split files (workaround for GitHub file size constraints)
if (!file.exists(here("data", "silva_nr99_v138.1_train_set.fa"))) {
  part_aa <- readLines(gzfile(here::here("data", "silva_nr99_v138.1_train_set.part-aa.gz")))
  part_ab <- readLines(gzfile(here::here("data", "silva_nr99_v138.1_train_set.part-ab.gz")))
  silva <- c(part_aa, part_ab)
  writeLines(silva, here("data", "silva_nr99_v138.1_train_set.fa"))
  rm(list = c("silva", "part_aa", "part_ab"))
  gc()
}
if (!file.exists(here("data", "silva_species_assignment_v138.1.fa"))) {
  part_aa <- readLines(gzfile(here::here("data", "silva_species_assignment_v138.1.part-aa.gz")))
  part_ab <- readLines(gzfile(here::here("data", "silva_species_assignment_v138.1.part-ab.gz")))
  silva <- c(part_aa, part_ab)
  writeLines(silva, here("data", "silva_species_assignment_v138.1.fa"))
  rm(list = c("silva", "part_aa", "part_ab"))
  gc()
}

# ETA: 10 minutes
taxa <- assignTaxonomy(sequence_table_no_chimera, here("data", "silva_nr99_v138.1_train_set.fa"))
taxa <- addSpecies(taxa, here("data", "silva_species_assignment_v138.1.fa"))

head(unname(taxa))
apply(taxa,2,function(x)mean(!is.na(x)))

# ----- 2.8. Tidy Up Sequence Names --------------------------------------------

sequences <- DNAStringSet(getSequences(sequence_table_no_chimera))

hash_names <- md5(getSequences(sequence_table_no_chimera))

colnames(sequence_table_no_chimera) <- hash_names
names(sequences) <- hash_names
reads <- rownames(taxa)
rownames(taxa) <- hash_names

reads <- reads[!is.na(taxa[, 6])]
taxa <- taxa[!is.na(taxa[, 6]), ]

sequence_table_no_chimera <- sequence_table_no_chimera[,colnames(sequence_table_no_chimera) %in% rownames(taxa)]
sequences <- sequences[names(sequences) %in% rownames(taxa)]

# ----- 2.9. Write OTU Data ----------------------------------------------------

data_OTUs <- data.frame(name = row.names(taxa),
                        sequence = reads,
                        length = nchar(row.names(taxa)),
                        kingdom = taxa[,1],
                        phylum = taxa[,2],
                        class = taxa[,3],
                        order = taxa[,4],
                        family = taxa[,5],
                        genus = taxa[,6],
                        species = taxa[,7])

write_csv(data_OTUs, here("data", "data_OTUs.csv"))

lines <- sapply(1:nrow(data_OTUs), function(i) {
  c(paste0(">", data_OTUs$name[i]), data_OTUs$sequence[i])})

writeLines(unlist(lines), here("data", "data_OTUs.fas"))

# ------------------------------------------------------------------------------
# ----- 3. Create Phyloseq Object ----------------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1. External Phylogeny Process ----------------------------------------

# The phylogenetic tree of OTUs was created as follows:
# 1. Alignment was created using MAFFT with default settings
# 2. Uninformative gap regions were removed by trimAL -gt 0.9
# 3. Trees were inferred from the alignment using FastTree with default settings

tree <- read.tree(here("data", "data_OTUs_tree.nwk"))

# ----- 3.2. Add Sample Metadata -----------------------------------------------

data_meta <- as.data.frame(read_csv(here("data", "data_metadata.csv")))
data_meta <- rbind(data_meta, c("Negative", rep(NA, 14), "ESS Negative"))

ID_to_sample <- data.frame(sample = sample_names(otu_table(t(sequence_table_no_chimera), taxa_are_rows = TRUE)),
                            ID = c(1:324, "Negative"))

data_meta <- left_join(data_meta, ID_to_sample)

spc <- read_csv(here("data", "data_drosophilidae.csv"))

data_meta <- left_join(data_meta, spc)

rownames(data_meta) <- data_meta$sample

# ----- 3.3. Generate Phyloseq Object ------------------------------------------

phyloseq_raw <- phyloseq(otu_table(t(sequence_table_no_chimera), taxa_are_rows = TRUE),
                      sample_data(data_meta),
                      tax_table(taxa),
                      phy_tree(tree),
                      sequences)


# ------------------------------------------------------------------------------
# ----- 4. Filter Phyloseq Object ----------------------------------------------
# ------------------------------------------------------------------------------

# ----- 4.1. Contaminant Removal -----------------------------------------------

# Removes taxa NAs
phyloseq_filtered <-prune_taxa(as.vector(!is.na(tax_table(phyloseq_raw)[,2])),phyloseq_raw)

# Removes chloroplasts
phyloseq_filtered<-prune_taxa(as.vector(tax_table(phyloseq_filtered)[,3]!="Chloroplast"),phyloseq_filtered)

# Removes archaea
phyloseq_filtered<-prune_taxa(as.vector(tax_table(phyloseq_filtered)[,1]!="Archaea"),phyloseq_filtered)

# Estimate and remove technical contaminants
sample_data(phyloseq_filtered)$is_negative <- sample_data(phyloseq_filtered)$Type == "Negative control"
contamdf.prev <- isContaminant(phyloseq_filtered, method="prevalence", neg="is_negative",threshold=0.1)
table(contamdf.prev$contaminant)
phyloseq_filtered<-prune_taxa(contamdf.prev$contaminant==FALSE,phyloseq_filtered)
sample_data(phyloseq_filtered)$size <- sample_sums(phyloseq_filtered)

save(phyloseq_filtered, file = here("data", "data_phyloseq_filtered.RData"))
load(here("data", "data_phyloseq_filtered.RData"))

# ----- 4.2. Minimum Read Depth Cutoff -----------------------------------------

cutoffs <- seq(1250, 1500, 50)

df_list <- list()

for (cutoff in cutoffs) {
  pruned_data <- prune_samples(sample_data(phyloseq_filtered)$size >= cutoff, phyloseq_filtered)
  df <- as.data.frame(table(sample_data(pruned_data)$Spc.Name.Short, sample_data(pruned_data)$Condition))
  df$cutoff <- cutoff
  df_list[[as.character(cutoff)]] <- df
}

df_all <- do.call(rbind, df_list)

df_all$Freq <- factor(df_all$Freq, levels = c("3", "2", "1"))
df_all$Var2 <- factor(df_all$Var2, levels = c("None", "Ringers", "DCV"))

ggplot(df_all) +
  geom_tile(aes(x = Var2, y = Var1, fill = Freq)) +
  facet_grid(cols = vars(cutoff)) +
  theme(axis.title = element_blank()) +
  scale_fill_manual(name = "Replicates", values = c("#1abc9c", "#f1c40f", "#e74c3c"))

# Minimum read depth cutoff of 1350 is the highest round value that doesn't
# reduce an experimental condition to 1 replicate

phyloseq_filtered <- prune_samples(sample_data(phyloseq_filtered)$size >= 1350, phyloseq_filtered)

# ----- 4.3. Minimum Prevalence Cutoff -----------------------------------------

phyloseq_cutoff <- tax_filter(phyloseq_filtered, min_prevalence = 0.07)

sort(unique(tax_table(phyloseq_cutoff)[, "Genus"]))

save(phyloseq_cutoff, file = here("data", "data_phyloseq_cutoff.RData"))


