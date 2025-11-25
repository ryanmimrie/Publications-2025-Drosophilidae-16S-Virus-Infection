# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure 1 =====================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script plots bacterial genera relative abundances across 32
# Drosophilidae host species and three experimental conditions: No injection,
# saline (Ringers) injection, and DCV infection.

# Due to contamination in samples of D. takahashii and Z. davidi, apparent in
# RNAseq runs of the same samples used for 16S analysis, these species have been
# removed.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(plotrix)
library(scales)
library(patchwork)
library(ggtree)
library(ape)
library(here)
library(Cairo)

# ----- 0.3. Load Data ---------------------------------------------------------

load(here("data", "data_phyloseq_cutoff.RData"))

data_viralLoad <- read_csv(here("data", "data_viralLoad.csv"))

tree_drosophilidae <- read.tree(here("data", "tree_hosts.nwk"))

tree_bacteria <- read.nexus(here("data", "tree_16S.nex"))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Phyloseq object to mean abundances --------------------------------

phyloseq_cutoff <- aggregate_taxa(phyloseq_cutoff, level="Genus")

data <- phyloseq_cutoff %>% tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
  dplyr::select(Genus, Sample, Abundance) %>% spread(Sample, Abundance)

data <- gather(data, key = "ID", value = "abundance", 2:ncol(data))

data$ID <- str_remove_all(data$ID, "_R1.fastq.gz") %>% as.numeric() %>% as.character() # NAs are the ESS negative

data <- left_join(data, phyloseq_cutoff@sam_data)

data <- na.omit(data)

data_sum <- data %>%
  group_by(Spc.Name.Short, Condition, Genus) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop")

data_sum$species_condition <- paste(data_sum$Spc.Name.Short, data_sum$Condition, sep = "_")

data_sum$Spc.Name.Short <- NULL
data_sum$Condition <- NULL

data_sum <- data_sum %>%
  separate(species_condition, into = c("species", "condition"), sep = "_")

data_sum$condition <- factor(data_sum$condition, levels = c("None", "Ringers", "DCV"))

data_sum <- filter(data_sum, !species %in% c("Dtak", "Zdav"))

data_sum$label <- ifelse(data_sum$condition == "None", "No Injection",
                         ifelse(data_sum$condition == "Ringers", "Saline Injection (Ringers)", "DCV Infection"))

data_sum$label <- factor(data_sum$label, levels = c("No Injection", "Saline Injection (Ringers)", "DCV Infection"))

data_sum$mean_abundance <- log10(data_sum$mean_abundance*100)

data_sum$mean_abundance <- ifelse(data_sum$mean_abundance == -Inf | data_sum$mean_abundance < -2, -2, data_sum$mean_abundance)

# ----- 1.2. Add host and bacteria species ladderised order --------------------

species <- data.frame(species = unique(data_sum$species),
                      full = c("D. affinis", "D. americana", "D. ananassae",
                               "D. arizonae", "D. baimaii", "D. buzzatii",
                               "D. euronotus", "D. flavomontana", "D. hydei",
                               "D. immigrans", "D. lacicola", "D. melanogaster",
                               "D. montana", "D. nasuta", "D. nebulosa",
                               "D. paramelanica", "D. persimilis", "D. prosaltans",
                               "D. pseudoobscura", "D. putrida", "D. saltans",
                               "D. santomea", "D. simulans", "D. sturtevanti",
                               "D. subobscura", "D. sucinea",
                               "D. teissieri", "D. virilis", "D. yakuba",
                               "S. lebanonensis", "S. pattersoni",
                               "Z. tuberculatus"))

data_sum <- left_join(data_sum, species)

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, str_remove_all(data_sum$full, " ")))

p_fly <- ggtree(tree_drosophilidae) +
  geom_tiplab() +
  coord_cartesian(xlim = c(0, 1))


data_sum$full <- factor(data_sum$full, levels = c("Z. tuberculatus", "D. putrida",
                                                  "D. immigrans", "D. nasuta", "D. euronotus",
                                                  "D. paramelanica", "D. hydei", "D. arizonae",
                                                  "D. buzzatii", "D. virilis", "D. americana",
                                                  "D. lacicola", "D. flavomontana", "D. montana",
                                                  "S. lebanonensis", "S. pattersoni", "D. nebulosa",
                                                  "D. sucinea", "D. sturtevanti", "D. prosaltans",
                                                  "D. saltans", "D. subobscura", "D. affinis",
                                                  "D. persimilis", "D. pseudoobscura", "D. baimaii",
                                                  "D. ananassae", "D. melanogaster",
                                                  "D. simulans", "D. teissieri", "D. santomea",
                                                  "D. yakuba"))

p_bac <- ggtree(tree_bacteria) +
  geom_tiplab() +
  coord_cartesian(xlim = c(0, 1))

bacteria <- c("Levilactobacillus", "Fructilactobacillus", "Companilactobacillus", "Lactobacillus",
  "Enterococcus", "Streptococcus", "Staphylococcus", "Bacillus",
  "Finegoldia", "Lawsonella", "Corynebacterium", "Rhodococcus",
  "Cutibacterium", "Afipia", "Bosea", "Brucella",
  "Brevundimonas", "Enhydrobacter", "Commensalibacter", "Acetobacter",
  "Gluconobacter", "Buttiauxella", "Providencia", "Serratia",
  "Escherichia-Shigella", "Pseudomonas", "Ralstonia", "Pelomonas",
  "Stenotrophomonas", "Dysgonomonas", "Bacteroides", "Sediminibacterium")

data_sum$Genus <- factor(data_sum$Genus, levels = rev(bacteria))

# ------------------------------------------------------------------------------
# ----- 2. Figure 1 ------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Figure 1 Heatmap --------------------------------------------------

plot <- ggplot(data_sum) +
  geom_tile(aes(x = full, y = Genus, fill = mean_abundance)) +
  coord_fixed(ratio = 1) +
  facet_wrap(~label, ncol = 3) +
  scale_fill_gradientn(colours = c('#000005', '#080716', '#110b2d', '#1e0848', '#300060',
                                   '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                   '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                   '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                       breaks = c(log10(1:10)-4, log10(1:10)-3, log10(1:10)-2, log10(1:10)-1, log10(1:10), log10(1:10)+1), limits = c(-2, 2),
                       labels = c(rep("", 9), "0.001%", rep("", 9), "<=0.01%", rep("", 9), "0.1%", rep("", 9), "1%", rep("", 9), "10%", rep("", 8), "100%", ""),
                       name = "Relative\nAbundance", guide = guide_colorbar(barheight = 10, barwidth = 0.75, title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(position = "right", expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 7, color = "#2e3440"),
        axis.text.x = element_text(angle = 90, hjust = 0, face = "italic", vjust = 0.5),
        axis.text.y = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(angle = 0, vjust = 0.5),
        legend.position = "right",
        strip.background = element_rect(fill = "#e5e9f0"))

ggsave(here("figures", "Figure_1_raw_hmp.png"), plot, dpi = 300, width = 8, height = 4)
ggsave(here("figures", "Figure_1_raw_hmp.svg"), plot, dpi = 300, width = 8, height = 4, device = cairo_pdf)


# ----- 2.2. Figure 1 Host Phylogeny -------------------------------------------

plot_drosophilidae <- ggtree(tree_drosophilidae) +
  geom_tiplab() +
  geom_treescale(width = 0.05) +
  coord_cartesian(xlim = c(0, 1))

ggsave(here("figures", "Figure_1_raw_phylo_dros.png"), plot_drosophilidae, dpi = 300, width = 3, height = 6)
ggsave(here("figures", "Figure_1_raw_phylo_dros.svg"), plot_drosophilidae, dpi = 300, width = 3, height = 6, device = cairo_pdf)

# ----- 2.3. Figure 1 Bacteria Phylogeny ---------------------------------------

plot_bacteria <- ggtree(tree_bacteria) +
  geom_tiplab() +
  geom_treescale(width = 0.05) +
  coord_cartesian(xlim = c(0, 1))

ggsave(here("figures", "Figure_1_raw_phylo_bact.png"), plot_bacteria, dpi = 300, width = 3, height = 6)
ggsave(here("figures", "Figure_1_raw_phylo_bact.svg"), plot_bacteria, dpi = 300, width = 3, height = 6, device = cairo_pdf)

# ----- 2.4. Figure 1 Viral Load Bars ------------------------------------------

data_viralLoad_sum <- data_viralLoad %>% group_by(species) %>%
  summarise(fold = mean(foldchange),
            se = std.error(foldchange))

data_viralLoad_sum <- filter(data_viralLoad_sum, !species %in% c("Z. davidi", "D. takahashii"))

data_viralLoad_sum$species <- factor(data_viralLoad_sum$species,
                                     levels = levels(data_sum$full))

shift <- 1.2

viral <- ggplot(data_viralLoad_sum) +
  geom_bar(aes(x = species, y = fold+shift, fill = fold), stat = "identity") +
  geom_errorbar(aes(x = species, ymin = fold+shift-se, ymax = fold+shift+se), width =0.5, size = 0.4, color = "#2e3440") +
  geom_hline(aes(yintercept = shift), linetype = "dotted") +
  theme_bw() +
  theme(text = element_text(size = 13, colour = "#2e3440"),
        axis.text.x = element_text(angle = 90, hjust = 0, face = "italic", vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(angle = 270, vjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_continuous(position = "right", name = "Viral\nLoad (Δ)", expand = c(0,0), limits = c(0, 9), breaks = c(shift, 3 + shift, 6 + shift),
                     labels = c(expression(10^0),
                                expression(10^3),
                                expression(10^6))) +
  scale_fill_gradientn(colours = c('#000005', '#080716', '#110b2d', '#1e0848', '#300060',
                                   '#43006a', '#57096e', '#6b116f', '#81176d', '#991d69',
                                   '#b02363', '#ca2d59', '#e03b50', '#ee504a', '#f66b4d',
                                   '#fa8657', '#fca368', '#fcc17d', '#fcdf96', '#fbffb2'),
                       name = "Viral\nLoad (Δ)", breaks = c(0, 3, 6),
                       labels = c(expression(10^0),
                                  expression(10^3),
                                  expression(10^6)),
                       guide = guide_colorbar(barheight = 2.5, barwidth = 1, title.position = "right", title.hjust = 0.5))

ggsave(here("figures", "Figure_1_raw_bars.png"), plot = viral, width = 6.3, height = 1.8)
ggsave(here("figures", "Figure_1_raw_bars.svg"), plot = viral, width = 6.27, height = 1.8, device = cairo_pdf)

