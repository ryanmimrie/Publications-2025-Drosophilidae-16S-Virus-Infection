# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure 5 & 6 =================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs and plots individual bacterial genera changes in
# abundance across Drosophilidae host species.

# Due to contamination in samples of D. takahashii and Z. davidi, apparent in
# RNAseq runs of the same samples used for 16S analysis, these species have been
# removed.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(MCMCglmm)
library(phyloseq)
library(microbiome)
library(here)
library(ape)
library(patchwork)
library(vegan)
#library(Cairo)
library(zCompositions)
library(compositions)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ----- 0.3. Load Data ---------------------------------------------------------

load(here("data", "data_phyloseq_cutoff.RData"))

data_metadata <- read_csv(here("data", "data_sampleMetadata.csv"))

data_fly <- read_csv(here("data", "data_drosophilidae.csv"))

tree_drosophilidae <- read.tree(here("data", "tree_hosts.nwk"))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1 Aggregate by Genus -------------------------------------------------

data_abundances <- aggregate_taxa(phyloseq_cutoff, level="Genus")

# ----- 1.2 Remove non-experimental samples ------------------------------------

data_abundances <- subset_samples(data_abundances, !(ID %in% c(as.character(307:324), "Negative")))
data_abundances <- subset_samples(data_abundances, !(Spc.Name %in% c("D. takahashii", "Z. davidi")))

# ----- 1.3. CLR Transform -----------------------------------------------------

data_presence <- as(otu_table(data_abundances), "matrix") %>% t() %>% as.data.frame()

data_presence$ID <- rownames(data_presence) %>% str_remove_all("_R1.fastq.gz") %>% as.numeric() %>% as.character()
data_presence <- gather(data_presence, key = "Genus", value = "presence", 1:(ncol(data_presence)-1))
data_presence$presence <- data_presence$presence > 0

data_clr <- as(otu_table(data_abundances), "matrix") %>% t() %>%
  cmultRepl(method = "CZM", output = "p-counts", z.warning = 1, z.delete = F)

data_clr <- clr(data_clr) %>% as.data.frame()

data_clr$ID <- rownames(data_clr) %>% str_remove_all("_R1.fastq.gz") %>% as.numeric() %>% as.character()

data_clr <- gather(data_clr, key = "Genus", value = "abundance", 1:(ncol(data_clr)-1))

data_abundances <- as.data.frame(as.matrix(sample_data(data_abundances)))

data_clr <- left_join(data_clr, data_abundances)

data_clr <- left_join(data_clr, data_presence)

data_clr$animal <- data_clr$Spc.Name %>% str_remove_all(" ")

data_clr$Genus <- ifelse(data_clr$Genus == "Escherichia-Shigella", "Escherichia", data_clr$Genus)

data_clr <- left_join(data_clr, data_fly)

# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ----- 1.3 Wrangle To Wide Format for MCMCglmm --------------------------------

data_clr_wide_presence <- select(data_clr, ID, Condition, animal, Genus, presence, Diet)
data_clr_wide_presence <- spread(data_clr_wide_presence, key = Genus, value = presence)

data_clr_wide_abundance <- select(data_clr, ID, Condition, animal, Genus, abundance, Diet)
data_clr_wide_abundance <- spread(data_clr_wide_abundance, key = Genus, value = abundance)

data_clr_wide_presence$spc <- data_clr_wide_presence$animal
data_clr_wide_abundance$spc <- data_clr_wide_abundance$animal

# ------------------------------------------------------------------------------
# ----- 2. Abundance Models ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Model Prior --------------------------------------------------------

prior_pr <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  #G2 = list(V = diag(1), nu = 0.002),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(1), nu = 0.002, fix = 1))

prior_ab <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  #G2 = list(V = diag(1), nu = 0.002),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(1), nu = 0.002))

# ----- 2.2 Iteration Multiplier -----------------------------------------------

itt <- 1 # Set to 1 for demonstration, 100 for publication

# ----- 2.3 Model Loop ---------------------------------------------------------

data_contrasts <- data.frame(Genus = NA,
                             Condition = NA,
                             Value_pr = NA,
                             Value_pr_low = NA,
                             Value_pr_high = NA,
                             Value_ab = NA,
                             Value_ab_low = NA,
                             Value_ab_high = NA,
                             Contrast_pr = NA,
                             Contrast_pr_low = NA,
                             Contrast_pr_high = NA,
                             Contrast_pr_prob = NA,
                             Contrast_pr_prob_low = NA,
                             Contrast_pr_prob_high = NA,
                             Contrast_ab = NA,
                             Contrast_ab_low = NA,
                             Contrast_ab_high = NA)[-1,]

holder <- data.frame(hold = NA)

for(Genus in unique(data_clr$Genus)){
  
  fixedeffs <- as.formula(paste(Genus, "~ Condition"))
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("None", "Ringers", "DCV"))
  data_clr_wide_abundance$Condition <- factor(data_clr_wide_abundance$Condition, levels = c("None", "Ringers", "DCV"))
  
  data_clr_wide_abundance_notzero <- data_clr_wide_abundance[data_clr_wide_presence[[Genus]] == TRUE, ]
  
  
  if (!file.exists(here("models", sprintf("%s_NoneContrast_presence.Rdata", Genus)))) {
    
    save(holder, file = here("models", sprintf("%s_NoneContrast_presence.Rdata", Genus)))
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_pr,
                         data = data_clr_wide_presence,
                         family = "threshold",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_pr, file = here("models", sprintf("%s_NoneContrast_presence.Rdata", Genus)))
    
  } else {load(here("models", sprintf("%s_NoneContrast_presence.Rdata", Genus)))}
  
  if (!file.exists(here("models", sprintf("%s_NoneContrast_abundance.Rdata", Genus)))) {
    
    save(holder, file = here("models", sprintf("%s_NoneContrast_abundance.Rdata", Genus)))
    
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_ab,
                         data = data_clr_wide_abundance_notzero,
                         family = "gaussian",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_ab, file = here("models", sprintf("%s_NoneContrast_abundance.Rdata", Genus)))
    
  } else {load(here("models", sprintf("%s_NoneContrast_abundance.Rdata", Genus)))}
  
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("Ringers", "None", "DCV"))
  data_clr_wide_abundance_notzero$Condition <- factor(data_clr_wide_abundance_notzero$Condition, levels = c("Ringers", "None", "DCV"))
  
  if (!file.exists(here("models", sprintf("%s_RingersContrast_presence.Rdata", Genus)))) {
    
    save(holder, file = here("models", sprintf("%s_RingersContrast_presence.Rdata", Genus)))
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_pr,
                         data = data_clr_wide_presence,
                         family = "threshold",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_pr, file = here("models", sprintf("%s_RingersContrast_presence.Rdata", Genus)))
    
  } else {load(here("models", sprintf("%s_RingersContrast_presence.Rdata", Genus)))}
  
  
  if (!file.exists(here("models", sprintf("%s_RingersContrast_abundance.Rdata", Genus)))) {
    
    save(holder, file = here("models", sprintf("%s_RingersContrast_abundance.Rdata", Genus)))
    
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_ab,
                         data = data_clr_wide_abundance_notzero,
                         family = "gaussian",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_ab, file = here("models", sprintf("%s_RingersContrast_abundance.Rdata", Genus)))
    
  } else {load(here("models", sprintf("%s_RingersContrast_abundance.Rdata", Genus)))}
  
}

