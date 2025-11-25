# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figures 7 & 8 ================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script models the relationship between individual bacterial
# abundances and viral load.

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

tree_drosophilidae <- read.tree(here("data", "tree_hosts.nwk"))

data_viralLoad <- read_csv(here("data", "data_viralLoad.csv"))

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

colnames(data_viralLoad) <- c("animal", "foldchange", "Block")
data_viralLoad$Block <- as.character(data_viralLoad$Block)
data_viralLoad$animal <- data_viralLoad$animal %>% str_remove_all(" ")

data_clr$animal <- data_clr$Spc.Name %>% str_remove_all(" ")

data_clr <- left_join(data_clr, data_viralLoad)

data_clr$Genus <- ifelse(data_clr$Genus == "Escherichia-Shigella", "Escherichia", data_clr$Genus)


# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ------------------------------------------------------------------------------
# ----- 2. Abundance Models ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Model Prior --------------------------------------------------------

prior_phylo_presence <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(2), nu = 0.002, fix = 2))

prior_phylo_abundance <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(2), nu = 0.002))

prior_nonphylo_presence <- list(R = list(V = diag(2), nu = 0.002, fix = 2))

prior_nonphylo_abundance <- list(R = list(V = diag(2), nu = 0.002))

# ----- 2.2 Iteration Multiplier -----------------------------------------------

itt <- 100 # Set to 1 for demonstration, 100 for publication

# ----- 2.3 Model Loop ---------------------------------------------------------

data_correlations <- data.frame(Genus = NA,
                                           Condition = NA,
                                           B_pr_phylo = NA,
                                           B_pr_phylo_low = NA,
                                           B_pr_phylo_high = NA,
                                           R_pr_phylo = NA,
                                           R_pr_phylo_low = NA,
                                           R_pr_phylo_high = NA,
                                           B_ab_phylo = NA,
                                           B_ab_phylo_low = NA,
                                           B_ab_phylo_high = NA,
                                           R_ab_phylo = NA,
                                           R_ab_phylo_low = NA,
                                           R_ab_phylo_high = NA,
                                           C_ab_phylo = NA,
                                           C_ab_phylo_low = NA,
                                           C_ab_phylo_high = NA,
                                           B_pr_nonphylo = NA,
                                           B_pr_nonphylo_low = NA,
                                           B_pr_nonphylo_high = NA,
                                           R_pr_nonphylo = NA,
                                           R_pr_nonphylo_low = NA,
                                           R_pr_nonphylo_high = NA,
                                           B_ab_nonphylo = NA,
                                           B_ab_nonphylo_low = NA,
                                           B_ab_nonphylo_high = NA,
                                           R_ab_nonphylo = NA,
                                           R_ab_nonphylo_low = NA,
                                           R_ab_nonphylo_high = NA,
                                           C_ab_nonphylo = NA,
                                           C_ab_nonphylo_low = NA,
                                           C_ab_nonphylo_high = NA)[-1,]

data_clr <- filter(data_clr, Genus %in% unique(data_clr$Genus)[1:16])

holder <- data.frame(hold = NA)

for(g in unique(data_clr$Genus)){
  for(c in c("None", "DCV")){
    
    data <- filter(data_clr, Genus == g, Condition == c)
    data <- select(data, abundance, presence, foldchange, animal, ID)
    
    data_pres <- select(data, presence, foldchange, animal, ID)
    data_abun <- filter(data, presence) %>% select(abundance, foldchange, animal, ID)
    
    data_pres <- gather(data_pres, key = "trait", value = "value", 1:2)
    data_abun <- gather(data_abun, key = "trait", value = "value", 1:2)
    
    data_pres$family <- ifelse(data_pres$trait == "presence", "threshold", "gaussian")
    data_abun$family <- "gaussian"
    
    data_pres$trait <- factor(data_pres$trait, levels = c("foldchange", "presence"))
    data_abun$trait <- factor(data_abun$trait, levels = c("foldchange", "abundance"))
    
    data_pres <- as.data.frame(data_pres)
    data_abun <- as.data.frame(data_abun)
    
    ab <- filter(data_abun, trait == "abundance")
    
    if (length(unique(ab$animal))>=1){

      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_phylo_presence.Rdata", g, c)))) {
      
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_phylo_presence.Rdata", g, c)))
        
        model_phylo_presence <- MCMCglmm(value ~ trait,
                                         random = ~us(trait):animal + ID,
                                         rcov = ~idh(trait):units,
                                         pedigree = tree_drosophilidae,
                                         family = NULL,
                                         prior = prior_phylo_presence,
                                         data = data_pres,
                                         nitt = 130000*itt,
                                         thin = 10*itt,
                                         burnin = 30000*itt)
      
        save(model_phylo_presence, file = here("models", sprintf("%s_%s_ViralLoad_phylo_presence.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_phylo_presence.Rdata", g, c)))}
      
      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_phylo_abundance.Rdata", g, c)))) {
        
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_phylo_abundance.Rdata", g, c)))
        
        model_phylo_abundance <- MCMCglmm(value ~ trait,
                                          random = ~us(trait):animal + ID,
                                          rcov = ~idh(trait):units,
                                          pedigree = tree_drosophilidae,
                                          family = NULL,
                                          prior = prior_phylo_abundance,
                                          data = data_abun,
                                          nitt = 130000*itt,
                                          thin = 10*itt,
                                          burnin = 30000*itt)
        
        save(model_phylo_abundance, file = here("models", sprintf("%s_%s_ViralLoad_phylo_abundance.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_phylo_abundance.Rdata", g, c)))}
      
      data_pres$family <- NULL
      data_abun$family <- NULL
      
      data_pres <- spread(data_pres, key = "trait", value = "value")
      data_abun <- spread(data_abun, key = "trait", value = "value")
      
      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_nonphylo_presence.Rdata", g, c)))) {
        
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_nonphylo_presence.Rdata", g, c)))
        
        model_nonphylo_presence <- MCMCglmm(cbind(foldchange, presence) ~ trait - 1,
                                            rcov = ~ us(trait):units,
                                            family = c("gaussian", "threshold"),
                                            prior = prior_nonphylo_presence,
                                            data = data_pres,
                                            nitt = 130000 * itt,
                                            thin = 10 * itt,
                                            burnin = 30000 * itt)
        
        save(model_nonphylo_presence, file = here("models", sprintf("%s_%s_ViralLoad_nonphylo_presence.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_nonphylo_presence.Rdata", g, c)))}
      
      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_nonphylo_abundance.Rdata", g, c)))) {
        
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_nonphylo_abundance.Rdata", g, c)))
        
        model_nonphylo_abundance <- MCMCglmm(cbind(foldchange, abundance) ~ trait - 1,
                                             rcov = ~ us(trait):units,
                                             family = c("gaussian", "gaussian"),
                                             prior = prior_nonphylo_abundance,
                                             data = data_abun,
                                             nitt = 130000 * itt,
                                             thin = 10 * itt,
                                             burnin = 30000 * itt)
        
        save(model_nonphylo_abundance, file = here("models", sprintf("%s_%s_ViralLoad_nonphylo_abundance.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_nonphylo_abundance.Rdata", g, c)))}
      
      
    }
    
  }

}