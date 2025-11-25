# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure 7 =====================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs and plots correlations in relative abundances
# between conditions across Drosophilidae host species.

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
library(Cairo)
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

data_clr <- filter(data_clr, presence == T)

data_clr <- left_join(data_clr, data_fly)

# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ----- 1.3 Wrangle To Wide Format for MCMCglmm --------------------------------

data_clr_wide <- select(data_clr, Block, animal, Genus, Condition, Diet, abundance)

data_clr_wide <- spread(data_clr_wide, key = Condition, value = abundance)

data_clr_wide_sum <- data_clr_wide %>% group_by(animal, Genus, Diet) %>%
  summarise(ab_none = mean(None, na.rm = T),
            ab_ring = mean(Ringers, na.rm = T),
            ab_dcv = mean(DCV, na.rm = T))

data_clr_wide_sum$spc <- data_clr_wide_sum$animal

# ------------------------------------------------------------------------------
# ----- 2. MCMCglmms -----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Iteration Multiplier ----------------------------------------------

itt <- 1 # Set to 1 for demonstration, 100 for publication

# ----- 2.2. Estimate Correlations ---------------------------------------------

prior <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(2), nu = 0.002))

if (!file.exists(here("models", "correlation_none_ringers.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_dcv <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_None_Ringers <- MCMCglmm(cbind(ab_none, ab_ring) ~ trait - 1,
                                 random = ~ us(trait):animal + us(trait):spc + Diet,
                                 rcov = ~ us(trait):units,
                                 pedigree = tree_drosophilidae,
                                 prior = prior,
                                 family = c("gaussian", "gaussian"),
                                 data = data,
                                 nitt = 130000*itt,
                                 thin = 10*itt,
                                 burnin = 30000*itt)
  
  save(model_None_Ringers, file = here("models", "correlation_none_ringers.Rdata"))
  
} else {
  load(here("models", "correlation_none_ringers.Rdata"))
}

if (!file.exists(here("models", "correlation_none_dcv.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_ring <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_None_DCV <- MCMCglmm(cbind(ab_dcv, ab_none) ~ trait - 1,
                             random = ~ us(trait):animal + us(trait):spc + Diet,
                             rcov = ~ us(trait):units,
                             pedigree = tree_drosophilidae,
                             prior = prior,
                             family = c("gaussian", "gaussian"),
                             data = data,
                             nitt = 130000*itt,
                             thin = 10*itt,
                             burnin = 30000*itt)
  
  save(model_None_DCV, file = here("models", "correlation_none_dcv.Rdata"))
  
} else {
  load(here("models", "correlation_none_dcv.Rdata"))
}


if (!file.exists(here("models", "correlation_ringers_dcv.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_none <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_Ringers_DCV <- MCMCglmm(cbind(ab_dcv, ab_ring) ~ trait - 1,
                                random = ~ us(trait):animal + us(trait):spc + Diet,
                                rcov = ~ us(trait):units,
                                pedigree = tree_drosophilidae,
                                prior = prior,
                                family = c("gaussian", "gaussian"),
                                data = data,
                                nitt = 130000*itt,
                                thin = 10*itt,
                                burnin = 30000*itt)
  
  save(model_Ringers_DCV, file = here("models", "correlation_ringers_dcv.Rdata"))
  
} else {
  load(here("models", "correlation_ringers_dcv.Rdata"))
}

cov_None_Ringers <- model_None_Ringers$VCV[,2] + model_None_Ringers$VCV[,6] + model_None_Ringers$VCV[,11]
varx_None_Ringers <- model_None_Ringers$VCV[,1] + model_None_Ringers$VCV[,5] + model_None_Ringers$VCV[,10]
vary_None_Ringers <- model_None_Ringers$VCV[,4] + model_None_Ringers$VCV[,8] + model_None_Ringers$VCV[,13]

cor_None_Ringers <- cov_None_Ringers / sqrt(varx_None_Ringers * vary_None_Ringers)
slope_None_Ringers <- cov_None_Ringers / varx_None_Ringers

cov_None_DCV <- model_None_DCV$VCV[,2] + model_None_DCV$VCV[,6] + model_None_DCV$VCV[,11]
varx_None_DCV <- model_None_DCV$VCV[,1] + model_None_DCV$VCV[,5] + model_None_DCV$VCV[,10]
vary_None_DCV <- model_None_DCV$VCV[,4] + model_None_DCV$VCV[,8] + model_None_DCV$VCV[,13]

cor_None_DCV <- cov_None_DCV / sqrt(varx_None_DCV * vary_None_DCV)
slope_None_DCV <- cov_None_DCV / varx_None_DCV

cov_Ringers_DCV <- model_Ringers_DCV$VCV[,2] + model_Ringers_DCV$VCV[,6] + model_Ringers_DCV$VCV[,11]
varx_Ringers_DCV <- model_Ringers_DCV$VCV[,1] + model_Ringers_DCV$VCV[,5] + model_Ringers_DCV$VCV[,10]
vary_Ringers_DCV <- model_Ringers_DCV$VCV[,4] + model_Ringers_DCV$VCV[,8] + model_Ringers_DCV$VCV[,13]

cor_Ringers_DCV <- cov_Ringers_DCV / sqrt(varx_Ringers_DCV * vary_Ringers_DCV)
slope_Ringers_DCV <- cov_Ringers_DCV / varx_Ringers_DCV

p1 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_none, y = ab_ring), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(slope_None_Ringers)) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(No Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(Saline Injection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))

p2 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_none, y = ab_dcv), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(slope_None_DCV)) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(No Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(DCV Infection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))


p3 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_ring, y = ab_dcv), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(slope_Ringers_DCV)) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(Saline Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(DCV Infection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))


plots <- p1 | p2 | p3

ggsave(here("figures", "Figure 7 raw.svg"), plot = plots, dpi = 300, width = 7, height = 2.5, device = cairo_pdf)
ggsave(here("figures", "Figure 7 raw.png"), plot = plots, dpi = 300, width = 7, height = 2.5)










