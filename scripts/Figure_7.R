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

# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ----- 1.3 Wrangle To Wide Format for MCMCglmm --------------------------------

data_clr_wide <- select(data_clr, Block, animal, Genus, Condition, abundance)

data_clr_wide <- spread(data_clr_wide, key = Condition, value = abundance)

data_clr_wide_sum <- data_clr_wide %>% group_by(animal, Genus) %>%
  summarise(ab_none = mean(None, na.rm = T),
            ab_ring = mean(Ringers, na.rm = T),
            ab_dcv = mean(DCV, na.rm = T))

# ------------------------------------------------------------------------------
# ----- 2. MCMCglmms -----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Iteration Multiplier ----------------------------------------------

itt <- 100 # Set to 1 for demonstration, 100 for publication

# ----- 2.2. Check for Phylogenetic Correlation --------------------------------

prior <- list(G = list(
  G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(3), nu = 0.002))

if (!file.exists(here("models", "correlation_phylo.Rdata"))) {

  model_phylo <- MCMCglmm(abundance ~ Condition,
                          random = ~ us(Condition):animal + ID,
                          rcov = ~idh(Condition):units,
                          pedigree = tree_drosophilidae,
                          prior = prior,
                          data = data_clr,
                          nitt = 130000*itt,
                          thin = 10*itt,
                          burnin = 30000*itt,
                          pr = TRUE)

  save(model_phylo, file = here("models", "correlation_phylo.Rdata"))
  
} else {load(here("models", "correlation_phylo.Rdata"))}

summary(model_phylo)

# ----- 2.3. Estimate Correlation Slope without Phylogeny ----------------------

prior <- list(R = list(V = diag(1), nu = 0.002))

if (!file.exists(here("models", "correlation_none_ringers_slope.Rdata"))) {
  
  data_curr <- data_clr_wide_sum
  data_curr$ab_dcv <- NULL
  data_curr <- as.data.frame(data_curr)
  
  model_None_Ringers_slope <- MCMCglmm(ab_none ~ 0 + ab_ring,
                                       rcov = ~units,
                                       prior = prior,
                                       data = na.omit(data_curr),
                                       nitt = 130000*itt,
                                       thin = 10*itt,
                                       burnin = 30000*itt,
                                       pr = TRUE)
  
  save(model_None_Ringers_slope, file = here("models", "correlation_none_ringers_slope.Rdata"))
  
} else {load(here("models", "correlation_none_ringers_slope.Rdata"))}

if (!file.exists(here("models", "correlation_none_dcv_slope.Rdata"))) {

  data_curr <- data_clr_wide_sum
  data_curr$ab_ring <- NULL
  data_curr <- as.data.frame(data_curr)
  
  model_None_DCV_slope <- MCMCglmm(ab_none ~ 0 + ab_dcv,
                                   rcov = ~units,
                                   prior = prior,
                                   data = na.omit(data_curr),
                                   nitt = 130000*itt,
                                   thin = 10*itt,
                                   burnin = 30000*itt,
                                   pr = TRUE)

  save(model_None_DCV_slope, file = here("models", "correlation_none_dcv_slope.Rdata"))
  
} else {load(here("models", "correlation_none_dcv_slope.Rdata"))}
  
if (!file.exists(here("models", "correlation_ringers_dcv_slope.Rdata"))) {

  data_curr <- data_clr_wide_sum
  data_curr$ab_none <- NULL
  data_curr <- as.data.frame(data_curr)
  
  model_Ringers_DCV_slope <- MCMCglmm(ab_ring ~ 0 + ab_dcv,
                                      rcov = ~units,
                                      prior = prior,
                                      data = na.omit(data_curr),
                                      nitt = 130000*itt,
                                      thin = 10*itt,
                                      burnin = 30000*itt,
                                      pr = TRUE)
  
  save(model_Ringers_DCV_slope, file = here("models", "correlation_ringers_dcv_slope.Rdata"))

} else {load(here("models", "correlation_ringers_dcv_slope.Rdata"))}

# ----- 2.4. Estimate Correlation Coefficient without Phylogeny ----------------

prior <- list(R = list(V = diag(2), nu = 0.002))

if (!file.exists(here("models", "correlation_none_ringers_coeff.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_dcv <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_None_Ringers_coeff <- MCMCglmm(cbind(ab_none, ab_ring) ~ trait - 1,
                                       rcov = ~ us(trait):units,
                                       prior = list(R = list(V = diag(2), nu = 2)),
                                       family = c("gaussian", "gaussian"),
                                       data = data,
                                       nitt = 130000*itt,
                                       thin = 10*itt,
                                       burnin = 30000*itt)
  
  save(model_None_Ringers_coeff, file = here("models", "correlation_none_ringers_coeff.Rdata"))
  
} else {
  load(here("models", "correlation_none_ringers_coeff.Rdata"))
}

if (!file.exists(here("models", "correlation_none_dcv_coeff.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_ring <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_None_DCV_coeff <- MCMCglmm(cbind(ab_none, ab_dcv) ~ trait - 1,
                                   rcov = ~ us(trait):units,
                                   prior = list(R = list(V = diag(2), nu = 2)),
                                   family = c("gaussian", "gaussian"),
                                   data = data,
                                   nitt = 130000*itt,
                                   thin = 10*itt,
                                   burnin = 30000*itt)
  
  save(model_None_DCV_coeff, file = here("models", "correlation_none_dcv_coeff.Rdata"))
  
} else {
  load(here("models", "correlation_none_dcv_coeff.Rdata"))
}


if (!file.exists(here("models", "correlation_ringers_dcv_coeff.Rdata"))) {
  
  data <- data_clr_wide_sum
  data$ab_none <- NULL
  data <- na.omit(data)
  data <- as.data.frame(data)
  
  model_Ringers_DCV_coeff <- MCMCglmm(cbind(ab_ring, ab_dcv) ~ trait - 1,
                                      rcov = ~ us(trait):units,
                                      prior = list(R = list(V = diag(2), nu = 2)),
                                      family = c("gaussian", "gaussian"),
                                      data = data,
                                      nitt = 130000*itt,
                                      thin = 10*itt,
                                      burnin = 30000*itt)
  
  save(model_Ringers_DCV_coeff, file = here("models", "correlation_ringers_dcv_coeff.Rdata"))
  
} else {
  load(here("models", "correlation_ringers_dcv_coeff.Rdata"))
}


p1 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_none, y = ab_ring), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(model_None_Ringers_slope$Sol[,1])) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(No Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(Saline Injection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))

p2 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_none, y = ab_dcv), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(model_None_DCV_slope$Sol[,1])) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(No Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(DCV Infection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))


p3 <- ggplot(data_clr_wide_sum) +
  geom_point(aes(x = ab_ring, y = ab_dcv), color = "#2c3e50", alpha = 1, size = 0.75) +
  geom_abline(slope = mean(model_Ringers_DCV_slope$Sol[,1])) +
  scale_x_continuous(name = "CLR-transformed Abundance\n(Saline Injection)", limits = c(-1.25, 10), expand = c(0,0)) +
  scale_y_continuous(name = "CLR-transformed Abundance\n(DCV Infection)",limits = c(-1.25, 10), expand = c(0,0)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(size = 8, color = "#2e3440"))


plots <- p1 | p2 | p3

ggsave(here("figures", "Figure 7 raw.svg"), plot = plots, dpi = 300, width = 7, height = 2.5, device = cairo_pdf)
ggsave(here("figures", "Figure 7 raw.png"), plot = plots, dpi = 300, width = 7, height = 2.5)




dimnames(model_None_Ringers_coeff$VCV)

mean(model_Ringers_DCV_slope$Sol[,1])
HPDinterval(model_Ringers_DCV_slope$Sol[,1])

mean(model_Ringers_DCV_coeff$VCV[,2] /
       sqrt(model_Ringers_DCV_coeff$VCV[,1] * model_Ringers_DCV_coeff$VCV[,4]))

HPDinterval(model_Ringers_DCV_coeff$VCV[,2] /
       sqrt(model_Ringers_DCV_coeff$VCV[,1] * model_Ringers_DCV_coeff$VCV[,4]))


