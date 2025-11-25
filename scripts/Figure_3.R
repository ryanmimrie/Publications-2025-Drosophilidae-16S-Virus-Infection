# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure 3 =====================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs and plots a PCA on microbiome genera relative
# abundances across Drosophilidae host species. Underneath these plots,
# estimates of phylogenetic heritability and repeatability are calculated from
# MCMCglmm models and plotted as mean and HPD interval, with an arbitrary MCMCp
# measure of significance against a threshold of 0.05 used as a test for non-
# zero estimates.

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

data_clr <- as(otu_table(data_abundances), "matrix") %>%
  cmultRepl(method = "CZM", output = "p-counts", z.warning = 1, z.delete = F)

data_clr <- clr(data_clr) %>% t()

# ------------------------------------------------------------------------------
# ----- 2. PCA Analysis & Plots ------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. PCA ---------------------------------------------------------------

pca_results <- prcomp(data_clr, center = TRUE, scale. = FALSE)

summary(pca_results)

# ----- 2.2. Scree Plot --------------------------------------------------------

data_scree <- data.frame(PC = seq_along(pca_results$sdev), 
                         Variance = (pca_results$sdev)^2 / sum((pca_results$sdev)^2))

ggplot(data_scree, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity") + 
  scale_x_continuous(name = "Principal Component") +
  scale_y_continuous(name = "Variance Explained") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())

# ----- 2.3. Join to Metadata --------------------------------------------------

pca_results <- as.data.frame(pca_results$x[,1:4])

pca_results$sample <- rownames(pca_results)

data_abundances <- as.data.frame(as.matrix(sample_data(data_abundances)))

data_abundances <- left_join(data_abundances, pca_results)

data_abundances$Condition <- factor(data_abundances$Condition, levels = c("None", "Ringers", "DCV"))

# ----- 2.4. Define Theme ------------------------------------------------------

theme_PCA <- theme(aspect.ratio = 1,
                   text = element_text(size = 8, color = "#2e3440"),
                   panel.grid = element_blank(),
                   legend.position = "none",
                   strip.background = element_rect(fill = "#e5e9f0"))

# ----- 2.5. PC1/PC2 Plot ------------------------------------------------------

p_PC12 <- ggplot(data_abundances) +
  geom_point(aes(x = PC1, y = PC2, color = Condition), size = 1) +
  stat_ellipse(aes(x = PC1, y = PC2, color = Condition), size = 0.75) +
  scale_x_continuous(name = "PC 1 [25.2%]") +
  scale_y_continuous(name = "PC 2 [8.9%]") +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#2980b9", "#2ecc71", "#f1c40f")) +
  theme_bw() +
  theme_PCA +
  theme(legend.position = "none") +
  ggtitle("A.")

p_PC34 <- ggplot(data_abundances) +
  geom_point(aes(x = PC3, y = PC4, color = Condition), size = 1) +
  stat_ellipse(aes(x = PC3, y = PC4, color = Condition), size = 0.75) +
  scale_x_continuous(name = "PC 3 [6.7%]") +
  scale_y_continuous(name = "PC 4 [4.5%]") +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#2980b9", "#2ecc71", "#f1c40f")) +
  theme_bw() +
  theme_PCA +
  ggtitle("")

pPC <- p_PC12 | p_PC34

# ------------------------------------------------------------------------------
# ----- 3. PCA Models ----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1 Thresholded MCMCp function for Heritability/Repeatability ----------

MCMCp <- function(posterior, threshold = 0.05){
  
  p_above <- mean(posterior > threshold)
  p_below <- mean(posterior <= threshold)
  
  p_value <- 2 * min(p_above, p_below)
  
  return(p_value)
  
}

# ----- 3.2 MCMC summary function for Heritability/Repeatability ---------------

summarise_MCMC <- function(quant, metric, data_list) {
  conditions <- c("Uninjured", "Saline Injection", "DCV Infection")
  map_dfr(seq_along(conditions), function(i) {
    data.frame(
      Metric = metric,
      Condition = conditions[i],
      Quant = quant,
      Mean = mean(data_list[[i]]),
      HPD_Low = HPDinterval(data_list[[i]])[1],
      HPD_High = HPDinterval(data_list[[i]])[2],
      p_value = MCMCp(data_list[[i]])
    )
  })
}

# ----- 3.3 Heritability Prior -------------------------------------------------

prior_her <- list(G = list(
  G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G2 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G3 = list(V = diag(1), nu = 0.002),
  G3 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(3), nu = 0.002))

# ----- 3.4 Repeatability Prior ------------------------------------------------

prior_rep <- list(G = list(
  G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G2 = list(V = diag(1), nu = 0.002),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(3), nu = 0.002))

# ----- 3.5 Iteration Multiplier -----------------------------------------------

itt <- 10 # Set to 1 for demonstration, 100 for publication

# ----- 3.6 Add columns for phylogenetic MCMCglmms -----------------------------

data_abundances$animal <- data_abundances$Spc.Name %>% str_remove_all(" ")

data_abundances <- as.data.frame(data_abundances)

data_abundances <- left_join(data_abundances, data_fly)

# ----- 3.7 Drop tips for species not in dataset -------------------------------

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_abundances$animal))

plot(tree_drosophilidae)

# ----- 3.8 PC1 Heritability Model ---------------------------------------------

if (!file.exists(here("models", "PC1_heritability.Rdata"))) {
  
  model_PC1_her <- MCMCglmm(PC1 ~ Condition,
                            random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_her,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC1_her, file = here("models", "PC1_heritability.Rdata"))

} else {load(here("models", "PC1_heritability.Rdata"))}

PC1_uninjured_her <- model_PC1_her$VCV[,1] / (model_PC1_her$VCV[,1] + model_PC1_her$VCV[,4])
PC1_saline_her <- model_PC1_her$VCV[,2] / (model_PC1_her$VCV[,2] + model_PC1_her$VCV[,5])
PC1_dcv_her <- model_PC1_her$VCV[,3] / (model_PC1_her$VCV[,3] + model_PC1_her$VCV[,6])

# ----- 3.9 PC1 Repeatability Model --------------------------------------------

if (!file.exists(here("models", "PC1_repeatability.Rdata"))) {
  
  model_PC1_rep <- MCMCglmm(PC1 ~ Condition,
                            random = ~ idh(Condition):animal + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_rep,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC1_rep, file = here("models", "PC1_repeatability.Rdata"))

} else {load(here("models", "PC1_repeatability.Rdata"))}

PC1_uninjured_rep <- model_PC1_rep$VCV[,1] / (model_PC1_rep$VCV[,1] + model_PC1_rep$VCV[,5])
PC1_saline_rep <- model_PC1_rep$VCV[,2] / (model_PC1_rep$VCV[,2] + model_PC1_rep$VCV[,6])
PC1_dcv_rep <- model_PC1_rep$VCV[,3] / (model_PC1_rep$VCV[,3] + model_PC1_rep$VCV[,7])

# ----- 3.10 PC2 Heritability Model --------------------------------------------

if (!file.exists(here("models", "PC2_heritability.Rdata"))) {

  model_PC2_her <- MCMCglmm(PC2 ~ Condition,
                            random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_her,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC2_her, file = here("models", "PC2_heritability.Rdata"))

} else {load(here("models", "PC2_heritability.Rdata"))}

PC2_uninjured_her <- model_PC2_her$VCV[,1] / (model_PC2_her$VCV[,1] + model_PC2_her$VCV[,4])
PC2_saline_her <- model_PC2_her$VCV[,2] / (model_PC2_her$VCV[,2] + model_PC2_her$VCV[,5])
PC2_dcv_her <- model_PC2_her$VCV[,3] / (model_PC2_her$VCV[,3] + model_PC2_her$VCV[6])



# ----- 3.11 PC2 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "PC2_repeatability.Rdata"))) {

  model_PC2_rep <- MCMCglmm(PC2 ~ Condition,
                            random = ~ idh(Condition):animal + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_rep,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC2_rep, file = here("models", "PC2_repeatability.Rdata"))

} else {load(here("models", "PC2_repeatability.Rdata"))}

PC2_uninjured_rep <- model_PC2_rep$VCV[,1] / (model_PC2_rep$VCV[,1] + model_PC2_rep$VCV[,5])
PC2_saline_rep <- model_PC2_rep$VCV[,2] / (model_PC2_rep$VCV[,2] + model_PC2_rep$VCV[,6])
PC2_dcv_rep <- model_PC2_rep$VCV[,3] / (model_PC2_rep$VCV[,3] + model_PC2_rep$VCV[,7])

# ----- 3.12 PC3 Heritability Model --------------------------------------------

if (!file.exists(here("models", "PC3_heritability.Rdata"))) {
  
  model_PC3_her <- MCMCglmm(PC3 ~ Condition,
                            random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_her,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC3_her, file = here("models", "PC3_heritability.Rdata"))

} else {load(here("models", "PC3_heritability.Rdata"))}

PC3_uninjured_her <- model_PC3_her$VCV[,1] / (model_PC3_her$VCV[,1] + model_PC3_her$VCV[,4])
PC3_saline_her <- model_PC3_her$VCV[,2] / (model_PC3_her$VCV[,2] + model_PC3_her$VCV[,5])
PC3_dcv_her <- model_PC3_her$VCV[,3] / (model_PC3_her$VCV[,3] + model_PC3_her$VCV[6])

# ----- 3.13 PC3 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "PC3_repeatability.Rdata"))) {
  
  model_PC3_rep <- MCMCglmm(PC3 ~ Condition,
                            random = ~ idh(Condition):animal + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_rep,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC3_rep, file = here("models", "PC3_repeatability.Rdata"))

} else {load(here("models", "PC3_repeatability.Rdata"))}

PC3_uninjured_rep <- model_PC3_rep$VCV[,1] / (model_PC3_rep$VCV[,1] + model_PC3_rep$VCV[,5])
PC3_saline_rep <- model_PC3_rep$VCV[,2] / (model_PC3_rep$VCV[,2] + model_PC3_rep$VCV[,6])
PC3_dcv_rep <- model_PC3_rep$VCV[,3] / (model_PC3_rep$VCV[,3] + model_PC3_rep$VCV[,7])

# ----- 3.14 PC4 Heritability Model --------------------------------------------

if (!file.exists(here("models", "PC4_heritability.Rdata"))) {
  
  model_PC4_her <- MCMCglmm(PC4 ~ Condition,
                            random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_her,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC4_her, file = here("models", "PC4_heritability.Rdata"))

} else {load(here("models", "PC4_heritability.Rdata"))}

PC4_uninjured_her <- model_PC4_her$VCV[,1] / (model_PC4_her$VCV[,1] + model_PC4_her$VCV[,4])
PC4_saline_her <- model_PC4_her$VCV[,2] / (model_PC4_her$VCV[,2] + model_PC4_her$VCV[,5])
PC4_dcv_her <- model_PC4_her$VCV[,3] / (model_PC4_her$VCV[,3] + model_PC4_her$VCV[,6])

# ----- 3.15 PC4 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "PC4_repeatability.Rdata"))) {

  model_PC4_rep <- MCMCglmm(PC4 ~ Condition,
                            random = ~ idh(Condition):animal + Diet + ID,
                            rcov = ~idh(Condition):units,
                            pedigree = tree_drosophilidae,
                            prior = prior_rep,
                            data = data_abundances,
                            nitt = 130000*itt,
                            thin = 10*itt,
                            burnin = 30000*itt,
                            pr = TRUE)
  
  save(model_PC4_rep, file = here("models", "PC4_repeatability.Rdata"))

} else {load(here("models", "PC4_repeatability.Rdata"))}

PC4_uninjured_rep <- model_PC4_rep$VCV[,1] / (model_PC4_rep$VCV[,1] + model_PC4_rep$VCV[,5])
PC4_saline_rep <- model_PC4_rep$VCV[,2] / (model_PC4_rep$VCV[,2] + model_PC4_rep$VCV[,6])
PC4_dcv_rep <- model_PC4_rep$VCV[,3] / (model_PC4_rep$VCV[,3] + model_PC4_rep$VCV[,7])

# ----- 3.16 Summarise Models --------------------------------------------------

quant_PCAs <- bind_rows(
  summarise_MCMC("Heritability", "PC 1 [25.2%]", list(PC1_uninjured_her, PC1_saline_her, PC1_dcv_her)),
  summarise_MCMC("Repeatability", "PC 1 [25.2%]", list(PC1_uninjured_rep, PC1_saline_rep, PC1_dcv_rep)),
  summarise_MCMC("Heritability", "PC 2 [8.9%]", list(PC2_uninjured_her, PC2_saline_her, PC2_dcv_her)),
  summarise_MCMC("Repeatability", "PC 2 [8.9%]", list(PC2_uninjured_rep, PC2_saline_rep, PC2_dcv_rep)),
  summarise_MCMC("Heritability", "PC 3 [6.7%]", list(PC3_uninjured_her, PC3_saline_her, PC3_dcv_her)),
  summarise_MCMC("Repeatability", "PC 3 [6.7%]", list(PC3_uninjured_rep, PC3_saline_rep, PC3_dcv_rep)),
  summarise_MCMC("Heritability", "PC 4 [4.5%]", list(PC4_uninjured_her, PC4_saline_her, PC4_dcv_her)),
  summarise_MCMC("Repeatability", "PC 4 [4.5%]", list(PC4_uninjured_rep, PC4_saline_rep, PC4_dcv_rep))
)

# ------------------------------------------------------------------------------
# ----- 4. Heritability/Repeatability Plots ------------------------------------
# ------------------------------------------------------------------------------

# ----- 4.1 Theme --------------------------------------------------------------

theme_herrep <- theme(aspect.ratio = 0.5,
                      text = element_text(size = 8, color = "#2e3440"),
                      axis.title = element_blank(),
                      panel.grid = element_blank(),
                      legend.title = element_text(size = 7),
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(0.3, "cm"),
                      legend.position = "none",
                      strip.background = element_rect(fill = "#e5e9f0"))

# ----- 4.2 PC1 Plot -----------------------------------------------------------

p_model_PC1 <- ggplot(filter(quant_PCAs, Metric == "PC 1 [25.2%]")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("B.")

# ----- 4.3 PC2 Plot -----------------------------------------------------------

p_model_PC2 <- ggplot(filter(quant_PCAs, Metric == "PC 2 [8.9%]")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")

# ----- 4.4 PC3 Plot -----------------------------------------------------------

p_model_PC3 <- ggplot(filter(quant_PCAs, Metric == "PC 3 [6.7%]")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")

# ----- 4.5 PC4 Plot -----------------------------------------------------------

p_model_PC4 <- ggplot(filter(quant_PCAs, Metric == "PC 4 [4.5%]")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")


# ------------------------------------------------------------------------------
# ----- 5. Figure 3 Output -----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1 Arrange Figure -----------------------------------------------------

p_3A <- p_PC12 | p_PC34

p_3B <- p_model_PC1 | p_model_PC2 | p_model_PC3 | p_model_PC4

p_3B

# ----- 5.2 Save Raw Outputs ---------------------------------------------------

ggsave(here("figures", "Figure 3A raw.svg"), plot = p_3A, dpi = 300, width = 7, height = 3.5, device = cairo_pdf)
ggsave(here("figures", "Figure 3A raw.png"), plot = p_3A, dpi = 300, width = 7, height = 3.5)

ggsave(here("figures", "Figure 3B raw.svg"), plot = p_3B, dpi = 300, width = 7, height = 2.5, device = cairo_pdf)
ggsave(here("figures", "Figure 3B raw.png"), plot = p_3B, dpi = 300, width = 7, height = 2.5)

