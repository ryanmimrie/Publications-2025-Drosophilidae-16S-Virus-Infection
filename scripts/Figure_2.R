# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure 2 =====================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script calculates and plots microbiome alpha diversities as
# violin plots across three experimental conditions. Underneath these plots,
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
#library(Cairo)

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

# ----- 1.3 Calculate Alpha Diversity Metrics ----------------------------------

data_diversity <- estimate_richness(data_abundances, measures = c("Observed", "Chao1", "Shannon", "Simpson"))

data_diversity$Pielou_evenness <- data_diversity$Shannon / log(data_diversity$Observed)

data_diversity$Simpson_evenness <- data_diversity$Simpson / max(data_diversity$Simpson)

data_diversity <- tibble::rownames_to_column(data_diversity, var = "SampleID")

data_diversity$ID <- str_remove_all(data_diversity$SampleID, "_R1.fastq.gz") %>% str_remove_all("X") %>% as.numeric()

data_diversity$SampleID <- NULL

data_abundances <- as.data.frame(as.matrix(sample_data(data_abundances)))

data_abundances$ID <- as.numeric(data_abundances$ID)

data_abundances <- left_join(as.data.frame(data_abundances), as.data.frame(data_diversity))

# ----- 1.4 Add columns for phylogenetic MCMCglmms -----------------------------

data_abundances <- filter(data_abundances, !Spc.Name.Short %in% c("Dtak", "Zdav"))

data_abundances$animal <- data_abundances$Spc.Name %>% str_remove_all(" ")

data_abundances <- as.data.frame(data_abundances)

data_abundances$Condition <- factor(data_abundances$Condition, levels = c("None", "Ringers", "DCV"))

data_abundances <- left_join(data_abundances, data_fly)

# ----- 1.5 Drop tips for species not in dataset -------------------------------

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_abundances$animal))

plot(tree_drosophilidae)

# ------------------------------------------------------------------------------
# ----- 2. Alpha Diversity Plots -----------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Wrangle Plot Object ------------------------------------------------

data_abundances_alphas <- select(data_abundances, ID, Spc.Name, Condition, Observed, Pielou_evenness, Shannon)

data_abundances_alphas <- gather(data_abundances_alphas, key = "Metric", value = "Value", 4:6)

data_abundances_alphas$Metric <- ifelse(data_abundances_alphas$Metric == "Observed", "Richness",
                                      ifelse(data_abundances_alphas$Metric == "Pielou_evenness", "Evenness", "Shannon Diversity"))

# ----- 2.2 Define Theme -------------------------------------------------------

theme_violin <- theme(aspect.ratio = 1,
                      text = element_text(size = 8, color = "#2e3440"),
                      axis.title = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid = element_blank(),
                      legend.title = element_text(size = 7),
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(0.3, "cm"),
                      strip.background = element_rect(fill = "#e5e9f0"))

# ----- 2.3 Richness Plot ------------------------------------------------------

# Significance points are added here, as each plot has a different y axis scale
# it is much easier to then duplicate these to the evenness and shannon plots
# manually in the svg file.

p_richness <- ggplot(filter(data_abundances_alphas, Metric == "Richness")) +
  geom_segment(aes(x = "None", xend = "Ringers", y = 34, y_end = 34)) +
  geom_segment(aes(x = "None", xend = "DCV", y = 32, y_end = 32)) +
  geom_segment(aes(x = "Ringers", xend = "DCV", y = 30, y_end = 30)) +
  geom_point(aes(x = "None",y = 34), size = 1) +
  geom_point(aes(x = "Ringers",y = 34), size = 1) +
  geom_point(aes(x = "None",y = 32), size = 1) +
  geom_point(aes(x = "DCV",y = 32), size = 1) +
  geom_point(aes(x = "Ringers",y = 30), size = 1) +
  geom_point(aes(x = "DCV",y = 30), size = 1) +
  geom_violin(aes(x = Condition, y = Value, fill = Condition, color = Condition), linewidth = 0.2, trim = F) +
  geom_boxplot(aes(x = Condition, y = Value), fill = "#34495e", color = "#34495e", width = 0.15, size = 0.5)+
  stat_summary(aes(x = Condition, y = Value), fun = median, geom = "point", color = "white", size = 1) +
  facet_wrap(~Metric, scales = "free") + 
  scale_fill_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#99D6FF", "#99FFC4", "#FFEB99")) +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#2980b9", "#2ecc71", "#f1c40f")) +
  scale_y_continuous(limits = c(0, 35)) +
  theme_bw() +
  theme_violin +
  theme(legend.position = "none") +
  ggtitle("A.")

# ----- 2.4 Evenness Plot ------------------------------------------------------

p_evenness <- ggplot(filter(data_abundances_alphas, Metric == "Evenness")) +
  geom_violin(aes(x = Condition, y = Value, fill = Condition, color = Condition), size = 0.2, trim = F) +
  geom_boxplot(aes(x = Condition, y = Value), fill = "#34495e", color = "#34495e", width = 0.15, size = 0.5)+
  stat_summary(aes(x = Condition, y = Value), fun = median, geom = "point", color = "white", size = 1) +
  facet_wrap(~Metric, scales = "free") + 
  scale_fill_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#99D6FF", "#99FFC4", "#FFEB99")) +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#2980b9", "#2ecc71", "#f1c40f")) +
  scale_y_continuous(limits = c(0, 1.35)) +
  theme_bw() +
  theme_violin +
  theme(legend.position = "none") +
  ggtitle("")

# ----- 2.5 Shannon Plot -------------------------------------------------------

p_shannon <- ggplot(filter(data_abundances_alphas, Metric == "Shannon Diversity")) +
  geom_violin(aes(x = Condition, y = Value, fill = Condition, color = Condition), size = 0.2, trim = F) +
  geom_boxplot(aes(x = Condition, y = Value), fill = "#34495e", color = "#34495e", width = 0.15, size = 0.5)+
  stat_summary(aes(x = Condition, y = Value), fun = median, geom = "point", color = "white", size = 1) +
  facet_wrap(~Metric, scales = "free") + 
  scale_fill_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#99D6FF", "#99FFC4", "#FFEB99")) +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#2980b9", "#2ecc71", "#f1c40f")) +
  scale_y_continuous(limits = c(0, 3.6)) +
  theme_bw() +
  theme_violin +
  ggtitle("")

# ------------------------------------------------------------------------------
# ----- 3. Alpha Diversity Models ----------------------------------------------
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

# ----- 3.6 Richness Heritability Model ----------------------------------------

if (!file.exists(here("models", "richness_heritability.Rdata"))) {
  
  model_richness_her <- MCMCglmm(Observed ~ Condition,
                                 random = ~ us(Condition):animal + us(Condition):Spc.Name + Diet + ID,
                                 rcov = ~idh(Condition):units,
                                 pedigree = tree_drosophilidae,
                                 prior = prior_her,
                                 data = data_abundances,
                                 nitt = 130000*itt,
                                 thin = 10*itt,
                                 burnin = 30000*itt,
                                 pr = TRUE)
  
  save(model_richness_her, file = here("models", "richness_heritability.Rdata"))

} else {load(here("models", "richness_heritability.Rdata"))}

richness_uninjured_her <- model_richness_her$VCV[,1] / (model_richness_her$VCV[,1] + model_richness_her$VCV[,10])
richness_saline_her <- model_richness_her$VCV[,5] / (model_richness_her$VCV[,5] + model_richness_her$VCV[,14])
richness_dcv_her <- model_richness_her$VCV[,9] / (model_richness_her$VCV[,9] + model_richness_her$VCV[,18])

# ----- 3.7 Richness Repeatability Model ---------------------------------------

if (!file.exists(here("models", "richness_repeatability.Rdata"))) {
  
  model_richness_rep <- MCMCglmm(Observed ~ Condition,
                                 random = ~ us(Condition):animal + Diet + ID,
                                 rcov = ~idh(Condition):units,
                                 pedigree = tree_drosophilidae,
                                 prior = prior_rep,
                                 data = data_abundances,
                                 nitt = 130000*itt,
                                 thin = 10*itt,
                                 burnin = 30000*itt,
                                 pr = TRUE)
  
  save(model_richness_rep, file = here("models", "richness_repeatability.Rdata"))
  
} else {load(here("models", "richness_repeatability.Rdata"))}

richness_uninjured_rep <- model_richness_rep$VCV[,1] / (model_richness_rep$VCV[,1] + model_richness_rep$VCV[,12])
richness_saline_rep <- model_richness_rep$VCV[,5] / (model_richness_rep$VCV[,5] + model_richness_rep$VCV[,13])
richness_dcv_rep <- model_richness_rep$VCV[,9] / (model_richness_rep$VCV[,9] + model_richness_rep$VCV[,14])

# ----- 3.8 Evenness Heritability Model ----------------------------------------

if (!file.exists(here("models", "evenness_heritability.Rdata"))) {
  
  model_evenness_her <- MCMCglmm(Pielou_evenness ~ Condition,
                                 random = ~ us(Condition):animal + us(Condition):Spc.Name + Diet + ID,
                                 rcov = ~idh(Condition):units,
                                 pedigree = tree_drosophilidae,
                                 prior = prior_her,
                                 data = data_abundances,
                                 nitt = 130000*itt,
                                 thin = 10*itt,
                                 burnin = 30000*itt,
                                 pr = TRUE)
  
  save(model_evenness_her, file = here("models", "evenness_heritability.Rdata"))
  
} else {load(here("models", "evenness_heritability.Rdata"))}

evenness_uninjured_her <- model_evenness_her$VCV[,1] / (model_evenness_her$VCV[,1] + model_evenness_her$VCV[,10])
evenness_saline_her <- model_evenness_her$VCV[,5] / (model_evenness_her$VCV[,5] + model_evenness_her$VCV[,14])
evenness_dcv_her <- model_evenness_her$VCV[,9] / (model_evenness_her$VCV[,9] + model_evenness_her$VCV[,18])

# ----- 3.9 Evenness Repeatability Model ---------------------------------------

if (!file.exists(here("models", "evenness_repeatability.Rdata"))) {
  
  model_evenness_rep <- MCMCglmm(Pielou_evenness ~ Condition,
                                 random = ~ us(Condition):animal + Diet + ID,
                                 rcov = ~idh(Condition):units,
                                 pedigree = tree_drosophilidae,
                                 prior = prior_rep,
                                 data = data_abundances,
                                 nitt = 130000*itt,
                                 thin = 10*itt,
                                 burnin = 30000*itt,
                                 pr = TRUE)
  
  save(model_evenness_rep, file = here("models", "evenness_repeatability.Rdata"))
  
} else {load(here("models", "evenness_repeatability.Rdata"))}

evenness_uninjured_rep <- model_evenness_rep$VCV[,1] / (model_evenness_rep$VCV[,1] + model_evenness_rep$VCV[,12])
evenness_saline_rep <- model_evenness_rep$VCV[,5] / (model_evenness_rep$VCV[,5] + model_evenness_rep$VCV[,13])
evenness_dcv_rep <- model_evenness_rep$VCV[,9] / (model_evenness_rep$VCV[,9] + model_evenness_rep$VCV[,14])

# ----- 3.10 Shannon Heritability Model ----------------------------------------

if (!file.exists(here("models", "shannon_heritability.Rdata"))) {
  
  model_shannon_her <- MCMCglmm(Shannon ~ Condition,
                                random = ~ us(Condition):animal + us(Condition):Spc.Name + Diet + ID,
                                rcov = ~idh(Condition):units,
                                pedigree = tree_drosophilidae,
                                prior = prior_her,
                                data = data_abundances,
                                nitt = 130000*itt,
                                thin = 10*itt,
                                burnin = 30000*itt,
                                pr = TRUE)
  
  save(model_shannon_her, file = here("models", "shannon_heritability.Rdata"))
  
} else {load(here("models", "shannon_heritability.Rdata"))}

shannon_uninjured_her <- model_shannon_her$VCV[,1] / (model_shannon_her$VCV[,1] + model_shannon_her$VCV[,10])
shannon_saline_her <- model_shannon_her$VCV[,5] / (model_shannon_her$VCV[,5] + model_shannon_her$VCV[,14])
shannon_dcv_her <- model_shannon_her$VCV[,9] / (model_shannon_her$VCV[,9] + model_shannon_her$VCV[,18])

# ----- 3.11 Shannon Repeatability Model ---------------------------------------

if (!file.exists(here("models", "shannon_repeatability.Rdata"))) {
  
  model_shannon_rep <- MCMCglmm(Shannon ~ Condition,
                                random = ~ us(Condition):animal + Diet + ID,
                                rcov = ~idh(Condition):units,
                                pedigree = tree_drosophilidae,
                                prior = prior_rep,
                                data = data_abundances,
                                nitt = 130000*itt,
                                thin = 10*itt,
                                burnin = 30000*itt,
                                pr = TRUE)
  
  save(model_shannon_rep, file = here("models", "shannon_repeatability.Rdata"))
  
} else {load(here("models", "shannon_repeatability.Rdata"))}

shannon_uninjured_rep <- model_shannon_rep$VCV[,1] / (model_shannon_rep$VCV[,1] + model_shannon_rep$VCV[,12])
shannon_saline_rep <- model_shannon_rep$VCV[,5] / (model_shannon_rep$VCV[,5] + model_shannon_rep$VCV[,13])
shannon_dcv_rep <- model_shannon_rep$VCV[,9] / (model_shannon_rep$VCV[,9] + model_shannon_rep$VCV[,14])

# ----- 3.12 Summarise Models --------------------------------------------------

quant_alphas <- bind_rows(
  summarise_MCMC("Heritability", "Richness", list(richness_uninjured_her, richness_saline_her, richness_dcv_her)),
  summarise_MCMC("Repeatability", "Richness", list(richness_uninjured_rep, richness_saline_rep, richness_dcv_rep)),
  summarise_MCMC("Heritability", "Evenness", list(evenness_uninjured_her, evenness_saline_her, evenness_dcv_her)),
  summarise_MCMC("Repeatability", "Evenness", list(evenness_uninjured_rep, evenness_saline_rep, evenness_dcv_rep)),
  summarise_MCMC("Heritability", "Shannon Diversity", list(shannon_uninjured_her, shannon_saline_her, shannon_dcv_her)),
  summarise_MCMC("Repeatability", "Shannon Diversity", list(shannon_uninjured_rep, shannon_saline_rep, shannon_dcv_rep))
)

quant_alphas$Metric <- factor(quant_alphas$Metric, levels = c("Richness", "Evenness", "Shannon Diversity"))

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

# ----- 4.2 Richness Model Plot ------------------------------------------------

p_model_richness <- ggplot(filter(quant_alphas, Metric == "Richness")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("B.")

# ----- 4.3 Evenness Model Plot ------------------------------------------------

p_model_evenness <- ggplot(filter(quant_alphas, Metric == "Evenness")) +
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

# ----- 4.4 Shannon Model Plot -------------------------------------------------

p_model_shannon <- ggplot(filter(quant_alphas, Metric == "Shannon Diversity")) +
  geom_point(aes(x = Mean, y = Condition, color = p_value < 0.05), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition, color = p_value < 0.05), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted",  size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_color_manual(name = "MCMC\np-value", labels = c("ns", "<0.05"), values = c("#bdc3c7", "#2c3e50")) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")

# ------------------------------------------------------------------------------
# ----- 5. Figure 2 Output -----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1 Arrange Figure -----------------------------------------------------

# Arranging in columns gives best alignment (still not perfect)

p_2 <- (p_richness / p_model_richness) | (p_evenness / p_model_evenness) | (p_shannon / p_model_shannon)

p_2

# ----- 5.2 Save Raw Outputs ---------------------------------------------------

ggsave(here("figures", "Figure 2 raw.png"), plot = p_2, dpi = 300, width = 7.5, height = 4)
ggsave(here("figures", "Figure 2 raw.svg"), plot = p_2, dpi = 300, width = 7.5, height = 4, device = cairo_pdf)

