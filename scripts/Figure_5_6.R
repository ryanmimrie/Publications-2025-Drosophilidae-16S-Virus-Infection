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

# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ----- 1.3 Wrangle To Wide Format for MCMCglmm --------------------------------

data_clr_wide_presence <- select(data_clr, ID, Condition, animal, Genus, presence)
data_clr_wide_presence <- spread(data_clr_wide_presence, key = Genus, value = presence)

data_clr_wide_abundance <- select(data_clr, ID, Condition, animal, Genus, abundance)
data_clr_wide_abundance <- spread(data_clr_wide_abundance, key = Genus, value = abundance)

# ------------------------------------------------------------------------------
# ----- 2. Abundance Models ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Model Prior --------------------------------------------------------

prior_pr <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002, fix = 1))

prior_ab <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

# ----- 2.2 Iteration Multiplier -----------------------------------------------

itt <- 100 # Set to 1 for demonstration, 100 for publication

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

for(Genus in unique(data_clr$Genus)){
  
  fixedeffs <- as.formula(paste(Genus, "~ Condition"))
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("None", "Ringers", "DCV"))
  data_clr_wide_abundance$Condition <- factor(data_clr_wide_abundance$Condition, levels = c("None", "Ringers", "DCV"))
  
  data_clr_wide_abundance_notzero <- data_clr_wide_abundance[data_clr_wide_presence[[Genus]] == TRUE, ]
  

  if (!file.exists(here("models", sprintf("%s_NoneContrast_presence.Rdata", Genus)))) {
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal,
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
  
  Pr_none <- as.mcmc(model_pr$Sol[,1])
  Pr_ring <- as.mcmc(model_pr$Sol[,1] + model_pr$Sol[,2])
  Pr_dcv <- as.mcmc(model_pr$Sol[,1] + model_pr$Sol[,3])
  
  Contrast_pr_ring <- as.mcmc(model_pr$Sol[,2])
  
  if (!file.exists(here("models", sprintf("%s_NoneContrast_abundance.Rdata", Genus)))) {
  
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal,
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
    
  Ab_none <- as.mcmc(model_ab$Sol[,1])
  Ab_ring <- as.mcmc(model_ab$Sol[,1] + model_ab$Sol[,2])
  Ab_dcv <- as.mcmc(model_ab$Sol[,1] + model_ab$Sol[,3])
  
  Contrast_ab_ring <- as.mcmc(model_ab$Sol[,2])
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("Ringers", "None", "DCV"))
  data_clr_wide_abundance_notzero$Condition <- factor(data_clr_wide_abundance_notzero$Condition, levels = c("Ringers", "None", "DCV"))
  
  if (!file.exists(here("models", sprintf("%s_RingersContrast_presence.Rdata", Genus)))) {
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal,
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
  
  Contrast_pr_dcv <- as.mcmc(model_pr$Sol[,3])
  
  
  if (!file.exists(here("models", sprintf("%s_RingersContrast_abundance.Rdata", Genus)))) {
    
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal,
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
  
  Contrast_ab_dcv <- as.mcmc(model_ab$Sol[,3])
  
  Contrast_pr_ring_probscale <- pnorm(Pr_ring) - pnorm(Pr_none)
  Contrast_pr_dcv_probscale <- pnorm(Pr_dcv) - pnorm(Pr_ring)
  
  data_contrasts <- rbind(data_contrasts,
                          data.frame(Genus = Genus,
                                     Condition = c("None", "Ringers", "DCV"),
                                     Value_pr = c(mean(Pr_none), mean(Pr_ring), mean(Pr_dcv)),
                                     Value_pr_low = c(HPDinterval(Pr_none)[1], HPDinterval(Pr_ring)[1], HPDinterval(Pr_dcv)[1]),
                                     Value_pr_high = c(HPDinterval(Pr_none)[2], HPDinterval(Pr_ring)[2], HPDinterval(Pr_dcv)[2]),
                                     Value_ab = c(mean(Ab_none), mean(Ab_ring), mean(Ab_dcv)),
                                     Value_ab_low = c(HPDinterval(Ab_none)[1], HPDinterval(Ab_ring)[1], HPDinterval(Ab_dcv)[1]),
                                     Value_ab_high = c(HPDinterval(Ab_none)[2], HPDinterval(Ab_ring)[2], HPDinterval(Ab_dcv)[2]),
                                     Contrast_pr = c(0, mean(Contrast_pr_ring), mean(Contrast_pr_dcv)),
                                     Contrast_pr_low = c(0, HPDinterval(Contrast_pr_ring)[1], HPDinterval(Contrast_pr_dcv)[1]),
                                     Contrast_pr_high = c(0, HPDinterval(Contrast_pr_ring)[2], HPDinterval(Contrast_pr_dcv)[2]),
                                     Contrast_pr_prob = c(0, mean(Contrast_pr_ring_probscale), mean(Contrast_pr_dcv_probscale)),
                                     Contrast_pr_prob_low = c(0, HPDinterval(Contrast_pr_ring_probscale)[1], HPDinterval(Contrast_pr_dcv_probscale)[1]),
                                     Contrast_pr_prob_high = c(0, HPDinterval(Contrast_pr_ring_probscale)[2], HPDinterval(Contrast_pr_dcv_probscale)[2]),
                                     Contrast_ab = c(0, mean(Contrast_ab_ring), mean(Contrast_ab_dcv)),
                                     Contrast_ab_low = c(0, HPDinterval(Contrast_ab_ring)[1], HPDinterval(Contrast_ab_dcv)[1]),
                                     Contrast_ab_high = c(0, HPDinterval(Contrast_ab_ring)[2], HPDinterval(Contrast_ab_dcv)[2])))
  
}



# ------------------------------------------------------------------------------
# ----- 3. Figure 5 ------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1. Wrangle Plot Elements ---------------------------------------------

data_contrasts$Condition <- factor(data_contrasts$Condition,
                                   levels = c("None", "Ringers", "DCV"))

data_contrasts$Genus <- ifelse(data_contrasts$Genus == "Escherichia", "Escherichia-Shigella", data_contrasts$Genus)

bacteria <- c("Levilactobacillus", "Fructilactobacillus", "Companilactobacillus", "Lactobacillus",
              "Enterococcus", "Streptococcus", "Staphylococcus", "Bacillus",
              "Finegoldia", "Lawsonella", "Corynebacterium", "Rhodococcus",
              "Cutibacterium", "Afipia", "Bosea", "Brucella",
              "Brevundimonas", "Enhydrobacter", "Commensalibacter", "Acetobacter",
              "Gluconobacter", "Buttiauxella", "Providencia", "Serratia",
              "Escherichia-Shigella", "Pseudomonas", "Ralstonia", "Pelomonas",
              "Stenotrophomonas", "Dysgonomonas", "Bacteroides", "Sediminibacterium")


data_contrasts$Genus <- factor(data_contrasts$Genus, levels = rev(bacteria))



p1 <- ggplot(filter(data_contrasts, Condition == "None")) +
  geom_errorbar(aes(xmin = pnorm(Value_pr_low), xmax = pnorm(Value_pr_high), y = Genus), width = 0, color = "#2c3e50") +
  geom_point(aes(x = pnorm(Value_pr), y = Genus), color = "#2c3e50") +
  scale_x_continuous(name = "Detection Probability") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p2 <- ggplot(filter(data_contrasts, Condition == "Ringers")) +
  geom_errorbar(aes(xmin = Contrast_pr_prob_low, xmax = Contrast_pr_prob_high, y = Genus), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_pr_prob, y = Genus), color = "#bdc3c7") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in Detection Probability", limits = c(-0.4, 0.4), expand = c(0,0), breaks = c(-0.25, 0, 0.25)) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p3 <- ggplot(filter(data_contrasts, Condition == "DCV")) +
  geom_errorbar(aes(xmin = Contrast_pr_prob_low, xmax = Contrast_pr_prob_high, y = Genus), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_pr_prob, y = Genus), color = "#bdc3c7") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in Detection Probability", limits = c(-0.4, 0.4), expand = c(0,0), breaks = c(-0.25, 0, 0.25)) +
  scale_y_discrete(position = "right") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p4 <- ggplot(filter(data_contrasts, Condition == "None")) +
  geom_errorbar(aes(xmin = Value_ab_low, xmax = Value_ab_high, y = Genus), width = 0, color = "#2c3e50") +
  geom_point(aes(x = Value_ab, y = Genus), color = "#2c3e50") +
  scale_x_continuous(name = "CLR-transformed Abundance") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p5 <- ggplot(filter(data_contrasts, Condition == "Ringers")) +
  geom_errorbar(aes(xmin = Contrast_ab_low, xmax = Contrast_ab_high, y = Genus), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_ab, y = Genus), color = "#bdc3c7") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in CLR-transformed Abundance", limits = c(-4, 4), expand = c(0,0)) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p6 <- ggplot(filter(data_contrasts, Condition == "DCV")) +
  geom_errorbar(aes(xmin = Contrast_ab_low, xmax = Contrast_ab_high, y = Genus), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_ab, y = Genus), color = "#bdc3c7") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in CLR-transformed Abundance", limits = c(-4, 4), expand = c(0,0)) +
  scale_y_discrete(position = "right") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

plot <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave(here("figures", "Figure 5_6 raw.svg"), plot = plot, dpi = 300, width = 7, height = 9, device = cairo_pdf)
ggsave(here("figures", "Figure 5_6 raw.png"), plot = plot, dpi = 300, width = 7, height = 9)

