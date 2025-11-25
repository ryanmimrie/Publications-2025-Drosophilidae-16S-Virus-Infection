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
library(Cairo)
library(zCompositions)
library(compositions)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ----- 0.3. Load Data ---------------------------------------------------------

load(here("data", "data_marginalseq_cutoff.RData"))

data_metadata <- read_csv(here("data", "data_sampleMetadata.csv"))

data_fly <- read_csv(here("data", "data_drosophilidae.csv"))

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

data_clr <- left_join(data_clr, data_fly)


# ----- 1.2 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ------------------------------------------------------------------------------
# ----- 2. Abundance Models ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Model Prior --------------------------------------------------------

prior_presence <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(2), nu = 0.002, fix = 2))

prior_abundance <- list(G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(2), nu = 0.002))

# ----- 2.2 Iteration Multiplier -----------------------------------------------

itt <- 1 # Set to 1 for demonstration, 100 for publication

# ----- 2.3 Model Loop ---------------------------------------------------------

data_correlations <- data.frame(Genus = NA,
                                Condition = NA,
                                B_pr_marginal = NA,
                                B_pr_marginal_low = NA,
                                B_pr_marginal_high = NA,
                                R_pr_marginal = NA,
                                R_pr_marginal_low = NA,
                                R_pr_marginal_high = NA,
                                B_ab_marginal = NA,
                                B_ab_marginal_low = NA,
                                B_ab_marginal_high = NA,
                                R_ab_marginal = NA,
                                R_ab_marginal_low = NA,
                                R_ab_marginal_high = NA,
                                C_ab_marginal = NA,
                                C_ab_marginal_low = NA,
                                C_ab_marginal_high = NA)[-1,]

data_clr <- filter(data_clr, Genus %in% unique(data_clr$Genus))

for(g in unique(data_clr$Genus)){
  for(c in c("None", "DCV")){
    
    data <- filter(data_clr, Genus == g, Condition == c)
    data <- select(data, abundance, presence, foldchange, animal, Diet, ID)
    
    data_pres <- select(data, presence, foldchange, animal, ID, Diet)
    data_abun <- filter(data, presence) %>% select(abundance, foldchange, animal, ID, Diet)
    
    data_pres <- gather(data_pres, key = "trait", value = "value", 1:2)
    data_abun <- gather(data_abun, key = "trait", value = "value", 1:2)
    
    data_pres$trait <- factor(data_pres$trait, levels = c("foldchange", "presence"))
    data_abun$trait <- factor(data_abun$trait, levels = c("foldchange", "abundance"))
    
    data_pres <- as.data.frame(data_pres)
    data_abun <- as.data.frame(data_abun)
    
    ab <- filter(data_abun, trait == "abundance")
    
    data_pres$family <- NULL
    data_abun$family <- NULL
    
    data_pres <- spread(data_pres, key = "trait", value = "value")
    data_abun <- spread(data_abun, key = "trait", value = "value")
    
    data_pres$spc <- data_pres$animal
    data_abun$spc <- data_abun$animal
    
    if (length(unique(ab$animal))>=1){
      
      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_presence.Rdata", g, c)))) {
        
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_presence.Rdata", g, c)))
        
        model_presence <- MCMCglmm(cbind(foldchange, presence) ~ trait - 1,
                                   random = ~ us(trait):animal + us(trait):spc + Diet,
                                   rcov = ~ us(trait):units,
                                   family = c("gaussian", "threshold"),
                                   prior = prior_presence,
                                   data = data_pres,
                                   nitt = 130000 * itt,
                                   thin = 10 * itt,
                                   burnin = 30000 * itt)
        
        save(model_presence, file = here("models", sprintf("%s_%s_ViralLoad_presence.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_presence.Rdata", g, c)))}
      
      if (!file.exists(here("models", sprintf("%s_%s_ViralLoad_abundance.Rdata", g, c)))) {
        
        save(holder, file = here("models", sprintf("%s_%s_ViralLoad_abundance.Rdata", g, c)))
        
        model_abundance <- MCMCglmm(cbind(foldchange, abundance) ~ trait - 1,
                                    random = ~ us(trait):animal + us(trait):spc + Diet,
                                    rcov = ~ us(trait):units,
                                             family = c("gaussian", "gaussian"),
                                             prior = prior_abundance,
                                             data = data_abun,
                                             nitt = 130000 * itt,
                                             thin = 10 * itt,
                                             burnin = 30000 * itt)
        
        save(model_abundance, file = here("models", sprintf("%s_%s_ViralLoad_abundance.Rdata", g, c)))
        
      } else {load(here("models", sprintf("%s_%s_ViralLoad_abundance.Rdata", g, c)))}
      
      
      # Presence marginal
      var_pr_marginal_viral <- model_presence$VCV[,1] + model_presence$VCV[,5] + model_presence$VCV[,10]
      var_pr_marginal_pr <- model_presence$VCV[,4] + model_presence$VCV[,8] + model_presence$VCV[,13]
      cov_pr_marginal_pr <- model_presence$VCV[,2] + model_presence$VCV[,6] + model_presence$VCV[,11]
      B_pr_marginal <- as.mcmc(cov_pr_marginal_pr / var_pr_marginal_pr)
      R_pr_marginal <- as.mcmc(cov_pr_marginal_pr /(sqrt(var_pr_marginal_viral * var_pr_marginal_pr)))
      
      # Abundance marginal
      var_ab_marginal_viral <- model_abundance$VCV[,1] + model_presence$VCV[,5] + model_presence$VCV[,10]
      var_ab_marginal_ab <- model_abundance$VCV[,4] + model_presence$VCV[,8] + model_presence$VCV[,13]
      cov_ab_marginal_ab <- model_abundance$VCV[,2] + model_presence$VCV[,6] + model_presence$VCV[,11]
      ab_marginal_viral <- model_abundance$Sol[,1]
      ab_marginal_ab <- model_abundance$Sol[,2]
      B_ab_marginal <- as.mcmc(cov_ab_marginal_ab / var_ab_marginal_ab)
      R_ab_marginal <- as.mcmc(cov_ab_marginal_ab /(sqrt(var_ab_marginal_viral * var_ab_marginal_ab)))
      C_ab_marginal <- ab_marginal_viral - (B_ab_marginal * ab_marginal_ab)
      
      data_correlations <- rbind(data_correlations,
                                 data.frame(Genus = g,
                                            Condition = c,
                                            B_pr_marginal = mean(B_pr_marginal),
                                            B_pr_marginal_low = HPDinterval(B_pr_marginal)[1],
                                            B_pr_marginal_high = HPDinterval(B_pr_marginal)[2],
                                            R_pr_marginal = mean(R_pr_marginal),
                                            R_pr_marginal_low = HPDinterval(R_pr_marginal)[1],
                                            R_pr_marginal_high = HPDinterval(R_pr_marginal)[2],
                                            B_ab_marginal = mean(B_ab_marginal),
                                            B_ab_marginal_low = HPDinterval(B_ab_marginal)[1],
                                            B_ab_marginal_high = HPDinterval(B_ab_marginal)[2],
                                            R_ab_marginal = mean(R_ab_marginal),
                                            R_ab_marginal_low = HPDinterval(R_ab_marginal)[1],
                                            R_ab_marginal_high = HPDinterval(R_ab_marginal)[2],
                                            C_ab_marginal = mean(C_ab_marginal),
                                            C_ab_marginal_low = HPDinterval(C_ab_marginal)[1],
                                            C_ab_marginal_high = HPDinterval(C_ab_marginal)[2]))
      
    }
    
  }
  
}

# ------------------------------------------------------------------------------
# ----- 3. Figure 7 ------------------------------------------------------------
# ------------------------------------------------------------------------------

bacteria <- c("Levilactobacillus", "Fructilactobacillus", "Companilactobacillus", "Lactobacillus",
              "Enterococcus", "Streptococcus", "Staphylococcus", "Bacillus",
              "Finegoldia", "Lawsonella", "Corynebacterium", "Rhodococcus",
              "Cutibacterium", "Afipia", "Bosea", "Brucella",
              "Brevundimonas", "Enhydrobacter", "Commensalibacter", "Acetobacter",
              "Gluconobacter", "Buttiauxella", "Providencia", "Serratia",
              "Escherichia-Shigella", "Pseudomonas", "Ralstonia", "Pelomonas",
              "Stenotrophomonas", "Dysgonomonas", "Bacteroides", "Sediminibacterium")

data_correlations$Condition <- factor(data_correlations$Condition, levels = c("None", "DCV"))

data_correlations$Genus <- ifelse(data_correlations$Genus == "Escherichia", "Escherichia-Shigella", data_correlations$Genus)

data_correlations$Genus <- factor(data_correlations$Genus, levels = rev(bacteria))

p1 <- ggplot(data_correlations) +
  geom_errorbar(aes(xmin = R_pr_marginal_low, xmax = R_pr_marginal_high, y = Genus, color = sign(R_pr_marginal_low) == sign(R_pr_marginal_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_pr_marginal, y = Genus, color = sign(R_pr_marginal_low) == sign(R_pr_marginal_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Correlation Coefficient (Viral Load vs. Bacterial Presence)") +
  scale_color_manual(values = c("lightgrey", "#2c3e50"), drop = F) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(text = element_text(size = 13, color = "#2e3440"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none",
        strip.background = element_rect(fill = "#e5e9f0")) +
  coord_cartesian(xlim = c(-1, 1))

# ------------------------------------------------------------------------------
# ----- 3. Figure 8 ------------------------------------------------------------
# ------------------------------------------------------------------------------

p2 <- ggplot(data_correlations) +
  geom_errorbar(aes(xmin = R_ab_marginal_low, xmax = R_ab_marginal_high, y = Genus, color = sign(R_ab_marginal_low) == sign(R_ab_marginal_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_ab_marginal, y = Genus, color = sign(R_ab_marginal_low) == sign(R_ab_marginal_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Correlation Coefficient (Viral Load vs. Bacterial Abundance)") +
  scale_color_manual(values = c("lightgrey", "#2c3e50"), drop = F) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(text = element_text(size = 13, color = "#2e3440"),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none",
        strip.background = element_rect(fill = "#e5e9f0")) +
  coord_cartesian(xlim = c(-1, 1))

plot <- p1 / p2

ggsave(here("figures", "Figure 8_9 raw.png"), plot = plot, dpi = 300, width = 7, height = 14)
ggsave(here("figures", "Figure 8_9 raw.svg"), plot = plot, dpi = 300, width = 7, height = 14, device = cairo_pdf)

# ------------------------------------------------------------------------------
# ----- 5. Figure 9 ------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1. Significant Presences ---------------------------------------------

# Significant presences:

# Streptococcus | None
# Bacillus | DCV

data_presence <- rbind(filter(data_clr, Genus == "Streptococcus", Condition == "None"),
                       filter(data_clr, Genus == "Bacillus", Condition == "DCV"))

data_presence$Genus <- ifelse(data_presence$Genus == "Streptococcus", "Streptococcus (Pre-infection)", "Bacillus (During infection)")

data_presence$Genus <- factor(data_presence$Genus, levels = c("Streptococcus (Pre-infection)", "Bacillus (During infection)"))


p1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(data = data_presence, aes(x = presence, y = foldchange), color = "#2C3E50") +
  facet_wrap(~Genus, nrow = 1) +
  scale_x_discrete(name = "Presence", labels = c("False", "True")) +
  scale_y_continuous(name = "Viral Load (Δ)", limits = c(-2, 10), expand = c(0,0),
                     breaks = seq(-2, 8, 2),
                     labels = c(expression(-10^2),
                                expression(10^0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6),
                                expression(10^8))) +
  theme_bw() +
  theme(text = element_text(size = 9, color = "#2e3440"),
        aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#e5e9f0")) +
  ggtitle("A.")

ggsave(here("figures", "Figure 10 raw.png"), plot = p1, dpi = 300, width = 5, height = 2.5)
ggsave(here("figures", "Figure 10 raw.svg"), plot = p1, dpi = 300, width = 5, height = 2.5)


