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

data_clr <- filter(data_clr, Genus %in% unique(data_clr$Genus))

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
      
      
      
      # Presence Phylo
      var_pr_phylo_viral <- model_phylo_presence$VCV[,1]
      var_pr_phylo_pr <- model_phylo_presence$VCV[,4]
      cov_pr_phylo_pr <- model_phylo_presence$VCV[,2]
      B_pr_phylo <- as.mcmc(cov_pr_phylo_pr / var_pr_phylo_pr)
      R_pr_phylo <- as.mcmc(cov_pr_phylo_pr /(sqrt(var_pr_phylo_viral * var_pr_phylo_pr)))
      
      # Abundance Phylo
      var_ab_phylo_viral <- model_phylo_abundance$VCV[,1]
      var_ab_phylo_ab <- model_phylo_abundance$VCV[,4]
      cov_ab_phylo_ab <- model_phylo_abundance$VCV[,2]
      ab_phylo_viral <- model_phylo_abundance$Sol[,1]
      ab_phylo_ab <- model_phylo_abundance$Sol[,2]
      B_ab_phylo <- as.mcmc(cov_ab_phylo_ab / var_ab_phylo_ab)
      R_ab_phylo <- as.mcmc(cov_ab_phylo_ab /(sqrt(var_ab_phylo_viral * var_ab_phylo_ab)))
      C_ab_phylo <- ab_phylo_viral - (B_ab_phylo * ab_phylo_ab)
      
      # Presence Nonphylo
      var_pr_nonphylo_viral <- model_nonphylo_presence$VCV[,1]
      var_pr_nonphylo_pr <- model_nonphylo_presence$VCV[,4]
      cov_pr_nonphylo_pr <- model_nonphylo_presence$VCV[,2]
      B_pr_nonphylo <- as.mcmc(cov_pr_nonphylo_pr / var_pr_nonphylo_pr)
      R_pr_nonphylo <- as.mcmc(cov_pr_nonphylo_pr /(sqrt(var_pr_nonphylo_viral * var_pr_nonphylo_pr)))
      
      # Abundance Nonphylo
      var_ab_nonphylo_viral <- model_nonphylo_abundance$VCV[,1]
      var_ab_nonphylo_ab <- model_nonphylo_abundance$VCV[,4]
      cov_ab_nonphylo_ab <- model_nonphylo_abundance$VCV[,2]
      ab_nonphylo_viral <- model_nonphylo_abundance$Sol[,1]
      ab_nonphylo_ab <- model_nonphylo_abundance$Sol[,2]
      B_ab_nonphylo <- as.mcmc(cov_ab_nonphylo_ab / var_ab_nonphylo_ab)
      R_ab_nonphylo <- as.mcmc(cov_ab_nonphylo_ab /(sqrt(var_ab_nonphylo_viral * var_ab_nonphylo_ab)))
      C_ab_nonphylo <- ab_nonphylo_viral - (B_ab_nonphylo * ab_nonphylo_ab)
      
      data_correlations <- rbind(data_correlations,
                                 data.frame(Genus = g,
                                            Condition = c,
                                            B_pr_phylo = mean(B_pr_phylo),
                                            B_pr_phylo_low = HPDinterval(B_pr_phylo)[1],
                                            B_pr_phylo_high = HPDinterval(B_pr_phylo)[2],
                                            R_pr_phylo = mean(R_pr_phylo),
                                            R_pr_phylo_low = HPDinterval(R_pr_phylo)[1],
                                            R_pr_phylo_high = HPDinterval(R_pr_phylo)[2],
                                            B_ab_phylo = mean(B_ab_phylo),
                                            B_ab_phylo_low = HPDinterval(B_ab_phylo)[1],
                                            B_ab_phylo_high = HPDinterval(B_ab_phylo)[2],
                                            R_ab_phylo = mean(R_ab_phylo),
                                            R_ab_phylo_low = HPDinterval(R_ab_phylo)[1],
                                            R_ab_phylo_high = HPDinterval(R_ab_phylo)[2],
                                            C_ab_phylo = mean(C_ab_phylo),
                                            C_ab_phylo_low = HPDinterval(C_ab_phylo)[1],
                                            C_ab_phylo_high = HPDinterval(C_ab_phylo)[2],
                                            B_pr_nonphylo = mean(B_pr_nonphylo),
                                            B_pr_nonphylo_low = HPDinterval(B_pr_nonphylo)[1],
                                            B_pr_nonphylo_high = HPDinterval(B_pr_nonphylo)[2],
                                            R_pr_nonphylo = mean(R_pr_nonphylo),
                                            R_pr_nonphylo_low = HPDinterval(R_pr_nonphylo)[1],
                                            R_pr_nonphylo_high = HPDinterval(R_pr_nonphylo)[2],
                                            B_ab_nonphylo = mean(B_ab_nonphylo),
                                            B_ab_nonphylo_low = HPDinterval(B_ab_nonphylo)[1],
                                            B_ab_nonphylo_high = HPDinterval(B_ab_nonphylo)[2],
                                            R_ab_nonphylo = mean(R_ab_nonphylo),
                                            R_ab_nonphylo_low = HPDinterval(R_ab_nonphylo)[1],
                                            R_ab_nonphylo_high = HPDinterval(R_ab_nonphylo)[2],
                                            C_ab_nonphylo = mean(C_ab_nonphylo),
                                            C_ab_nonphylo_low = HPDinterval(C_ab_nonphylo)[1],
                                            C_ab_nonphylo_high = HPDinterval(C_ab_nonphylo)[2]))
      
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
  geom_errorbar(aes(xmin = R_pr_phylo_low, xmax = R_pr_phylo_high, y = Genus, color = sign(R_pr_phylo_low) == sign(R_pr_phylo_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_pr_phylo, y = Genus, color = sign(R_pr_phylo_low) == sign(R_pr_phylo_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Phylogenetic Correlation Coefficient") +
  scale_color_manual(values = c("lightgrey", "#2c3e50"), drop = F) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(text = element_text(size = 13, color = "#2e3440"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#e5e9f0")) +
  coord_cartesian(xlim = c(-1, 1))

p2 <- ggplot(data_correlations) +
  geom_errorbar(aes(xmin = R_pr_nonphylo_low, xmax = R_pr_nonphylo_high, y = Genus, color = sign(R_pr_nonphylo_low) == sign(R_pr_nonphylo_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_pr_nonphylo, y = Genus, color = sign(R_pr_nonphylo_low) == sign(R_pr_nonphylo_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Non-phylogenetic Correlation Coefficient") +
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

plot <- p1 | p2

ggsave(here("figures", "Figure 8 raw.png"), plot = plot, dpi = 300, width = 7, height = 7)
ggsave(here("figures", "Figure 8 raw.svg"), plot = plot, dpi = 300, width = 7, height = 7)

# ------------------------------------------------------------------------------
# ----- 3. Figure 8 ------------------------------------------------------------
# ------------------------------------------------------------------------------

p1 <- ggplot(data_correlations) +
  geom_errorbar(aes(xmin = R_ab_phylo_low, xmax = R_ab_phylo_high, y = Genus, color = sign(R_ab_phylo_low) == sign(R_ab_phylo_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_ab_phylo, y = Genus, color = sign(R_ab_phylo_low) == sign(R_ab_phylo_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Phylogenetic Correlation Coefficient") +
  scale_color_manual(values = c("lightgrey", "#2c3e50"), drop = F) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(text = element_text(size = 13, color = "#2e3440"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = "#e5e9f0")) +
  coord_cartesian(xlim = c(-1, 1))

p2 <- ggplot(data_correlations) +
  geom_errorbar(aes(xmin = R_ab_nonphylo_low, xmax = R_ab_nonphylo_high, y = Genus, color = sign(R_ab_nonphylo_low) == sign(R_ab_nonphylo_high)), width = 0, linewidth = 1) +
  geom_point(aes(x = R_ab_nonphylo, y = Genus, color = sign(R_ab_nonphylo_low) == sign(R_ab_nonphylo_high)), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_y_discrete(drop = F, position = "right") +
  scale_x_continuous(name = "Non-phylogenetic Correlation Coefficient") +
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

plot <- p1 | p2

ggsave(here("figures", "Figure 9 raw.png"), plot = plot, dpi = 300, width = 7, height = 7)
ggsave(here("figures", "Figure 9 raw.svg"), plot = plot, dpi = 300, width = 7, height = 7)

# ------------------------------------------------------------------------------
# ----- 5. Figure 9 ------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1. Significant Presences ---------------------------------------------

# Significant presences:

# Streptococcus | None
# Bacillus | DCV
# Pseudomonas | DCV

data_presence <- rbind(filter(data_clr, Genus == "Streptococcus", Condition == "None"),
                       filter(data_clr, Genus == "Bacillus", Condition == "DCV"),
                       filter(data_clr, Genus == "Pseudomonas", Condition == "DCV"))

data_presence$Genus <- ifelse(data_presence$Genus == "Streptococcus", "Streptococcus (Pre-infection)",
                              ifelse(data_presence$Genus == "Bacillus", "Bacillus (During infection)", "Pseudomonas (During infection)"))

data_presence$Genus <- factor(data_presence$Genus, levels = c("Streptococcus (Pre-infection)", "Bacillus (During infection)", "Pseudomonas (During infection)"))


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


# ----- 5.2. Significant Abundances --------------------------------------------

# Significant abundances:
# Ralstonia | None
# Acetobacter | DCV

data_abundance <- rbind(filter(data_clr, Genus == "Ralstonia", Condition == "None", presence == T),
                        filter(data_clr, Genus == "Acetobacter", Condition == "DCV", presence == T))

data_abundance$Genus <- ifelse(data_abundance$Genus == "Ralstonia", "Ralstonia (Pre-infection)", "Acetobacter (During infection)")

data_abundance$Genus <- factor(data_abundance$Genus, levels = c("Ralstonia (Pre-infection)", "Acetobacter (During infection)"))

trendlines <- data.frame(Genus = c("Ralstonia (Pre-infection)", "Acetobacter (During infection)"),
                         B = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$B_ab_nonphylo,
                               filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$B_ab_nonphylo),
                         B_low = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$B_ab_nonphylo_low,
                                   filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$B_ab_nonphylo_low),
                         B_high = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$B_ab_nonphylo_high,
                                    filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$B_ab_nonphylo_high),
                         C = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$C_ab_nonphylo,
                               filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$C_ab_nonphylo),
                         C_low = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$C_ab_nonphylo_low,
                                   filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$C_ab_nonphylo_low),
                         C_high = c(filter(data_correlations, Genus == "Ralstonia", Condition == "None")$C_ab_nonphylo_high,
                                    filter(data_correlations, Genus == "Acetobacter", Condition == "DCV")$C_ab_nonphylo_high))

trendlines$Genus <- factor(trendlines$Genus, levels = c("Ralstonia (Pre-infection)", "Acetobacter (During infection)"))

x_vals <- seq(-1, 12, length.out = 100)

ribbon_df <- trendlines %>%
  rowwise() %>%
  do({
    genus <- .$Genus
    tibble(Genus = genus,
           abundance = x_vals,
           y_low = .$B_low * x_vals + .$C_low,
           y_high = .$B_high * x_vals + .$C_high)}) %>%
  ungroup()

data_abundance_spc <- data_abundance %>% group_by(Genus, animal) %>%
  summarise(mean_ab = mean(abundance),
            mean_vir = mean(foldchange))

p2a <- ggplot() +
  geom_ribbon(data = filter(ribbon_df, Genus == "Ralstonia (Pre-infection)"), aes(x = abundance, ymin = y_low, ymax = y_high), fill = "#BDC3C7", alpha = 0.5) +
  geom_abline(data = filter(trendlines, Genus == "Ralstonia (Pre-infection)"), aes(slope = B, intercept = C), color = "#0064FF") +
  geom_point(data = filter(data_abundance_spc, Genus == "Ralstonia (Pre-infection)", animal != "D.simulans"), aes(x = mean_ab, y = mean_vir), alpha = 0.5) +
  facet_wrap(~Genus) +
  scale_y_continuous(name = "Viral Load (Δ)", breaks = seq(-2, 8, 2),
                     labels = c(expression(-10^2),
                                expression(10^0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6),
                                expression(10^8))) +
  scale_x_continuous(name = "CLR-transformed Abundance") +
  theme_bw() +
  theme(text = element_text(size = 9, color = "#2e3440"),
        aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#e5e9f0")) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 8)) +
  ggtitle("B.")

p2b <- ggplot() +
  geom_ribbon(data = filter(ribbon_df, Genus == "Acetobacter (During infection)"), aes(x = abundance, ymin = y_low, ymax = y_high), fill = "#BDC3C7", alpha = 0.5) +
  geom_abline(data = filter(trendlines, Genus == "Acetobacter (During infection)"), aes(slope = B, intercept = C), color = "#0064FF") +
  geom_point(data = filter(data_abundance_spc, Genus == "Acetobacter (During infection)", !animal %in% c("D.euronotus", "D.sucinea")), aes(x = mean_ab, y = mean_vir), alpha = 0.5) +
  facet_wrap(~Genus) +
  scale_y_continuous(name = "Viral Load (Δ)", breaks = seq(-2, 8, 2),
                     labels = c(expression(-10^2),
                                expression(10^0),
                                expression(10^2),
                                expression(10^4),
                                expression(10^6),
                                expression(10^8))) +
  scale_x_continuous(name = "CLR-transformed Abundance") +
  theme_bw() +
  theme(text = element_text(size = 9, color = "#2e3440"),
        aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "#e5e9f0"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(5, 9), ylim = c(0, 8))

# ----- 5.3. Save Raw Figure ---------------------------------------------------

plot2 <- p2a | p2b

ggsave(here("figures", "Figure 10a raw.png"), plot = p1, dpi = 300, width = 7, height = 2.5)
ggsave(here("figures", "Figure 10a raw.svg"), plot = p1, dpi = 300, width = 7, height = 2.5)

ggsave(here("figures", "Figure 10b raw.png"), plot = plot2, dpi = 300, width = 7, height = 2.5)
ggsave(here("figures", "Figure 10b raw.svg"), plot = plot2, dpi = 300, width = 7, height = 2.5)
