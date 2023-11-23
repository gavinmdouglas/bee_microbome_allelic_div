rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(parallel)

setwd("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/strainfacts/accessory")
# Explore if the distribution of different alleles per gene departs from a
# basic model where every strain increases the probability that a different allele will be present proportionally.

bcftools_hap_relabun <- readRDS("output/bcftools_based_RDS/bcftools_hap_relabun.rds")

bcftools_strain_info <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/strainfacts/core/output/bcftools.based.filt.comm.rds")

bcftools_sp_allele_info <- list()

for (sp in names(bcftools_strain_info)) {
  
  num_strains <- nrow(bcftools_strain_info[[sp]])
  num_strains_per_sample <- colSums(bcftools_strain_info[[sp]] > 0)
  
  all_num_alleles <- as.numeric()
  sample_prevalence <- as.numeric()
  mean_num_alleles_when_present <- as.numeric()
  max_num_alleles_when_present <- as.numeric()
  min_num_alleles_when_present <- as.numeric()
  mean_ratio_num_alleles_strains <- as.numeric()
  
  full_p_out <- as.numeric()
  LLR <- as.numeric()
  
  cor_present_rho_out <- as.numeric()
  cor_present_p_out <- as.numeric()
  
  for (g in names(bcftools_hap_relabun[[sp]])) {
    
    hap_comm <- bcftools_hap_relabun[[sp]][[g]]
    
    num_alleles <- length(unique(hap_comm$strain))
    all_num_alleles <- c(all_num_alleles, num_alleles)
    
    num_alleles_per_sample <- table(hap_comm$sample)
    
    mean_num_alleles_when_present <- c(mean_num_alleles_when_present, mean(num_alleles_per_sample))
    min_num_alleles_when_present <- c(min_num_alleles_when_present, min(num_alleles_per_sample))
    max_num_alleles_when_present <- c(max_num_alleles_when_present, max(num_alleles_per_sample))
    
    # Only keep samples where the strains were called as present, ignore others.
    hap_comm <- hap_comm[which(hap_comm$sample %in% colnames(bcftools_strain_info[[sp]])), ]
    num_alleles_per_sample <- table(hap_comm$sample)
    
    if (nrow(hap_comm) == 0 | num_alleles <= 2) {
      
      sample_prevalence <- c(sample_prevalence, 0)
      full_p_out <- c(full_p_out, NA)
      LLR <- c(LLR, NA)
      
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA)
      
      mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, NA)
      
      next
      
    }
    
    strain_samples_missing_alleles <- names(num_strains_per_sample)[which(! names(num_strains_per_sample) %in% names(num_alleles_per_sample))]
    
    if (length(strain_samples_missing_alleles) > 0) {
      missing_vec <- rep(0, length(strain_samples_missing_alleles))
      names(missing_vec) <- strain_samples_missing_alleles
      num_alleles_per_sample_w_filler <- c(num_alleles_per_sample, missing_vec)
      num_alleles_per_sample_w_filler <- num_alleles_per_sample_w_filler[names(num_strains_per_sample)]
    } else {
      num_alleles_per_sample_w_filler <- num_alleles_per_sample
    }
    
    sample_prevalence <- c(sample_prevalence, length(which(num_alleles_per_sample_w_filler > 0)))
    
    mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, mean(num_alleles_per_sample_w_filler / num_strains_per_sample))
    
    # Multinomial test based on all samples (including those where the gene is absent).
    exp_prob <- sum(num_alleles_per_sample_w_filler) / sum(num_strains_per_sample)
    exp_lik <- exp_prob * num_strains_per_sample
    exp_prob_dist <- exp_lik / sum(exp_lik)
    
    num_reps <- 10000
    multinomial_out <- XNomial::xmonte(obs = num_alleles_per_sample_w_filler, expr = exp_prob_dist, ntrials = num_reps, detail = 0)
    full_p <- ((multinomial_out$pLLR * num_reps) + 1) / (num_reps + 1)
    full_p_out <- c(full_p_out, full_p)
    LLR <- c(LLR, multinomial_out$observedLLR)
    
    if (length(num_alleles_per_sample) > 5) {
      cor_present_out <- cor.test(num_alleles_per_sample, num_strains_per_sample[names(num_alleles_per_sample)], method = "spearman", exact = FALSE)
      cor_present_rho_out <- c(cor_present_rho_out, cor_present_out$estimate)
      cor_present_p_out <- c(cor_present_p_out, cor_present_out$p.value)
    } else {
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA) 
    }
  }
  
  sp_allele_info <- data.frame(Species = sp,
                               Gene = names(bcftools_hap_relabun[[sp]]),
                               num_alleles = all_num_alleles,
                               sample_prev = sample_prevalence,
                               mean_num_alleles_when_present = mean_num_alleles_when_present,
                               min_num_alleles_when_present = min_num_alleles_when_present,
                               max_num_alleles_when_present = max_num_alleles_when_present,
                               mean_ratio_num_alleles_strains_when_present = mean_ratio_num_alleles_strains,
                               multinomial_p = full_p_out,
                               multinomial_LLR = LLR,
                               rho_present = cor_present_rho_out,
                               cor_present_p = cor_present_p_out)
  
  bcftools_sp_allele_info[[sp]] <- sp_allele_info
  
}

bcftools_allele_info <- do.call(rbind, bcftools_sp_allele_info)





varscan_hap_relabun <- readRDS("output/varscan_based_RDS/varscan_hap_relabun.rds")

varscan_strain_info <- readRDS("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/comp_mapping/strainfacts/core/output/varscan.based.filt.comm.rds")

varscan_sp_allele_info <- list()

for (sp in names(varscan_strain_info)) {
  
  num_strains <- nrow(varscan_strain_info[[sp]])
  num_strains_per_sample <- colSums(varscan_strain_info[[sp]] > 0)
  
  all_num_alleles <- as.numeric()
  sample_prevalence <- as.numeric()
  mean_num_alleles_when_present <- as.numeric()
  min_num_alleles_when_present <- as.numeric()
  max_num_alleles_when_present <- as.numeric()
  mean_ratio_num_alleles_strains <- as.numeric()
  
  full_p_out <- as.numeric()
  LLR <- as.numeric()
  
  cor_present_rho_out <- as.numeric()
  cor_present_p_out <- as.numeric()
  
  for (g in names(varscan_hap_relabun[[sp]])) {
    
    hap_comm <- varscan_hap_relabun[[sp]][[g]]
    
    num_alleles <- length(unique(hap_comm$strain))
    all_num_alleles <- c(all_num_alleles, num_alleles)
    
    num_alleles_per_sample <- table(hap_comm$sample)
    
    mean_num_alleles_when_present <- c(mean_num_alleles_when_present, mean(num_alleles_per_sample))
    min_num_alleles_when_present <- c(min_num_alleles_when_present, min(num_alleles_per_sample))
    max_num_alleles_when_present <- c(max_num_alleles_when_present, max(num_alleles_per_sample))
    
    # Only keep samples where the strains were called as present, ignore others.
    hap_comm <- hap_comm[which(hap_comm$sample %in% colnames(varscan_strain_info[[sp]])), ]
    num_alleles_per_sample <- table(hap_comm$sample)
    
    if (nrow(hap_comm) == 0 | num_alleles <= 2) {
      
      sample_prevalence <- c(sample_prevalence, 0)
      full_p_out <- c(full_p_out, NA)
      LLR <- c(LLR, NA)
      
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA)
      
      mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, NA)
      
      next
      
    }
    
    strain_samples_missing_alleles <- names(num_strains_per_sample)[which(! names(num_strains_per_sample) %in% names(num_alleles_per_sample))]
    
    if (length(strain_samples_missing_alleles) > 0) {
      missing_vec <- rep(0, length(strain_samples_missing_alleles))
      names(missing_vec) <- strain_samples_missing_alleles
      num_alleles_per_sample_w_filler <- c(num_alleles_per_sample, missing_vec)
      num_alleles_per_sample_w_filler <- num_alleles_per_sample_w_filler[names(num_strains_per_sample)]
    } else {
      num_alleles_per_sample_w_filler <- num_alleles_per_sample
    }
    
    sample_prevalence <- c(sample_prevalence, length(which(num_alleles_per_sample_w_filler > 0)))
    
    mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains, mean(num_alleles_per_sample_w_filler / num_strains_per_sample))
    
    # Multinomial test based on all samples (including those where the gene is absent).
    exp_prob <- sum(num_alleles_per_sample_w_filler) / sum(num_strains_per_sample)
    exp_lik <- exp_prob * num_strains_per_sample
    exp_prob_dist <- exp_lik / sum(exp_lik)
    
    num_reps <- 10000
    multinomial_out <- XNomial::xmonte(obs = num_alleles_per_sample_w_filler, expr = exp_prob_dist, ntrials = num_reps, detail = 0)
    full_p <- ((multinomial_out$pLLR * num_reps) + 1) / (num_reps + 1)
    full_p_out <- c(full_p_out, full_p)
    LLR <- c(LLR, multinomial_out$observedLLR)
    
    if (length(num_alleles_per_sample) > 5) {
      cor_present_out <- cor.test(num_alleles_per_sample, num_strains_per_sample[names(num_alleles_per_sample)], method = "spearman", exact = FALSE)
      cor_present_rho_out <- c(cor_present_rho_out, cor_present_out$estimate)
      cor_present_p_out <- c(cor_present_p_out, cor_present_out$p.value)
    } else {
      cor_present_rho_out <- c(cor_present_rho_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA) 
    }
  }
  
  sp_allele_info <- data.frame(Species = sp,
                               Gene = names(varscan_hap_relabun[[sp]]),
                               num_alleles = all_num_alleles,
                               sample_prev = sample_prevalence,
                               mean_num_alleles_when_present = mean_num_alleles_when_present,
                               min_num_alleles_when_present = min_num_alleles_when_present,
                               max_num_alleles_when_present = max_num_alleles_when_present,
                               mean_ratio_num_alleles_strains_when_present = mean_ratio_num_alleles_strains,
                               multinomial_p = full_p_out,
                               multinomial_LLR = LLR,
                               rho_present = cor_present_rho_out,
                               cor_present_p = cor_present_p_out)
  
  varscan_sp_allele_info[[sp]] <- sp_allele_info
  
}

varscan_allele_info <- do.call(rbind, varscan_sp_allele_info)



ggplot(data = varscan_allele_info, aes(x = sample_prev, y = min_num_alleles_when_present)) +
  geom_point() +
  facet_wrap(Species ~ .)



bcftools_allele_info_full_p <- ggplot(data = bcftools_allele_info, mapping = aes(y = Species, x = full_p)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("bcftools")


varscan_allele_info_full_p <- ggplot(data = varscan_allele_info, mapping = aes(y = Species, x = full_p)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("varscan")

plot_grid(bcftools_allele_info_full_p, varscan_allele_info_full_p)



bcftools_allele_info_rho_present <- ggplot(data = bcftools_allele_info, mapping = aes(y = Species, x = rho_present)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("bcftools")


varscan_allele_info_rho_present <- ggplot(data = varscan_allele_info, mapping = aes(y = Species, x = rho_present)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("varscan")

plot_grid(bcftools_allele_info_rho_present, varscan_allele_info_rho_present)




bcftools_allele_info_mean_ratio_num_alleles_strains_when_present <- ggplot(data = bcftools_allele_info, mapping = aes(y = Species, x = mean_ratio_num_alleles_strains_when_present)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("bcftools")


varscan_allele_info_mean_ratio_num_alleles_strains_when_present <- ggplot(data = varscan_allele_info, mapping = aes(y = Species, x = mean_ratio_num_alleles_strains_when_present)) +
  geom_quasirandom(col = "black", groupOnX = FALSE) +
  geom_boxplot(fill = "grey", alpha = 0.8, outlier.shape = NA) +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  ggtitle("varscan")

plot_grid(bcftools_allele_info_mean_ratio_num_alleles_strains_when_present, varscan_allele_info_mean_ratio_num_alleles_strains_when_present)



test <- lm(formula = multinomial_LLR ~ num_alleles + sample_prev + Species,
           data = varscan_allele_info)

tmp <- nls(formula = multinomial_LLR ~ SSlogis(num_alleles + sample_prev + Species),
           data = varscan_allele_info)


allele_info <- list(bcftools_allele_info = bcftools_allele_info,
                    varscan_allele_info = varscan_allele_info)

saveRDS(object = allele_info, file = "/home/gdouglas/tmp/allele_info.rds")