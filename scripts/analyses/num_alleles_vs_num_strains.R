rm(list = ls(all.names = TRUE))

# Explore if the distribution of different alleles per gene departs from a
# basic model where every strain increases the probability that a different allele will be present proportionally.

compare_strains_vs_all_alleles <- function(strain_tab, all_alleles) {

  num_strains <- nrow(strain_tab)
  num_strains_per_sample <- colSums(strain_tab > 0)
  
  all_num_alleles <- numeric()
  sample_prevalence <- numeric()
  mean_num_alleles_when_present <- numeric()
  max_num_alleles_when_present <- numeric()
  min_num_alleles_when_present <- numeric()
  mean_ratio_num_alleles_strains <- numeric()

  full_p_out <- as.numeric()
  LLR <- as.numeric()
  
  cor_present_tau_out <- as.numeric()
  cor_present_p_out <- as.numeric()
  
  for (g in names(all_alleles)) {
    
    allele_comm <- all_alleles[[g]]
    
    num_alleles <- length(unique(allele_comm$strain))
    all_num_alleles <- c(all_num_alleles, num_alleles)
    
    num_alleles_per_sample <- table(allele_comm$sample)
    
    mean_num_alleles_when_present <- c(mean_num_alleles_when_present, mean(num_alleles_per_sample))
    min_num_alleles_when_present <- c(min_num_alleles_when_present, min(num_alleles_per_sample))
    max_num_alleles_when_present <- c(max_num_alleles_when_present, max(num_alleles_per_sample))
    
    # Only keep samples where the strains were called as present, ignore others.
    allele_comm <- allele_comm[which(allele_comm$sample %in% colnames(strain_tab)), ]
    num_alleles_per_sample <- table(allele_comm$sample)
    
    if (nrow(allele_comm) == 0 | num_alleles <= 2) {
      
      sample_prevalence <- c(sample_prevalence, 0)
      full_p_out <- c(full_p_out, NA)
      LLR <- c(LLR, NA)
      
      cor_present_tau_out <- c(cor_present_tau_out, NA)
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
    
    sample_prevalence <- c(sample_prevalence,
                           length(which(num_alleles_per_sample_w_filler > 0)))
    
    mean_ratio_num_alleles_strains <- c(mean_ratio_num_alleles_strains,
                                        mean(num_alleles_per_sample_w_filler / num_strains_per_sample))
    
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
      cor_present_out <- cor.test(num_alleles_per_sample,
                                  num_strains_per_sample[names(num_alleles_per_sample)],
                                  method = "kendall", exact = FALSE)
      cor_present_tau_out <- c(cor_present_tau_out, cor_present_out$estimate)
      cor_present_p_out <- c(cor_present_p_out, cor_present_out$p.value)
    } else {
      cor_present_tau_out <- c(cor_present_tau_out, NA)
      cor_present_p_out <- c(cor_present_p_out, NA) 
    }
  }

  strain_vs_allele_info <- data.frame(
                               Species = sp,
                               Gene = names(all_alleles),
                               num_alleles = all_num_alleles,
                               sample_prev = sample_prevalence,
                               mean_num_alleles_when_present = mean_num_alleles_when_present,
                               min_num_alleles_when_present = min_num_alleles_when_present,
                               max_num_alleles_when_present = max_num_alleles_when_present,
                               mean_ratio_num_alleles_strains_when_present = mean_ratio_num_alleles_strains,
                               multinomial_p = full_p_out,
                               multinomial_LLR = LLR,
                               tau_present = cor_present_tau_out,
                               cor_present_p = cor_present_p_out,
                               prop_samples_w_one_strain = length(which(num_strains_per_sample == 1)) / length(num_strains_per_sample),
                               mean_number_of_strains = mean(num_strains_per_sample))
  
  return(strain_vs_allele_info)

}
  
setwd('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/')

allele_relabun <- readRDS('accessory/RDS/strainfacts_accessory_allele_relabun.rds')
strain_relabun <- readRDS('core/RDS/strainfacts_core.genome_comm.rds')

raw_info <- list()

for (sp in names(strain_relabun)) {

  # Read in file of accessory genes to consider for this species (ignore others).
  sp_acc_gene_file <- paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory/',
                            sp,
                            '.txt',
                            sep = '')

  sp_acc_genes <- read.table(file = sp_acc_gene_file, header = FALSE, sep = '', stringsAsFactors = FALSE)$V1

  for (d in names(strain_relabun[[sp]])) {

    if (! sp %in% names(allele_relabun[[d]])) { next }

    # Subset allele set to specified set.
    gene_subset <- intersect(names(allele_relabun[[d]][[sp]]), sp_acc_genes)

    if (length(gene_subset) < 50) { next }

    allele_list <- allele_relabun[[d]][[sp]][gene_subset]

    strain_vs_allele_out <- compare_strains_vs_all_alleles(strain_tab = strain_relabun[[sp]][[d]],
                                                           allele_list)

    strain_vs_allele_out$Dataset <- d

    raw_info[[paste(sp, d, sep = '_')]] <- strain_vs_allele_out

  }  
}

allele_info <- do.call(rbind, raw_info)

write.table(x = allele_info,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/num_strain_vs_num_allele_breakdown.tsv',
            col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)
