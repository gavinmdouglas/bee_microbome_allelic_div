rm(list = ls(all.names = TRUE))

library(parallel)

compute_within_sample_mean_strain_divergence <- function(divergence_tab,
                                                         relabun_tab) {

  within_sample_divergence <- numeric()
  for (sampleid in colnames(relabun_tab)) {
    sample_strains <- sort(rownames(relabun_tab)[which(relabun_tab[, sampleid] > 0)])
    
    if (length(sample_strains) > 1) {
      
      per_within_sample_divergence <- numeric()
      
      for (i in 1:(length(sample_strains) - 1)) {
        
        strain1 <- sample_strains[i]
        
        for (j in (i + 1):length(sample_strains)) {
          
          strain2 <- sample_strains[j]
          
          if (strain1 == strain2) { next }
          
          strain_combo <- paste(strain1, strain2, sep = '_')
          
          per_within_sample_divergence <- c(per_within_sample_divergence,
                                            divergence_tab[strain_combo, 'percent_identity'])
          
        }
        
      }
      
      within_sample_divergence <- c(within_sample_divergence,
                                    mean(per_within_sample_divergence))
    } else {
      
      within_sample_divergence <- c(within_sample_divergence, NA)
      
    }
    
  }
  
  return(data.frame(sample = colnames(relabun_tab),
                    mean_within_divergence = within_sample_divergence))
  
}


compute_permuted_within_sample_mean_strain_divergence <- function(random_seed,
                                                                  divergence_tab,
                                                                  relabun_tab) {
  # Permute row names:
  set.seed(random_seed)
  rownames(relabun_tab) <- sample(rownames(relabun_tab))
  
  within_sample_divergence <- numeric()
  for (sampleid in colnames(relabun_tab)) {
    sample_strains <- sort(rownames(relabun_tab)[which(relabun_tab[, sampleid] > 0)])
    
    if (length(sample_strains) > 1) {
      
      per_within_sample_divergence <- numeric()
      
      for (i in 1:(length(sample_strains) - 1)) {
        
        strain1 <- sample_strains[i]
        
        for (j in (i + 1):length(sample_strains)) {
          
          strain2 <- sample_strains[j]
          
          if (strain1 == strain2) { next }
          
          strain_combo <- paste(strain1, strain2, sep = '_')
          
          per_within_sample_divergence <- c(per_within_sample_divergence,
                                            divergence_tab[strain_combo, 'percent_identity'])
          
        }
        
      }
      
      within_sample_divergence <- c(within_sample_divergence,
                                    mean(per_within_sample_divergence))
    } else {
      
      within_sample_divergence <- c(within_sample_divergence, NA)
      
    }
    
  }
  
  # Return just the overall mean of this replicate.
  return(mean(within_sample_divergence, na.rm = TRUE))
  
}

# Compare whether strains within samples are more or less divergent than those between samples.
strain_divergence <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strainfacts_and_ref_strain_hamming.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
strain_divergence <- strain_divergence[which(strain_divergence$comparable_length >= 1000), ]

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

num_reps <- 1000

raw_out <- list()

for (sp in names(strain_relabun)) {
  print(sp)
  for (dataset in names(strain_relabun[[sp]])) {
    print(dataset)
    sp_dataset <- paste(sp, dataset, sep = '.')

    strain_divergence_subset <- strain_divergence[which(strain_divergence$annot == sp_dataset), ]
    rownames(strain_divergence_subset) <- paste(strain_divergence_subset$seq1, strain_divergence_subset$seq2, sep = '_')

    strain_relabun_subset <- strain_relabun[[sp]][[dataset]]
    rownames(strain_relabun_subset) <- paste(dataset, rownames(strain_relabun_subset), sep = '.')

    sample_obs_mean_identity <- compute_within_sample_mean_strain_divergence(divergence_tab = strain_divergence_subset,
                                                                             relabun_tab = strain_relabun_subset)

    overall_obs_mean_identity <- mean(sample_obs_mean_identity$mean_within_divergence, na.rm = TRUE)
    
    permuted_mean_identities <- parallel::mclapply(X = 1:num_reps,
                                                   FUN = function(x) {
                                                     compute_permuted_within_sample_mean_strain_divergence(
                                                       random_seed = x,
                                                       divergence_tab = strain_divergence_subset,
                                                       relabun_tab = strain_relabun_subset)
                                                   },
                                                   mc.cores = 40)

    permuted_mean_identities <- unlist(permuted_mean_identities)

    p_lower <- (length(which(permuted_mean_identities <= overall_obs_mean_identity)) + 1) / (num_reps + 1)
    p_higher <- (length(which(permuted_mean_identities >= overall_obs_mean_identity)) + 1) / (num_reps + 1)

    raw_out[[paste(sp, dataset, sep = '||')]] <- data.frame(species = sp,
                                                                   dataset = dataset,
                                                                   obs_mean = overall_obs_mean_identity,
                                                                   permuted_mean = mean(permuted_mean_identities),
                                                                   p_lower = p_lower,
                                                                   p_higher = p_higher,
                                                                   num_strains = nrow(strain_relabun_subset),
                                                                   num_samples = ncol(strain_relabun_subset),
                                                                   mean_strains_per_sample = mean(colSums(strain_relabun_subset > 0)))
    
  }

}

combined_out <- do.call(rbind, raw_out)
rownames(combined_out) <- NULL
combined_out <- combined_out[which(! is.na(combined_out$obs_mean)), ]

write.table(x = combined_out,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/within_sample_strain_permutation_results.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')
