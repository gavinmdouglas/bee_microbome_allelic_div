rm(list = ls(all.names = TRUE))

library(reshape2)
library(parallel)

compute_within_sample_mean_allele_divergence <- function(divergence_tab,
                                                         relabun_tab) {
  
  unique_samples <- unique(relabun_tab$sample)
  
  within_sample_divergence <- numeric()
  for (sampleid in unique_samples) {
    sample_strains <- sort(relabun_tab[which(relabun_tab$sample == sampleid), 'strain'])
    
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
  
  return(data.frame(sample = unique_samples,
                    mean_within_divergence = within_sample_divergence))
  
}


compute_permuted_within_sample_mean_allele_divergence <- function(random_seed,
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

allele_divergence <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/strainfacts_allele_hamming.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
allele_divergence <- allele_divergence[which(allele_divergence$comparable_length >= 100), ]

allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

prepped_input <- list()

for (dataset in names(allele_relabun)) {
  print(dataset)
  for (sp in names(allele_relabun[[dataset]])) {
    print(sp)
    for (gene in names(allele_relabun[[dataset]][[sp]])) {
      
      gene_dataset <- paste(gene, dataset, sep = '.')
      
      if (! gene_dataset %in% allele_divergence$annot) { next }
      
      allele_divergence_subset <- allele_divergence[which(allele_divergence$annot == gene_dataset), ]
      rownames(allele_divergence_subset) <- paste(allele_divergence_subset$seq1, allele_divergence_subset$seq2, sep = '_')
      
      allele_relabun_subset <- allele_relabun[[dataset]][[sp]][[gene]]
      allele_relabun_subset$strain <- paste('allele', as.character(allele_relabun_subset$strain), sep = '_')
      rownames(allele_relabun_subset) <- paste(dataset, rownames(allele_relabun_subset), sep = '.')
      
      # Then convert the allele tab to wide format.
      allele_relabun_subset_wide <- reshape2::dcast(data = allele_relabun_subset,
                                                    formula = strain ~ sample,
                                                    fill = 0,
                                                    value.var = 'community')
      rownames(allele_relabun_subset_wide) <- allele_relabun_subset_wide$strain
      allele_relabun_subset_wide <- allele_relabun_subset_wide[, -1, drop = FALSE]
      
      prepped_input[[gene_dataset]] <- list(allele_divergence_subset = allele_divergence_subset,
                                            allele_relabun_subset = allele_relabun_subset,
                                            allele_relabun_subset_wide = allele_relabun_subset_wide,
                                            sp = sp,
                                            dataset = dataset,
                                            gene = gene)
      
    }
  }
}

raw_out <- parallel::mclapply(prepped_input,
                              function(x) {
                                
                                sample_obs_mean_identity <- compute_within_sample_mean_allele_divergence(divergence_tab = x$allele_divergence_subset,
                                                                                                         relabun_tab = x$allele_relabun_subset)
                                
                                overall_obs_mean_identity <- mean(sample_obs_mean_identity$mean_within_divergence, na.rm = TRUE)
                                
                                
                                permuted_mean_identities <- numeric()
                                
                                num_reps <- 1000
                                
                                for (i in 1:num_reps) {
                                  permuted_mean_identities <- c(
                                    permuted_mean_identities,
                                    compute_permuted_within_sample_mean_allele_divergence(
                                      random_seed = i,
                                      divergence_tab = x$allele_divergence_subset,
                                      relabun_tab = x$allele_relabun_subset_wide))
                                }

                                p_lower <- (length(which(permuted_mean_identities <= overall_obs_mean_identity)) + 1) / (num_reps + 1)
                                p_higher <- (length(which(permuted_mean_identities >= overall_obs_mean_identity)) + 1) / (num_reps + 1)
                                
                                return(data.frame(species = x$sp,
                                                  dataset = x$dataset,
                                                  gene = x$gene,
                                                  obs_mean = overall_obs_mean_identity,
                                                  permuted_mean = mean(permuted_mean_identities),
                                                  p_lower = p_lower,
                                                  p_higher = p_higher,
                                                  num_alleles = nrow(x$allele_relabun_subset_wide),
                                                  num_samples = ncol(x$allele_relabun_subset_wide),
                                                  mean_alleles_per_sample = mean(colSums(x$allele_relabun_subset_wide > 0))))
                                
                                },
                              mc.cores = 55)

combined_out <- do.call(rbind, raw_out)
rownames(combined_out) <- NULL
combined_out <- combined_out[which(! is.na(combined_out$obs_mean)), ]

write.table(x = combined_out,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/within_sample_allele_permutation_results.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')
