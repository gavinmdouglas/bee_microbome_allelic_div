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


compute_between_sample_mean_strain_divergence <- function(divergence_tab,
                                                          relabun_tab) {
  
  between_sample_divergence <- numeric()
  sample1_set <- character()
  sample2_set <- character()

  if (ncol(relabun_tab) > 1) {
    for (col_i in 1:(ncol(relabun_tab) - 1)) {
      sample_i <- colnames(relabun_tab)[col_i]
      sample_i_strains <- sort(rownames(relabun_tab)[which(relabun_tab[, sample_i] > 0)])
      
      for (col_j in (col_i + 1):ncol(relabun_tab)) {
        sample_j <- colnames(relabun_tab)[col_j]
        sample_j_strains <- sort(rownames(relabun_tab)[which(relabun_tab[, sample_j] > 0)])

        per_between_sample_divergence <- numeric()
        
        for (strain1 in sample_i_strains) {
          
          for (strain2 in sample_j_strains) {
            
            if (strain1 == strain2) { next }
            
            strain_combo <- paste(sort(c(strain1, strain2)), collapse = '_')
            
            per_between_sample_divergence <- c(per_between_sample_divergence,
                                               divergence_tab[strain_combo, 'percent_identity'])
          }
          
        }
        
        between_sample_divergence <- c(between_sample_divergence,
                                       mean(per_between_sample_divergence))
        sample1_set <- c(sample1_set, sample_i)
        sample2_set <- c(sample2_set, sample_j)
      }
    }
  }

  return(data.frame(sample1 = sample1_set,
                    sample2 = sample2_set,
                    mean_between_divergence = between_sample_divergence))
}

# Compare whether strains within samples are more or less divergent than those between samples.
strain_divergence <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strainfacts_and_ref_strain_hamming.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
strain_divergence <- strain_divergence[which(strain_divergence$comparable_length >= 1000), ]

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

# First compute for strains.
within_div_raw <- list()
between_div_raw <- list()
for (sp in names(strain_relabun)) {

  for (dataset in names(strain_relabun[[sp]])) {
    
    sp_dataset <- paste(sp, dataset, sep = '.')
    
    strain_divergence_subset <- strain_divergence[which(strain_divergence$annot == sp_dataset), ]
    rownames(strain_divergence_subset) <- paste(strain_divergence_subset$seq1, strain_divergence_subset$seq2, sep = '_')
    
    strain_relabun_subset <- strain_relabun[[sp]][[dataset]]
    rownames(strain_relabun_subset) <- paste(dataset, rownames(strain_relabun_subset), sep = '.')
    
    within_div_raw[[paste(sp, dataset, sep = '||')]] <- compute_within_sample_mean_strain_divergence(divergence_tab = strain_divergence_subset,
                                                                                                     relabun_tab = strain_relabun_subset)
    within_div_raw[[paste(sp, dataset, sep = '||')]]$species <- sp
    within_div_raw[[paste(sp, dataset, sep = '||')]]$dataset <- dataset
  
    
    between_div_raw[[paste(sp, dataset, sep = '||')]] <- compute_between_sample_mean_strain_divergence(divergence_tab = strain_divergence_subset,
                                                                                                       relabun_tab = strain_relabun_subset)
    between_div_raw[[paste(sp, dataset, sep = '||')]]$species <- sp
    between_div_raw[[paste(sp, dataset, sep = '||')]]$dataset <- dataset
    
    
  }

}

within_div <- do.call(rbind, within_div_raw)
rownames(within_div) <- NULL
within_div <- within_div[which(! is.na(within_div$mean_within_divergence)), ]

between_div <- do.call(rbind, between_div_raw)
rownames(between_div) <- NULL
between_div <- between_div[which(! is.na(between_div$mean_between_divergence)), ]

write.table(x = within_div,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/within_sample_strain_divergence.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')

write.table(x = between_div,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/between_sample_strain_divergence.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')
