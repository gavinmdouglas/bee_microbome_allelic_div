rm(list = ls(all.names = TRUE))

# Compare strain and allele presence/absence profiles, based on Jaccard distances.
# Calculate the minimum Jaccard distance vs strain value per allele (and return the mean).
# Also, return the mean of this statistic after randomizing the allele names per sample (and each sample scrambled independently).
# Do this scrambling for 1000 replicates per gene, and also compute the proportion of replicates with equal or lower metric values.

library(ggplot2)
library(parallel)
library(proxy)
library(reshape2)
library(Rfast)

compute_strain_vs_allele_jaccard <- function(strain_abun,
                                             allele_abun,
                                             intersecting_samples_only = TRUE) {
  
  allele_abun$strain <- paste0('allele_', allele_abun$strain)
  
  allele_abun_wide <- reshape2::dcast(data = allele_abun,
                                      formula =  sample ~ strain,
                                      value.var = 'community',
                                      fill = 0)
  rownames(allele_abun_wide) <- allele_abun_wide$sample
  allele_abun_wide <- allele_abun_wide[, -1, drop = FALSE]
  
  # Ignore allele samples without strain information (and drop any alleles missing after that).
  intersecting_samples <- intersect(rownames(allele_abun_wide), rownames(strain_abun))
  if (length(intersecting_samples) == 0) {
    return(list(overall_mean_of_min_jaccard = NA,
                overall_min_jaccard = NA,
                strain_vs_allele_jaccard = NA,
                num_samples = NA,
                num_alleles = NA,
                num_strains = NA,
                mean_sample_allelic_jaccard = NA,
                mean_sample_strain_jaccard = NA,
                mean_permuted_strain_vs_allele_jaccard = NA,
                permuted_mean_min_jaccard_P = NA))
  }
  allele_abun_wide <- allele_abun_wide[intersecting_samples, , drop = FALSE]
  allele_abun_wide <- allele_abun_wide[, which(colSums(allele_abun_wide) > 0), drop = FALSE]
  
  if (intersecting_samples_only) {
    strain_abun <- strain_abun[intersecting_samples, , drop = FALSE]
    strain_abun <- strain_abun[, which(colSums(strain_abun) > 0), drop = FALSE]
  } else {
    # Add in all samples with strains and lacking alleles into allele matrix, as all 0's (if needed).
    if (nrow(strain_abun) > nrow(allele_abun_wide)) {
      missing_samples <- setdiff(rownames(strain_abun), rownames(allele_abun_wide))
      dummy_allele <- matrix(0, nrow = length(missing_samples), ncol = ncol(allele_abun_wide))
      rownames(dummy_allele) <- missing_samples
      colnames(dummy_allele) <- colnames(allele_abun_wide)
      allele_abun_wide <- rbind(allele_abun_wide, dummy_allele)
    }
  }
  
  strain_vs_allele_jaccard <- proxy::dist(x = strain_abun,
                                          y = allele_abun_wide,
                                          method = 'Jaccard',
                                          by_rows = FALSE)
  
  min_jaccard_by_allele <- Rfast::colMins(x = strain_vs_allele_jaccard,
                                          value = TRUE)
  
  overall_mean_jaccard <- mean(min_jaccard_by_allele)
  overall_min_jaccard <- min(min_jaccard_by_allele)
  
  # Also get mean pairwise Jaccard distance of samples based on allele presence/absence.
  mean_sample_allelic_jaccard <- mean(proxy::dist(x = allele_abun_wide,
                                                  method = 'Jaccard',
                                                  by_rows = TRUE))

  # And same based on strains (which is important to compute each time, as it will change
  # depending on the sample subset).
  mean_sample_strain_jaccard <- mean(proxy::dist(x = strain_abun,
                                                 method = 'Jaccard',
                                                 by_rows = TRUE))
  
  # Finally, perform permutation test, where all frequencies are scrambled across each
  # sample independently.
  num_alleles <- ncol(allele_abun_wide)
  
  if (num_alleles > 1) {
  
    permuted_mean_min_jaccard <- as.numeric()
    
    num_reps <- 1000
    for (rep_i in 1:num_reps) {
     
      permuted_allele_abun_wide <- allele_abun_wide
      for (row_j in 1:nrow(permuted_allele_abun_wide)) {
        permuted_allele_abun_wide[row_j, ] <- as.numeric(sample(x = allele_abun_wide[row_j, ], size = num_alleles))
      }
      
      permuted_strain_vs_allele_jaccard <- proxy::dist(x = strain_abun,
                                                       y = permuted_allele_abun_wide,
                                                       method = 'Jaccard',
                                                       by_rows = FALSE)
      
      permuted_mean_min_jaccard <- c(permuted_mean_min_jaccard,
                                     mean(Rfast::colMins(x = permuted_strain_vs_allele_jaccard,
                                                         value = TRUE)))
  
    }
    
    overall_permuted_mean_min_jaccard <- mean(permuted_mean_min_jaccard)

    permuted_mean_min_jaccard_P <- (length(which(permuted_mean_min_jaccard <= overall_mean_jaccard)) + 1) / (num_reps + 1)
    
  } else {
   
    # Trivial if there is only one allele.
    overall_permuted_mean_min_jaccard <- overall_mean_jaccard
    
    permuted_mean_min_jaccard_P <- 1
     
  }
  
  return(list(overall_mean_of_min_jaccard = overall_mean_jaccard,
              overall_min_jaccard = overall_min_jaccard,
              num_samples = length(intersecting_samples),
              num_alleles = ncol(allele_abun_wide),
              num_strains = ncol(strain_abun),
              strain_vs_allele_jaccard = strain_vs_allele_jaccard,
              mean_sample_allelic_jaccard = mean_sample_allelic_jaccard,
              mean_sample_strain_jaccard = mean_sample_strain_jaccard,
              mean_permuted_strain_vs_allele_jaccard = overall_permuted_mean_min_jaccard,
              permuted_mean_min_jaccard_P = permuted_mean_min_jaccard_P))
  
}

all_strain_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
all_genes_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

raw_mean_min_list <- list()

for (species in names(all_strain_abun)) {
  
  print(species)
  
  species_strain_abun <- t(all_strain_abun[[species]])
  
  all_out <- parallel::mclapply(X = names(all_genes_abun[[species]]),
                                FUN = function(GENE) {
                                  compute_strain_vs_allele_jaccard(strain_abun = species_strain_abun,
                                                                   allele_abun = all_genes_abun[[species]][[GENE]])
                                },
                                mc.cores = 40)
  
  raw_mean_min <- vapply(all_out, function(x) { x$overall_mean_of_min_jaccard }, numeric(1))
  raw_overall_min <- vapply(all_out, function(x) { x$overall_min_jaccard }, numeric(1))
  raw_num_samples <- vapply(all_out, function(x) { x$num_samples }, integer(1))
  raw_num_alleles <- vapply(all_out, function(x) { x$num_alleles }, integer(1))
  raw_num_strains <- vapply(all_out, function(x) { x$num_strains }, integer(1))
  raw_mean_sample_allelic_jaccard <- vapply(all_out, function(x) { x$mean_sample_allelic_jaccard }, numeric(1))
  raw_mean_sample_strain_jaccard <- vapply(all_out, function(x) { x$mean_sample_strain_jaccard }, numeric(1))
  raw_mean_permuted_strain_vs_allele_jaccard <- vapply(all_out, function(x) { x$mean_permuted_strain_vs_allele_jaccard }, numeric(1))
  raw_permuted_mean_min_jaccard_P <- vapply(all_out, function(x) { x$permuted_mean_min_jaccard_P }, numeric(1))
  
  raw_mean_min_list[[species]] <- data.frame(species = species,
                                             gene = names(all_genes_abun[[species]]),
                                             overall_min = raw_overall_min,
                                             mean_min = raw_mean_min,
                                             num_samples = raw_num_samples,
                                             num_alleles = raw_num_alleles,
                                             num_strains = raw_num_strains,
                                             mean_sample_allelic_jaccard = raw_mean_sample_allelic_jaccard,
                                             mean_sample_strain_jaccard = raw_mean_sample_strain_jaccard,
                                             permuted_mean_min = raw_mean_permuted_strain_vs_allele_jaccard,
                                             permuted_mean_min_P = raw_permuted_mean_min_jaccard_P)
}

mean_min_out <- do.call(rbind, raw_mean_min_list)
mean_min_out <- mean_min_out[which(! is.na(mean_min_out$overall_min)), ]

write.table(x = mean_min_out,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
