rm(list = ls(all.names = TRUE))

# Compare strain and allele presence/absence profiles, based on Jaccard distances.
# Calculate the minimum Jaccard distance vs strain value per allele (and return the mean).
# Also, return the mean of this statistic after randomizing the allele names per sample (and each sample scrambled independently).
# Do this scrambling for 1000 replicates per gene, and also compute the proportion of replicates with equal or lower metric values.

# Note that this should be done for both strains based on StrainFacts and based on StrainGST.
# The same subset of samples should be compared for each, to make sure this is a fair comparison.

# Note that allele-level inferences were only made with StrainFacts, as although StrainGE/GST can give you SNV calls, these are not
# phased into alleles.

library(ggplot2)
library(parallel)
library(proxy)
library(reshape2)
library(Rfast)

strain_vs_allele_jaccard_and_permute <- function(strain_abun,
                                                 allele_abun) {
  # Compute key strain vs allele jaccard
  strain_vs_allele_jaccard <- proxy::dist(x = strain_abun,
                                          y = allele_abun,
                                          method = 'Jaccard',
                                          by_rows = FALSE)
  
  min_jaccard_by_allele <- Rfast::colMins(x = strain_vs_allele_jaccard,
                                          value = TRUE)
  
  overall_mean_jaccard <- mean(min_jaccard_by_allele)
  overall_min_jaccard <- min(min_jaccard_by_allele)
  
  # Get Jaccard based on strains (which is important to compute each time,
  # as it will change depending on the sample subset).
  # Possibly included as a confounding variable to help control
  # for the overall strain variation.
  mean_sample_strain_jaccard <- mean(proxy::dist(x = strain_abun,
                                                 method = 'Jaccard',
                                                 by_rows = TRUE))
  
  # Finally, perform permutation test, where all frequencies are scrambled across each
  # sample independently.
  num_alleles <- ncol(allele_abun)

  if (num_alleles > 1) {
    
    permuted_mean_min_jaccard <- as.numeric()
    
    num_reps <- 1000
    for (rep_i in 1:num_reps) {
      
      permuted_allele_abun_wide <- allele_abun
      for (row_j in 1:nrow(permuted_allele_abun_wide)) {
        permuted_allele_abun_wide[row_j, ] <- as.numeric(sample(x = allele_abun[row_j, ], size = num_alleles))
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
  
  return(list(overall_mean_jaccard=overall_mean_jaccard,
              overall_min_jaccard=overall_min_jaccard,
              mean_sample_strain_jaccard=mean_sample_strain_jaccard,
              overall_permuted_mean_min_jaccard=overall_permuted_mean_min_jaccard,
              permuted_mean_min_jaccard_P=permuted_mean_min_jaccard_P))
}

compute_strain_vs_allele_jaccard <- function(strainfacts_abun,
                                             straingst_abun,
                                             allele_abun) {
  
  allele_abun$strain <- paste0('allele_', allele_abun$strain)
  
  allele_abun_wide <- reshape2::dcast(data = allele_abun,
                                      formula =  sample ~ strain,
                                      value.var = 'community',
                                      fill = 0)

  rownames(allele_abun_wide) <- allele_abun_wide$sample
  allele_abun_wide <- allele_abun_wide[, -1, drop = FALSE]
  
  # Ignore allele samples without strain information (and drop any alleles missing after that).
  intersecting_samples <- intersect(rownames(allele_abun_wide), rownames(strainfacts_abun))
  intersecting_samples <- intersect(intersecting_samples, rownames(straingst_abun))

  if (length(intersecting_samples) == 0) {
    return(list(num_samples = NA,
                num_alleles = NA,
                num_strainfacts_strains = NA,
                num_straingst_strains = NA,
                
                mean_sample_allelic_jaccard = NA,
                
                strainfacts_mean_of_min_jaccard = NA,
                strainfacts_min_jaccard = NA,
                mean_strainfacts_vs_allele_jaccard = NA,
                mean_strainfacts_jaccard = NA,
                mean_permuted_strainfacts_vs_allele_jaccard = NA,
                permuted_strainfacts_mean_min_jaccard_P = NA,
                
                straingst_mean_of_min_jaccard = NA,
                straingst_min_jaccard = NA,
                mean_straingst_vs_allele_jaccard = NA,
                mean_straingst_jaccard = NA,
                mean_permuted_straingst_vs_allele_jaccard = NA,
                permuted_straingst_mean_min_jaccard_P = NA))
  }

  allele_abun_wide <- allele_abun_wide[intersecting_samples, , drop = FALSE]
  allele_abun_wide <- allele_abun_wide[, which(colSums(allele_abun_wide) > 0), drop = FALSE]
  
  strainfacts_abun <- strainfacts_abun[intersecting_samples, , drop = FALSE]
  strainfacts_abun <- strainfacts_abun[, which(colSums(strainfacts_abun) > 0), drop = FALSE]
  
  straingst_abun <- straingst_abun[intersecting_samples, , drop = FALSE]
  straingst_abun <- straingst_abun[, which(colSums(straingst_abun) > 0), drop = FALSE]
  
  # Get Jaccard summaries for StrainGE vs alleles and StrainFacts vs alleles separately.
  strainfacts_jaccard_out <- strain_vs_allele_jaccard_and_permute(strain_abun = strainfacts_abun,
                                                                   allele_abun = allele_abun_wide)

  straingst_jaccard_out <- strain_vs_allele_jaccard_and_permute(strain_abun = straingst_abun,
                                                               allele_abun = allele_abun_wide)
  
  # Also get mean pairwise Jaccard distance of samples based on allele presence/absence.
  mean_sample_allelic_jaccard <- mean(proxy::dist(x = allele_abun_wide,
                                                  method = 'Jaccard',
                                                  by_rows = TRUE))

  return(list(num_samples = length(intersecting_samples),
              num_alleles = ncol(allele_abun_wide),
              num_strainfacts_strains = ncol(strainfacts_abun),
              num_straingst_strains = ncol(straingst_abun),
              
              mean_sample_allelic_jaccard = mean_sample_allelic_jaccard,
              
              strainfacts_mean_of_min_jaccard = strainfacts_jaccard_out$overall_mean_jaccard,
              strainfacts_min_jaccard = strainfacts_jaccard_out$overall_min_jaccard,
              mean_strainfacts_vs_allele_jaccard = strainfacts_jaccard_out$mean_sample_strain_jaccard,
              mean_strainfacts_jaccard = strainfacts_jaccard_out$mean_sample_strain_jaccard,
              mean_permuted_strainfacts_vs_allele_jaccard = strainfacts_jaccard_out$overall_permuted_mean_min_jaccard,
              permuted_strainfacts_mean_min_jaccard_P = strainfacts_jaccard_out$permuted_mean_min_jaccard_P,
              
              straingst_mean_of_min_jaccard = straingst_jaccard_out$overall_mean_jaccard,
              straingst_min_jaccard = straingst_jaccard_out$overall_min_jaccard,
              mean_straingst_vs_allele_jaccard = straingst_jaccard_out$mean_sample_strain_jaccard,
              mean_straingst_jaccard = straingst_jaccard_out$mean_sample_strain_jaccard,
              mean_permuted_straingst_vs_allele_jaccard = straingst_jaccard_out$overall_permuted_mean_min_jaccard,
              permuted_straingst_mean_min_jaccard_P = straingst_jaccard_out$permuted_mean_min_jaccard_P))

}

all_straingst_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')
all_strainfacts_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
all_genes_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

raw_mean_min_list <- list()

for (species in names(all_strainfacts_abun)) {
  
  for (dataset in names(all_strainfacts_abun[[species]])) {
  
    print(species)
    
    if (! dataset %in% names(all_genes_abun)) { next }
    if (! species %in% names(all_genes_abun[[dataset]])) { next }
    if (! species %in% names(all_straingst_abun)) { next }
    if (! dataset %in% names(all_straingst_abun[[species]])) { next }
    
    # Skip dataset/species combinations where there are fewer than 50 genes with inferred alleles.
    if (length(all_genes_abun[[dataset]][[species]]) < 50) { next }
    
    species_strainfacts_abun <- t(all_strainfacts_abun[[species]][[dataset]])
    species_straingst_abun <- t(all_straingst_abun[[species]][[dataset]])
    
    all_out <- parallel::mclapply(X = names(all_genes_abun[[dataset]][[species]]),
                                  FUN = function(GENE) {
                                    compute_strain_vs_allele_jaccard(strainfacts_abun = species_strainfacts_abun,
                                                                     straingst_abun = species_straingst_abun,
                                                                     allele_abun = all_genes_abun[[dataset]][[species]][[GENE]])
                                  },
                                  mc.cores = 40)
    
    raw_num_samples <- vapply(all_out, function(x) { x$num_samples }, integer(1))
    raw_num_alleles <- vapply(all_out, function(x) { x$num_alleles }, integer(1))
    raw_num_strainfacts_strains <- vapply(all_out, function(x) { x$num_strainfacts_strains }, integer(1))
    raw_num_straingst_strains <- vapply(all_out, function(x) { x$num_straingst_strains }, integer(1))
    
    raw_mean_sample_allelic_jaccard <- vapply(all_out, function(x) { x$mean_sample_allelic_jaccard }, numeric(1))
    
    raw_strainfacts_mean_of_min_jaccard <- vapply(all_out, function(x) { x$strainfacts_mean_of_min_jaccard }, numeric(1))
    raw_strainfacts_min_jaccard <- vapply(all_out, function(x) { x$strainfacts_min_jaccard }, numeric(1))
    raw_mean_strainfacts_vs_allele_jaccard <- vapply(all_out, function(x) { x$mean_strainfacts_vs_allele_jaccard }, numeric(1))
    raw_mean_strainfacts_jaccard <- vapply(all_out, function(x) { x$mean_strainfacts_jaccard }, numeric(1))
    raw_mean_permuted_strainfacts_vs_allele_jaccard <- vapply(all_out, function(x) { x$mean_permuted_strainfacts_vs_allele_jaccard }, numeric(1))
    raw_permuted_strainfacts_mean_min_jaccard_P <- vapply(all_out, function(x) { x$permuted_strainfacts_mean_min_jaccard_P }, numeric(1))

    raw_straingst_mean_of_min_jaccard <- vapply(all_out, function(x) { x$straingst_mean_of_min_jaccard }, numeric(1))
    raw_straingst_min_jaccard <- vapply(all_out, function(x) { x$straingst_min_jaccard }, numeric(1))
    raw_mean_straingst_vs_allele_jaccard <- vapply(all_out, function(x) { x$mean_straingst_vs_allele_jaccard }, numeric(1))
    raw_mean_straingst_jaccard <- vapply(all_out, function(x) { x$mean_straingst_jaccard }, numeric(1))
    raw_mean_permuted_straingst_vs_allele_jaccard <- vapply(all_out, function(x) { x$mean_permuted_straingst_vs_allele_jaccard }, numeric(1))
    raw_permuted_straingst_mean_min_jaccard_P <- vapply(all_out, function(x) { x$permuted_straingst_mean_min_jaccard_P }, numeric(1))

    raw_mean_min_list[[paste(species, dataset, sep = '_')]] <- data.frame(species = species,
                                                                          dataset = dataset,
                                                                          gene = names(all_genes_abun[[dataset]][[species]]),
                                                                          num_samples = raw_num_samples,
                                                                          num_alleles = raw_num_alleles,
                                                                          num_strainfacts_strains = raw_num_strainfacts_strains,
                                                                          num_straingst_strains = raw_num_straingst_strains,
                                                                          
                                                                          mean_sample_allelic_jaccard = raw_mean_sample_allelic_jaccard,
                    
                                                                          strainfacts_overall_min = raw_strainfacts_min_jaccard,
                                                                          strainfacts_mean_min = raw_strainfacts_mean_of_min_jaccard,
                                                                          mean_strainfacts_vs_allele_jaccard = raw_mean_strainfacts_vs_allele_jaccard,
                                                                          mean_strainfacts_jaccard = raw_mean_strainfacts_jaccard,
                                                                          mean_permuted_strainfacts_vs_allele_jaccard = raw_mean_permuted_strainfacts_vs_allele_jaccard,
                                                                          permuted_strainfacts_mean_min_jaccard_P = raw_permuted_strainfacts_mean_min_jaccard_P,
                                                                          
                                                                          straingst_overall_min = raw_straingst_min_jaccard,
                                                                          straingst_mean_min = raw_straingst_mean_of_min_jaccard,
                                                                          mean_straingst_vs_allele_jaccard = raw_mean_straingst_vs_allele_jaccard,
                                                                          mean_straingst_jaccard = raw_mean_straingst_jaccard,
                                                                          mean_permuted_straingst_vs_allele_jaccard = raw_mean_permuted_straingst_vs_allele_jaccard,
                                                                          permuted_straingst_mean_min_jaccard_P = raw_permuted_straingst_mean_min_jaccard_P)
    
  }
}

mean_min_out <- do.call(rbind, raw_mean_min_list)
mean_min_out <- mean_min_out[which(! is.na(mean_min_out$strainfacts_overall_min)), ]
mean_min_out <- mean_min_out[which(! is.na(mean_min_out$straingst_overall_min)), ]

write.table(x = mean_min_out,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
