rm(list = ls(all.names = TRUE))

# Per-species, get basic diversity measures per species.
# Already have the mean number of strains per sample (i.e., richness),
# but also just get the overall number of strains in total, as well as the
# mean evenness and Shannon index.

# Also keep track of mean number of strains per sample for StrainGST profiles (restricted to the same samples).

# For strains.
strainfacts_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
straingst_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')

# Make sure the same samples are compared.
strain_species <- intersect(names(strainfacts_relabun), names(straingst_relabun))
strain_samples <- list()

for (sp in strain_species) {
  tmp_datasets <- intersect(names(strainfacts_relabun[[sp]]), names(straingst_relabun[[sp]]))
  if (length(tmp_datasets) == 0) { next }
  
  for (tmp_d in tmp_datasets) {
    tmp_intersecting_samp <- intersect(colnames(strainfacts_relabun[[sp]][[tmp_d]]), colnames(straingst_relabun[[sp]][[tmp_d]]))
    if (length(tmp_intersecting_samp) == 0) { next }
    
    if (! sp %in% names(strain_samples))   {
      strain_samples[[sp]] <- list()
    }
    
    strain_samples[[sp]][[tmp_d]] <- tmp_intersecting_samp
  }
}

summary_raw <- list()

for (sp in names(strain_samples)) {
  
  for (d in names(strain_samples[[sp]])) {
    
    num_strainfacts_samples_w_strains <- ncol(strainfacts_relabun[[sp]][[d]])
    num_straingst_samples_w_strains <- ncol(straingst_relabun[[sp]][[d]])
    num_intersecting_samples <- length(strain_samples[[sp]][[d]])
   
    straingst_profile <- straingst_relabun[[sp]][[d]][, strain_samples[[sp]][[d]], drop = FALSE]
    straingst_profile <- straingst_profile[which(rowSums(straingst_profile) > 0), , drop = FALSE]
    straingst_profile <- straingst_profile / colSums(straingst_profile)
    
    strainfacts_profile <- strainfacts_relabun[[sp]][[d]][, strain_samples[[sp]][[d]], drop = FALSE]
    strainfacts_profile <- strainfacts_profile[which(rowSums(strainfacts_profile) > 0), , drop = FALSE]
    strainfacts_profile <- strainfacts_profile / colSums(strainfacts_profile)
   
    straingst_richness <- colSums(straingst_profile > 0)
    straingst_simpsons_e <- (1 - colSums(straingst_profile ** 2)) / colSums(straingst_profile > 0)
    straingst_shannon <- vegan::diversity(straingst_profile, index = 'shannon', MARGIN = 2, base = exp(1))

    strainfacts_richness <- colSums(strainfacts_profile > 0)
    strainfacts_simpsons_e <- (1 - colSums(strainfacts_profile ** 2)) / colSums(strainfacts_profile > 0)
    strainfacts_shannon <- vegan::diversity(strainfacts_profile, index = 'shannon', MARGIN = 2, base = exp(1))

   summary_raw[[paste(sp, d, sep = '_')]] <- data.frame(species = sp,
                                                        dataset = d,
                                                        
                                                        num_straingst_samples_w_strains = num_straingst_samples_w_strains,
                                                        num_strainfacts_samples_w_strains = num_strainfacts_samples_w_strains,
                                                        num_intersecting_samples = num_intersecting_samples,
                                                        
                                                        straingst_total_strains = nrow(straingst_profile),
                                                        strainfacts_total_strains = nrow(strainfacts_profile),
                                                        
                                                        straingst_mean_richness = mean(straingst_richness),
                                                        straingst_sd_richness = sd(straingst_richness),
                                                        mean_straingst_simpsons_e = mean(straingst_simpsons_e),
                                                        sd_straingst_simpsons_e = sd(straingst_simpsons_e),
                                                        mean_straingst_shannon = mean(straingst_shannon),
                                                        sd_straingst_shannon = sd(straingst_shannon),
                                                        
                                                        strainfacts_mean_richness = mean(strainfacts_richness),
                                                        strainfacts_sd_richness = sd(strainfacts_richness),
                                                        mean_strainfacts_simpsons_e = mean(strainfacts_simpsons_e),
                                                        sd_strainfacts_simpsons_e = sd(strainfacts_simpsons_e),
                                                        mean_strainfacts_shannon = mean(strainfacts_shannon),
                                                        sd_strainfacts_shannon = sd(strainfacts_shannon))

  }
}

strain_summary <- do.call(rbind, summary_raw)

write.table(x = strain_summary,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_diversity.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
