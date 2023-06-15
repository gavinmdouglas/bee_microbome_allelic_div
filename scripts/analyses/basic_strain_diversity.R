rm(list = ls(all.names = TRUE))

# Per-species, get basic diversity measures per species.
# Already have the mean number of strains per sample (i.e., richness),
# but also just get the overall number of strains in total, as well as the
# mean evenness and Shannon index.

# For strains.
strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

species_strain_diversity <- data.frame(species = names(strain_relabun),
                                       total_samples_w_strains = NA,
                                       total_strains = NA,
                                       mean_strain_SimpsonsE = NA,
                                       sd_strain_SimpsonsE = NA,
                                       mean_strain_Shannon = NA,
                                       sd_strain_Shannon = NA)
rownames(species_strain_diversity) <- names(strain_relabun)

for (sp in names(strain_relabun)) {
   strain_profile <- strain_relabun[[sp]] / colSums(strain_relabun[[sp]])
   
   simpsons_e <- (1 - colSums(strain_profile ** 2)) / colSums(strain_profile > 0)
   shannon <- vegan::diversity(strain_relabun[[sp]], index = 'shannon', MARGIN = 2, base = exp(1))
   
   species_strain_diversity[sp, 'species'] <- sp
   species_strain_diversity[sp, 'total_samples_w_strains'] <- ncol(strain_profile)
   species_strain_diversity[sp, 'total_strains'] <- nrow(strain_profile)
   
   species_strain_diversity[sp, 'mean_strain_SimpsonsE'] <- mean(simpsons_e)
   species_strain_diversity[sp, 'sd_strain_SimpsonsE'] <- sd(simpsons_e)
   
   species_strain_diversity[sp, 'mean_strain_Shannon'] <- mean(shannon)
   species_strain_diversity[sp, 'sd_strain_Shannon'] <- sd(shannon)
}

write.table(x = species_strain_diversity,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_basic_diversity.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
