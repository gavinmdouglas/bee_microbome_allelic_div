rm(list = ls(all.names = TRUE))

# Create simple file mapping gene ids to the allele ids that were called as present (based on relative abundance) across the samples.

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

strain_ids_present <- data.frame(species = names(strain_relabun),
                                 Ellegaard2019 = NA,
                                 Ellegaard2020 = NA,
                                 Sun2022 = NA,
                                 Wu2021 = NA,
                                 Zhang2022 = NA)
rownames(strain_ids_present) <- names(strain_relabun)

for (species in names(strain_relabun)) {
  for (d in names(strain_relabun[[species]])) {
    strain_ids_present[species, d] <- paste(gsub('^strain_', '', rownames(strain_relabun[[species]][[d]])), collapse = ',')
  }
}

write.table(x = strain_ids_present,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strains_called_present.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
