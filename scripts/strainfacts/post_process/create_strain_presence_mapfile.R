rm(list = ls(all.names = TRUE))

# Create simple file mapping gene ids to the allele ids that were called as present (based on relative abundance) across the samples.

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

strain_ids_present <- data.frame(species = names(strain_relabun),
                                 strain_numbers_present = NA)
rownames(strain_ids_present) <- names(strain_relabun)

for (species in names(strain_relabun)) {
  strain_ids_present[species, 'strain_numbers_present'] <- paste(gsub('^strain_', '', rownames(strain_relabun[[species]])), collapse = ',')
}

write.table(x = strain_ids_present,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strains_called_present.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
