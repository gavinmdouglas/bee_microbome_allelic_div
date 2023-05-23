# Filter CheckM output table to identify genome accession that pass quality thresholds.
# Write these retained accessions out by species for further downstream steps.

rm(list = ls(all.names = TRUE))

checkm_output <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/CheckM_output.tsv.gz',
                            header = TRUE, sep ="\t", stringsAsFactors = FALSE, row.names = 1)

# Read in accession ids by species
accessions_by_species <- list()
species_names <- gsub(".tsv$", "", list.files("../adding_new_microbiota_genomes/accessions_to_process/"))
accession_filepaths <- list.files("../adding_new_microbiota_genomes/accessions_to_process/", full.names = TRUE)
for (i in 1:length(accession_filepaths)) {
  accessions_by_species[[species_names[i]]] <- gsub("GCF_", "GCA_", read.table(accession_filepaths[i],
                                                          stringsAsFactors = FALSE, sep = "\t", header = TRUE)$accession)
}

checkm_output_highqual <- checkm_output[which(checkm_output$Completeness > 95 & checkm_output$Contamination < 2), ]

# Write out these accession files (note that these were intermediate files and not actually the *final* genome sets)
# per species. This is because at the next stage, FastANI pairwise comparisons were used to correct many of the
# species labels (and six genomes were removed). 
for (s in species_names) {
  species_accessions <- accessions_by_species[[s]]
  
  # Re-format species genome accessions to match "GCA.GCF" format.
  species_accessions <- gsub('^GC._', 'GCA.GCF_', species_accessions)

  highqual_species_accessions <- species_accessions[which(species_accessions %in% rownames(checkm_output_highqual))]

  outfile <- paste0('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_draft/', s, '.txt')

  write.table(x = checkm_output_highqual[highqual_species_accessions, 'accession'],
              file = outfile,
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

