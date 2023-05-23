rm(list = ls(all.names = TRUE))

# Quick reference of the number and percentage of genomes filtered out at CheckM step.

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
checkm_filt_tallies <- data.frame(matrix(NA, nrow = length(species_names), ncol = 2))
rownames(checkm_filt_tallies) <- species_names
colnames(checkm_filt_tallies) <- c("Downloaded", "Pass_Complete95_Contam2")


for (s in species_names) {
  species_accessions <- accessions_by_species[[s]]
  
  # Re-format species genome accessions to match "GCA.GCF" format.
  species_accessions <- gsub('^GC._', 'GCA.GCF_', species_accessions)
  
  checkm_filt_tallies[s, "Downloaded"] <- length(species_accessions)
  checkm_filt_tallies[s, "Pass_Complete95_Contam2"] <- length(which(species_accessions %in% rownames(checkm_output_highqual)))
}

# Total number of downloaded genomes.
sum(checkm_filt_tallies$Downloaded)

# Number of CheckM-passed genomes (and percentage).
sum(checkm_filt_tallies$Pass_Complete95_Contam2)

round((sum(checkm_filt_tallies$Pass_Complete95_Contam2) / sum(checkm_filt_tallies$Downloaded)) * 100, 2)
