# Clean-up strainfacts output for strains.
# Filter out really rare strain abundances (set to 0) and write out new table.

rm(list = ls(all.names = TRUE))

library(reshape2)

filt_cutoff <- 0.01

strain_comm <- list()

strain_comm_files <- list.files(path = "/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/output/core/comm/",
                                  full.names = TRUE, pattern = ".comm.tsv.gz")

for (strain_comm_file in strain_comm_files) {
  
  species <- basename(strain_comm_file)
  
  species <- gsub(".comm.tsv.gz", "", species)
  
  print(species)
  
  sample_list_file <- paste("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/prepped_input/",
                            species, "/samples.tsv.gz", sep = "")
  
  sample_map <- read.table(file = sample_list_file, stringsAsFactors = FALSE, sep = "\t", header = FALSE, row.names = 1)
  
  comm_output <- read.table(file = strain_comm_file, header = TRUE, sep = "\t")
  
  comm_output[which(comm_output$community < filt_cutoff), "community"] <- 0
  
  if (length(unique(comm_output$sample)) != nrow(sample_map)) {
    stop('Mismatch in expected number of samples and those in sample mapfile!') 
  }
  
  comm_output$sample <- as.character(comm_output$sample)

  if (length(which(! comm_output$sample %in% rownames(sample_map))) > 0) {
    stop('Error - sample missing in mapfile.') 
  }
  
  comm_output$sample <- sample_map[comm_output$sample, 1]
  
  comm_output_wide <- dcast(data = comm_output, formula = strain ~ sample, value.var = "community")
  rownames(comm_output_wide) <- paste("strain", as.character(comm_output_wide$strain), sep = "_")
  comm_output_wide <- comm_output_wide[, -1]
  
  comm_output_wide <- comm_output_wide[which(rowSums(comm_output_wide) > 0), , drop = FALSE]
  
  print(min(colSums(comm_output_wide)))

  comm_output_wide <- data.frame(sweep(comm_output_wide, 2, colSums(comm_output_wide), '/')) * 100

  write.table(file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                        species, '.tsv', sep = ''),
              x = comm_output_wide,
              col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    
  strain_comm[[species]] <- comm_output_wide

}

saveRDS(object = strain_comm,
        file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
