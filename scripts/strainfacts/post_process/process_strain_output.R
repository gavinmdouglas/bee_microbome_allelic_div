# Clean-up StrainFacts output for strain-level inferences.
# Filter out really rare strain abundances (set to 0) and write out new table.

rm(list = ls(all.names = TRUE))

library(reshape2)

filt_cutoff <- 0.01

datasets = c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

strain_comm <- list()

for (dataset in datasets) {
  
  output_dir <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output', dataset, sep = '/')
  
  strain_comm_files <- list.files(path = output_dir,
                                  full.names = TRUE, pattern = ".comm.tsv.gz")

  for (strain_comm_file in strain_comm_files) {

    species <- basename(strain_comm_file)

    species <- gsub(".comm.tsv.gz", "", species)

    if (! species %in% names(strain_comm)) {
      strain_comm[[species]] <- list()
    }
        
    sample_list_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_core_input',
                              dataset,
                              species,
                             'samples.tsv', sep = '/')
    
    sample_map <- read.table(file = sample_list_file, stringsAsFactors = FALSE, sep = "\t", header = FALSE, row.names = 1)
    
    comm_output <- read.table(file = strain_comm_file, header = TRUE, sep = "\t")
    
    if (length(which(comm_output$community < filt_cutoff)) > 0) {
      comm_output[which(comm_output$community < filt_cutoff), "community"] <- 0
    }
    
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
    
    print(c(dataset, species, min(colSums(comm_output_wide))))

    comm_output_wide <- data.frame(sweep(comm_output_wide, 2, colSums(comm_output_wide), '/')) * 100
  
    write.table(file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/community_relabun/',
                             dataset, '.', species, '.tsv', sep = ''),
                x = comm_output_wide,
                col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)

    strain_comm[[species]][[dataset]] <- comm_output_wide
  
  }

}

saveRDS(object = strain_comm,
        file = '/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/RDS/strainfacts_core.genome_comm.rds')
