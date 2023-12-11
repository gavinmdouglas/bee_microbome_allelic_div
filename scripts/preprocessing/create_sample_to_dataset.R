rm(list = ls(all.names = TRUE))

dataset_to_samples <- list()

datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

for (dataset in datasets) {

  SRR_file <- paste("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/", dataset, "_SRRs.txt.gz", sep = "")
  
  dataset_to_samples[[dataset]] <- data.frame(sample = read.table(SRR_file, stringsAsFactors = FALSE, header = FALSE)$V1,
                                              dataset = dataset)
  
  
    
}

dataset_map <- do.call(rbind, dataset_to_samples)

write.table(x = dataset_map, file = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/sample_to_dataset.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
