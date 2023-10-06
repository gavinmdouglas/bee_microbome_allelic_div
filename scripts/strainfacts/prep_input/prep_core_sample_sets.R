# Get sample sets where species called as present based on having at least 90 percent core.

rm(list = ls(all.names = TRUE))

species_presence_90per <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dataset_samples <- list()
datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')
for (dataset in datasets) {
  SRR_path <- paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/', dataset, '_SRRs.txt.gz', sep = '')
  dataset_samples[[dataset]] <- read.table(SRR_path, header = FALSE, stringsAsFactors = FALSE)$V1
}


for (sp in colnames(species_presence_90per)) {

  samples_present <- sort(rownames(species_presence_90per)[which(species_presence_90per[, sp] > 0)], decreasing = FALSE)

  if (length(samples_present) <= 2) { next }

  outfile <- paste("/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/core_pre_info/samples_w_90percent_core_TEST/",
                   sp,
                   ".txt", sep = "")

  write.table(x = samples_present, file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Also get sample set for each dataset individually.
  for (dataset in datasets) {
    dataset_samples_present <- intersect(samples_present, dataset_samples[[dataset]])
    if (length(dataset_samples_present) <= 2) { next }
    dataset_outfile <- paste("/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/core_pre_info/samples_w_90percent_core_by_dataset/",
                             dataset,
                             '/',
                             sp,
                             ".txt", sep = "")
    write.table(x = dataset_samples_present, file = dataset_outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
}
