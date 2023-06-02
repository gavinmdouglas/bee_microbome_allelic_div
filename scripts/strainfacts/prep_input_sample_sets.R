# Get sample sets where species called as present based on having at least 90 percent core.

rm(list = ls(all.names = TRUE))

species_presence_90per <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

for (sp in colnames(species_presence_90per)) {

  samples_present <- sort(rownames(species_presence_90per)[which(species_presence_90per[, sp] > 0)], decreasing = FALSE)

  if (length(samples_present) <= 2) { next }

  outfile <- paste("/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/strainfacts_running/core/pre_info/samples_w_90percent_core/",
                   sp,
                   ".txt", sep = "")

  write.table(x = samples_present, file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)

}
