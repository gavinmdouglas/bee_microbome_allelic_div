rm(list = ls(all.names = TRUE))

library(parallel)
library(vegan)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

# Parse StrainFacts (subsamples / replicate) output and determine centroid replicate per gene/dataset.
# Write this information to a file, so that just these centroid replicates
# will be parsed when parsing all of the output files.

comm_file_path <- '/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_out_subsampled/subsample20/comm'
prepped_path <- '/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input_subsampled20/'

all_comm_files <- list.files(path = comm_file_path,
                            pattern = ".tsv.gz$",
                            full.names = TRUE)

# Get single file per gene (collapsing all replicates).
all_comm_files_rep_rm <- sub('\\.rep.\\.', '.rep.*.', all_comm_files)
all_comm_files_rep_rm <- sub('\\.rep10\\.', '.rep.*.', all_comm_files_rep_rm)
unique_all_comm_genes <- all_comm_files_rep_rm[which(! duplicated(all_comm_files_rep_rm))]

centroid_replicate_files <- parallel::mclapply(unique_all_comm_genes, identify_centroid_comm_replicate, mc.cores = 40)

centroid_replicate_files <- unlist(centroid_replicate_files)
centroid_replicate_files <- centroid_replicate_files[which(! is.na(centroid_replicate_files))]

write.table(x = centroid_replicate_files,
            file = '/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_out_subsampled_processed/subsample20/centroid_rep_comm_files.txt',
            col.names = FALSE, row.names = FALSE, quote = FALSE)
