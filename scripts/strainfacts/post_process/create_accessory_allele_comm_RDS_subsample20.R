rm(list = ls(all.names = TRUE))

# Create gene allele rel. abun. tables converted to lists of RDS objects that are quicker to read in and use.
library(parallel)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

all_species <- colnames(read.table('/scratch/gdouglas/projects/honey_bee/species_presence_core_90percent.tsv.gz',
                                   header = TRUE, stringsAsFactors = FALSE, row.names = 1))

# Initialize output list.
allele_relabun <- list()

datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

# Read in previously-identified centroid outputs per replicate.
centroid_comms <- read.table(file = '/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_out_subsampled_processed/subsample20/centroid_rep_comm_files.txt',
                             stringsAsFactors = FALSE, header = FALSE)$V1

for (d in datasets) {

  print(d)
  
  allele_relabun[[d]] <- list()
  
  for (sp in all_species) {

    print(sp)
    
    centroid_comms_subset <- grep(d, centroid_comms, value = TRUE)
    centroid_comms_subset <- grep(sp, centroid_comms_subset, value = TRUE)
    
    if (length(centroid_comms_subset) == 0) { next }
    
    allele_relabun[[d]][[sp]] <- parallel::mclapply(X = centroid_comms_subset, FUN = preprocess_subsampled_comm, mc.cores = 40)
    
    genes_subset <- sapply(centroid_comms_subset,
                           function(x) { 
                             strsplit(x, '\\.')[[1]][4]
                           }
    )
    names(genes_subset) <- NULL
    genes_subset <- gsub('_metagenotype$', '', genes_subset)

    names(allele_relabun[[d]][[sp]]) <- genes_subset
    
  }
}

saveRDS(object = allele_relabun,
        file = "/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_out_subsampled_processed/subsample20/strainfacts_accessory_allele_relabun.rds")
