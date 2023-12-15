rm(list = ls(all.names = TRUE))

# Create gene allele rel. abun. tables converted to lists of RDS objects that are quicker to read in and use.
library(parallel)

preprocess_comm <- function(gene, dataset) {

  gene_samples_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input/',
                             dataset,
                             '/samples/',
                             gene,
                             '_samples.tsv',
                             sep = '')
  
  gene_samples <- read.table(gene_samples_file, sep = "\t", row.names = 1, header = FALSE)
  
  comms_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_output/',
                      dataset,
                      '/comm/',
                      gene,
                      '.comm.tsv.gz',
                      sep = '')

  allele_freq <- read.table(comms_file, header = TRUE, sep = "\t")
  
  if (length(unique(allele_freq$sample)) != nrow(gene_samples)) {
    stop('Mismatch in expected number of samples and those in sample mapfile!') 
  }
  
  allele_freq$sample <- as.character(allele_freq$sample)
  
  if (length(which(! allele_freq$sample %in% rownames(gene_samples))) > 0) {
    stop('Error - sample missing in mapfile.') 
  }
  
  allele_freq$sample <- gene_samples[allele_freq$sample, 1]
  
  orig_samples <- sort(unique(allele_freq$sample))
  
  # Only keep strain calls that were at least >= 1%
  allele_freq <- allele_freq[which(allele_freq$community >= 0.01), ]
  
  if (! identical(orig_samples, sort(unique(allele_freq$sample)))) {
     stop("ERROR - samples dropped after excluding low abun strains.")
  }
  
  # Total-sum scale by sample remaining rel. abun.
  sample_sums <- aggregate(x = community ~ sample,
                           data = allele_freq,
                           FUN = sum)
  rownames(sample_sums) <- sample_sums$sample
  allele_freq$community <- (allele_freq$community / sample_sums[allele_freq$sample, "community"]) * 100
  
  return(allele_freq)

}

all_species <- colnames(read.table('/scratch/gdouglas/projects/honey_bee/species_presence_core_90percent.tsv.gz',
                          header = TRUE, stringsAsFactors = FALSE, row.names = 1))

datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

allele_relabun <- list()

for (d in datasets) {
  
  print(d)

  allele_relabun[[d]] <- list()
  
  dataset_fit_outfolder <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_output/',
                                 d,
                                 '/fit/', sep = '')
  
  dataset_samples_infolder <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input/',
                                    d,
                                    '/samples/', sep = '')
  
  all_fit <- list.files(dataset_fit_outfolder, pattern = ".fit$")
  out_genes <- gsub(".fit$", "", all_fit)
  
  all_samples <- list.files(dataset_samples_infolder, pattern = "_samples.tsv$")
  in_genes <- gsub("_samples.tsv$", "", all_samples)
  
  all_genes <- intersect(out_genes, in_genes)
  
  for (sp in all_species) {
    
    print(sp)
    
    sp_genes <- grep(sp, all_genes, value = TRUE)
    
    if (length(sp_genes) == 0) { next }
    
    allele_relabun[[d]][[sp]] <- parallel::mclapply(X = sp_genes, FUN = preprocess_comm, dataset = d, mc.cores = 40)

    names(allele_relabun[[d]][[sp]]) <- sp_genes
  
  }

}

saveRDS(object = allele_relabun,
        file = "/scratch/gdouglas/projects/honey_bee/strainfacts_working/accessory_output_processed/strainfacts_accessory_allele_relabun.rds")
