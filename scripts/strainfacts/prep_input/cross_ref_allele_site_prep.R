rm(list = ls(all.names = TRUE))

# Realized that INDELs had been included in StrainFacts input, but that this was not desired afterwards.
# This script parses the input sites files output at the time of prepping the the StrainFacts input,
# and identifies and genes that included INDELs in the variants being used.

# Note that only filtered accessory genes are considered, to reduce the number of re-run jobs.

# Note that this script would not be needed to reproduce results (as the scripts to prep strainfacts input now ignore INDELs)
filt_genes <- character()

gene_id_files <- list.files(path = '/scratch/gdouglas/projects/honey_bee/gene_sets/accessory',
                            pattern = '.txt.gz', full.names = TRUE)

for (gene_id_file in gene_id_files) {
  filt_genes <- c(filt_genes,
                  read.table(file = gene_id_file, header = FALSE, stringsAsFactors = FALSE)$V1)
}

datasets = c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

for (d in datasets) {
  
  genes_to_reun <- character()
  
  dataset_site_files <-  list.files(path = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input', d, 'sites', sep = '/'),
                                    pattern = '.tsv.gz', full.names = TRUE)

  for (dataset_site_file in dataset_site_files) {
    
    site_gene <- gsub("_sites.tsv.gz$", "", basename(dataset_site_file))
    
    if (! site_gene %in% filt_genes) { next }
    
    sites_tab <- read.table(file = dataset_site_file, header = FALSE, sep = '_', stringsAsFactors = FALSE)
    
    base_nchar <- c(sapply(sites_tab$V2, nchar), sapply(sites_tab$V3, nchar))
    
    if (length(which(base_nchar != 1)) > 0) {
      genes_to_reun <- c(genes_to_reun, site_gene)
    }
    
  }
    
  gene_outfile <- paste("/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/accessory_genes_to_rerun/",
                        d, ".txt", sep = "")

  write.table(x = genes_to_reun, file = gene_outfile, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
