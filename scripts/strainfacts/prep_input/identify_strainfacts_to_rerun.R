rm(list = ls(all.names = TRUE))

# Identify gene allele outputs that differ from re-run input after removing INDELs,
# even though these genes don't have INDELs in the final set, apparently including INDELs when filtering by samples
# changed the final set.
# Also check whether the input sites are exactly the same or not, which could also theoretically change.

# Ignore genes that have already been re-run because they had INDELs in the final output.

# Note that this script would not need to be used to reproduce the analyses, as the python script will now ignore INDELs.

all_species <- colnames(read.table('/scratch/gdouglas/projects/honey_bee/species_presence_core_90percent.tsv.gz',
                          header = TRUE, stringsAsFactors = FALSE, row.names = 1))

datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

allele_relabun <- list()

# Only consider focal genes.
focal_genes <- character()
for (GENE_FILE in list.files("/scratch/gdouglas/projects/honey_bee/gene_sets/accessory/", full.names = TRUE)) {
  tmp_genes <- read.table(GENE_FILE, stringsAsFactors = FALSE, header=FALSE)$V1
  focal_genes <- c(focal_genes, tmp_genes)
}

# Genes already rerun
already_rerun <- list()
for (d in datasets) {
  RERUN_GENE_FILE = paste("/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/accessory_genes_to_rerun/", d, ".txt", sep = "")
  already_rerun[[d]] <- read.table(RERUN_GENE_FILE, stringsAsFactors = FALSE, header=FALSE)$V1
}

to_rerun <- list()

for (d in datasets) {
  
  to_rerun[[d]] <- character()
  
  print(d)

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
  
  for (gene in all_genes) {
    
    if (gene %in% already_rerun[[d]]) { next }
    if (! gene %in% focal_genes) { stop("Non-focal gene!") }
    
    new_sample_tab_file <- paste(dataset_samples_infolder, gene, '_samples.tsv', sep = "")
    new_sample_tab <- read.table(new_sample_tab_file, header = FALSE, sep = '\t')
    
    new_site_tab_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input/',
                               d, '/sites/', gene, "_sites.tsv", sep = '')
    new_site_tab <- read.table(file = new_site_tab_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
     
    old_sample_tab_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input_OLD/prepped_accessory_input/',
                                 d, '/samples/', gene, '_samples.tsv', sep = "")
    old_sample_tab <- read.table(old_sample_tab_file, header = FALSE, sep = '\t')
    
    old_site_tab_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input_OLD/prepped_accessory_input/',
                               d, '/sites/', gene, "_sites.tsv", sep = '')
    old_site_tab <- read.table(file = old_site_tab_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
    
    if (! identical(new_site_tab, old_site_tab)) {
      to_rerun[[d]] <- c(to_rerun[[d]], gene)
    } else if (! identical(new_sample_tab, old_sample_tab)) {
      to_rerun[[d]] <- c(to_rerun[[d]], gene)
    }
    
  }
  
}

for (d in names(to_rerun)) {
  
  outfile <- paste("/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/additional_accessory_genes_to_rerun/", d, ".txt", sep = "")

  write.table(x = to_rerun[[d]], file = outfile, col.names = FALSE, row.names = FALSE, quote = FALSE)  
}

