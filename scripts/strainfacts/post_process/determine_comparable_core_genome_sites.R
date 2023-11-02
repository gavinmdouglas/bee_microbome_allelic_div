rm(list = ls(all.names = TRUE))

# Determine sites in core genome alignment(s) that should be compared.
# All other sites should be converted to N's and/or removed.

all_species <- read.table('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/species_w_strains.txt',
                          header = FALSE, sep = ' ', stringsAsFactors = FALSE)$V1

for (sp in all_species) {
 
  included_datasets <- read.table(paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                                  sp,
                                  '.datasets.intersect.txt',
                                  sep = ''),
                                  stringsAsFactors = FALSE, header = FALSE)$V1

  all_possible_datasets <- read.table(paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
					                                      sp,
									                                        '.datasets.txt',
									                                        sep = ''),
				                                        stringsAsFactors = FALSE, header = FALSE)$V1

  dataset_sites <- list()
  for (d in all_possible_datasets) {

    dataset_folder <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_core_input/',
                            d,
                            '/',
                            sp,
                            sep = '')

	  dataset_any_genes <- read.table(paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                               sp,
			       '.',
			       d,
                               '.genes.txt',
                               sep = ''),
                               stringsAsFactors = FALSE, header = FALSE)$V1


    invariant_site_file <- paste(dataset_folder, 'passable_invariant_sites.tsv', sep = '/')
    variant_site_file <- paste(dataset_folder, 'sites.tsv', sep = '/')

    invariant_sites <- read.table(file = invariant_site_file,
                                  header = FALSE,
                                  sep = '\t',
                                  stringsAsFactors = FALSE)
    colnames(invariant_sites) <- c('gene', 'pos')

    variant_sites_raw <- read.table(file = variant_site_file,
                                header = FALSE,
                                sep = '\t',
                                stringsAsFactors = FALSE)

    variant_sites_raw <- strsplit(x = variant_sites_raw$V2, split = '\\|')

    variant_sites <- data.frame(gene = sapply(variant_sites_raw, function(x) { x[1] }),
                                pos = sapply(variant_sites_raw, function(x) { x[2] }))
    
    invariant_sites <- invariant_sites[which(invariant_sites$gene %in% dataset_any_genes), ]
    variant_sites <- variant_sites[which(variant_sites$gene %in% dataset_any_genes), ]
    
    invariant_sites_pasted <- paste(invariant_sites$gene, invariant_sites$pos, sep = '|')
    variant_sites_pasted <- paste(variant_sites$gene, variant_sites$pos, sep = '|')

    dataset_sites[[d]] <- c(invariant_sites_pasted, variant_sites_pasted)
    
    dataset_sites_split <- strsplit(x = dataset_sites[[d]], split = '\\|')
    
    dataset_sites_df <- data.frame(gene = sapply(dataset_sites_split, function(x) { x[1] }),
                                   pos = sapply(dataset_sites_split, function(x) { x[2] }))
    
    write.table(x = dataset_sites_df,
                file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                             sp,
                             '.',
                             d,
                             '.sites.txt', sep = ''),
                sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    
  }
  
  intersecting_sites <- dataset_sites[[included_datasets[1]]]
  if (length(included_datasets) > 1) {
    for (d in included_datasets[-1]) {
      intersecting_sites <- intersect(intersecting_sites, dataset_sites[[d]])
    }
  }
  
  intersecting_sites <- sort(intersecting_sites)
  
  intersecting_sites_split <- strsplit(x = intersecting_sites, split = '\\|')

  intersecting_sites_df <- data.frame(gene = sapply(intersecting_sites_split, function(x) { x[1] }),
                                      pos = sapply(intersecting_sites_split, function(x) { x[2] }))

  write.table(x = intersecting_sites_df,
              file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                           sp,
                           '.sites.intersect.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

}
