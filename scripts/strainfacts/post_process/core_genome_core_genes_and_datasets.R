rm(list = ls(all.names = TRUE))

# Determine core genes that should be used for core genome alignment per species.
# Key step is to throw out datasets that cause the number of genes to include to be below 50.

# Also, realized that just focusing on dataset-specific comparisons is very useful as well,
# so I also output the genes per dataset as well.

all_species <- read.table('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/species_w_strains.txt',
                          header = FALSE, sep = ' ', stringsAsFactors = FALSE)$V1

for (sp in all_species) {
 
  ref_genes <- gsub('.fa',
                    '',
                    list.files(paste('/scratch/gdouglas/projects/honey_bee/ref_genomes/core_gene_fastas',
                                     sp,
                                     sep = '/')))
  
  sp_mapfiles <- list.files('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/core_genome_mapfiles',
                            full.names = TRUE,
                            pattern = sp)
  
  sp_datasets <- sort(sapply(sp_mapfiles,
                        function(x) { gsub('\\..*$', '', basename(x)) }))
  names(sp_datasets) <- NULL
  write.table(x = sp_datasets,
              file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                           sp,
                           '.datasets.txt',
                           sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  max_num_intersects <- 9999999
  
  included_datasets <- character()
  gene_set <- ref_genes
  
  # First get union of all dataset genes that interesect with references,
  # and write out individual dataset gene sets (again that intersect with ref).
  sp_any_dataset_gene <- character()
  for (sp_mapfile in sp_mapfiles) {
    
    dataset <- basename(sp_mapfile)
    dataset <- gsub('\\..*$', '', dataset)

    dataset_genes <- read.table(file = sp_mapfile,
                                stringsAsFactors = FALSE,
                                header = TRUE)$gene
    
    dataset_genes_ref_intersect <- intersect(dataset_genes, ref_genes)
    
    sp_any_dataset_gene <- unique(c(sp_any_dataset_gene, dataset_genes_ref_intersect))
    
    write.table(x = dataset_genes_ref_intersect,
                file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                             sp,
                             '.',
                             dataset,
                             '.genes.txt', sep = ''),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  write.table(x = sp_any_dataset_gene,
              file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                           sp, '.genes.any.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
  while(length(sp_mapfiles) > 0) {

    dataset_intersects <- integer()
    datasets <- character()

    for (sp_mapfile in sp_mapfiles) {
      
      dataset <- basename(sp_mapfile)
      
      dataset <- gsub('\\..*$', '', dataset)
      
      datasets <- c(datasets, dataset)

      dataset_genes <- read.table(file = sp_mapfile,
                                  stringsAsFactors = FALSE,
                                  header = TRUE)$gene
      
      dataset_intersects <- c(dataset_intersects,
                              length(intersect(gene_set,
                                               dataset_genes)))

    }

    max_num_intersects <- max(dataset_intersects)
    if (max_num_intersects < 50) { break }

    max_intersect_i <- which.max(dataset_intersects)

    included_datasets <- c(included_datasets, datasets[max_intersect_i])

    best_dataset_genes <- read.table(file = sp_mapfiles[max_intersect_i],
                                     stringsAsFactors = FALSE,
                                     header = TRUE)$gene
    
    gene_set <- intersect(gene_set, best_dataset_genes)
    
    sp_mapfiles <- sp_mapfiles[-max_intersect_i]
    
  }

  write.table(x = gene_set,
              file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                           sp,
                           '.genes.intersect.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  write.table(x = included_datasets,
              file = paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/core_output_processed/strainfacts_vs_ref_core_gene_intersect/',
                           sp,
                           '.datasets.intersect.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

}
