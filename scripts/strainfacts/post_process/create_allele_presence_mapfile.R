rm(list = ls(all.names = TRUE))

# Create simple file mapping gene ids to the allele ids that were called as present (based on relative abundance) across the samples.
all_genes_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')


for(d in names(all_genes_abun)) {

  print(d)

  raw_allele_ids <- list()

  for (species in names(all_genes_abun[[d]])) {

    focal_alleles <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory/', species, '.txt', sep = ''),
                                stringsAsFactors = FALSE, header = FALSE)$V1
    
    print(species)

    tmp_abun <- all_genes_abun[[d]][[species]]

    intersecting_genes <- intersect(names(tmp_abun), focal_alleles)
    
    tmp_abun <- tmp_abun[intersecting_genes]
    
    raw_allele_ids[[species]] <- data.frame(gene = names(tmp_abun),
                                            allele_numbers_present = NA)

    for (i in 1:length(tmp_abun)) {

      raw_allele_ids[[species]][i, 'allele_numbers_present'] <- paste(sort(unique(tmp_abun[[i]]$strain)),
                                                                      collapse = ',')

    }

  }

  allele_ids <- do.call(rbind, raw_allele_ids)
  rownames(allele_ids) <- NULL
  
  outfile <- paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/alleles_called_present_focal/', d, '.tsv', sep = '')
  
  write.table(x = allele_ids,
              file = outfile,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}

