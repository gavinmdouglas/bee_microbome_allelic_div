rm(list = ls(all.names = TRUE))

# Create simple file mapping gene ids to the allele ids that were called as present (based on relative abundance) across the samples.

all_genes_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

raw_allele_ids <- list()

for (species in names(all_genes_abun)) {
  
  print(species)
  
  raw_allele_ids[[species]] <- data.frame(gene = names(all_genes_abun[[species]]),
                                          allele_numbers_present = NA)
  
  for (i in 1:length(all_genes_abun[[species]])) {

    raw_allele_ids[[species]][i, 'allele_numbers_present'] <- paste(sort(unique(all_genes_abun[[species]][[i]]$strain)),
                                                                    collapse = ',')

  }

}

allele_ids <- do.call(rbind, raw_allele_ids)
rownames(allele_ids) <- NULL

write.table(x = allele_ids,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/alleles_called_present.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
