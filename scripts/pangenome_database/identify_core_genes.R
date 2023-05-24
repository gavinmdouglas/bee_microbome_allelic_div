### Call core genes for each species and make a single table with all genes mapped to.
### For species with < 10 genomes, choose core genes based on CheckM lineage core genes.

rm(list = ls(all.names = TRUE))

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt",
                      stringsAsFactors = FALSE)$V1

custom_core_needed <- c()

core_genes <- list()

for (sp in species) {

  sp_panaroo_outfile <- paste("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_genome_panaroo/", sp, "/gene_presence_absence_roary.csv.gz", sep = "")
  
  sp_panaroo <- read.table(sp_panaroo_outfile, header = TRUE, sep = ",", comment.char = "", quote = "")
  
  rownames(sp_panaroo) <- paste(sp, sp_panaroo$Gene, sep = "_")
  
  sp_num_genomes <- ncol(sp_panaroo) - 14
  
  if (sp_num_genomes < 10) {
    custom_core_needed <- c(custom_core_needed, sp)
    next
  }
  
  core_genes[[sp]] <- rownames(sp_panaroo)[which(sp_panaroo$No..isolates == sp_num_genomes)]

  # Write genes out as text files.
  write.table(x = core_genes[[sp]],
              file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/gene_sets/', sp, '.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)

}

saveRDS(object = core_genes,
        file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/panaroo_only_potential_core.rds')


# Then, also parse out core genes based on core marker genes used by CheckM for each species.
panaroo_passed_core <- list()

for (sp in custom_core_needed) {
  
  panaroo_passed_core[[sp]] <- c()
  
  sp_checkm_marker_hits_file <- paste("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/CheckM_marker_genes/", sp, "_marker_hits.tsv.gz", sep = "")
  sp_checkm_marker_hit_info <- read.table(sp_checkm_marker_hits_file, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  
  # Remove marker genes that have multiple hits in the same genome.
  sp_checkm_marker_hits <- sp_checkm_marker_hit_info$Gene.Id[grep(",", sp_checkm_marker_hit_info$Gene.Id, invert = TRUE)]
  
  sp_panaroo_outfile <- paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence_roary.csv.gz", sep = "")
  sp_panaroo <- read.table(sp_panaroo_outfile, header = TRUE, sep = ",", comment.char = "", quote = "")
  rownames(sp_panaroo) <- paste(sp, sp_panaroo$Gene, sep = "_")
  
  sp_num_genomes <- ncol(sp_panaroo) - 14
  
  colnames(sp_panaroo[, 15:ncol(sp_panaroo)])
  
  sp_potential_core <- list()

  for (gene in rownames(sp_panaroo)) {
    cds_ids <- sp_panaroo[gene, 15:ncol(sp_panaroo)]

    match_count = 0
    for (cds in cds_ids) {
      if (cds %in% sp_checkm_marker_hits) {
        match_count <- match_count + 1
      }
    }

    if (match_count == sp_num_genomes) {
      panaroo_passed_core[[sp]] <- c(panaroo_passed_core[[sp]], gene)
    }

  }

  # Write genes out as text files.
  write.table(x = panaroo_passed_core[[sp]],
              file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/gene_sets/', sp, '.txt', sep = ''),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}


saveRDS(object = panaroo_passed_core,
        file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/checkm_panaroo_passed_potential_core.rds')
