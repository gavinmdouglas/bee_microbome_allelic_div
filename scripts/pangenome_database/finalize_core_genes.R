rm(list = ls(all.names = TRUE))

# Get core genes after filtering out genes that were < 200 bp in length and/or that were in multi-sequence clusters.

panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)

trimmed_genes <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.bed",
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

singleton_members <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/singleton_members.txt.gz",
                                stringsAsFactors = FALSE, header = FALSE, sep = "\t")$V1

passing_genes <- intersect(trimmed_genes, singleton_members)

core_genes_filt <- core_genes

for (sp in names(core_genes)) {

  core_genes_filt[[sp]] <- intersect(core_genes[[sp]], passing_genes)
  
  print(sp)

  print(length(core_genes_filt[[sp]]) / length(core_genes[[sp]]))
  
  outfile <- paste("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/filt_gene_sets/",
                   sp,
                   ".txt", sep = "")
  
  write.table(x = core_genes_filt[[sp]], file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)

}

saveRDS(object = core_genes_filt, file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/core_genes.singletons.above_len.rds")
