rm(list = ls(all.names = TRUE))

# Identify set of single-copy accessory genes per species.
# Also, make sure that all core genes identified are ubiquitous and single-copy.
# Last, only consider genes above length cut-off.

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

core_genes <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/core_INTERMEDIATE/RDS_working/core_genes.singletons.above_len.rds')

# Read in genes that passed trimming step.
trimmed_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.bed',
                            header = FALSE, sep = '\t', stringsAsFactors = FALSE)

for (sp in all_species) {

  sp_panaroo_file <- paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_panaroo/', sp, '/gene_presence_absence.csv.gz', sep = '')
  sp_panaroo <- read.table(file = sp_panaroo_file, header = TRUE, sep = ',', stringsAsFactors = FALSE, quote = '', comment.char = '')

  rownames(sp_panaroo) <- paste(sp, sp_panaroo$Gene, sep = '_')
  
  sp_panaroo <- sp_panaroo[, 4:ncol(sp_panaroo), drop = FALSE]
  
  # Make sure that all core genes are in this table.
  if (length(setdiff(core_genes[[sp]], rownames(sp_panaroo))) > 0) {
    stop('Some core genes missing for ', sp)
  }
  
  # Exclude all genes below length cut-off (or otherwise filtered out for read mapping).
  sp_panaroo <- sp_panaroo[intersect(rownames(sp_panaroo), trimmed_genes$V1), , drop = FALSE]
  
  # Ignore any multi-copy genes.
  multicopy_row_i <- numeric()
  for (genomeid in colnames(sp_panaroo)) {
    multicopy_row_i <- unique(c(multicopy_row_i, grep(';', sp_panaroo[, genomeid])))
  }
  if (length(multicopy_row_i) > 0) {
    sp_panaroo <- sp_panaroo[-multicopy_row_i, , drop = FALSE]
  }

  # Make sure that all core genes are ubiquitous in Panaroo table.
  ubiq_genes <- rownames(sp_panaroo)[which(rowSums(sp_panaroo != '') == ncol(sp_panaroo))]

  final_core_genes <- intersect(core_genes[[sp]], ubiq_genes)

  # Now identify accessory genes, so exclude core genes (even non-final core genes).
  sp_panaroo <- sp_panaroo[setdiff(rownames(sp_panaroo), core_genes[[sp]]), , drop = FALSE]
  
  # Make sure no accessory genes are ubiquitous after this step.
  # (as long as there are at least two genomes).
  if (ncol(sp_panaroo) > 1) {
    ubiq_gene_i <- which(rowSums(sp_panaroo != '') == ncol(sp_panaroo))
    if (length(ubiq_gene_i) > 0) {
      sp_panaroo <- sp_panaroo[-ubiq_gene_i, , drop = FALSE]
    }
  }

  write.table(x = final_core_genes,
              file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/core/', sp, '.txt', sep = ''),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

  write.table(x = rownames(sp_panaroo),
              file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory/', sp, '.txt', sep = ''),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}
