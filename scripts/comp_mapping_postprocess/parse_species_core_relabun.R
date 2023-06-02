rm(list = ls(all.names = TRUE))

# Parse out read depths of core genes for species across all samples (for instances where species were identified as present).

species_presence <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz',
                               header = TRUE, sep = "\t", row.names = 1)

bedgraph_files <- list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/mean_depth_per_site', full.names = TRUE)
all_samples <- basename(gsub('.mean.bedGraph.gz$', '', bedgraph_files))

tmp_in1 <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/mean_depth_per_site/SRR10810002.mean.bedGraph.gz',
                      header = FALSE, sep = '\t', stringsAsFactors = FALSE)
gene_mean_depths <- data.frame(matrix(NA, nrow = nrow(tmp_in1), ncol = length(all_samples)))
colnames(gene_mean_depths) <- all_samples
rownames(gene_mean_depths) <- tmp_in1$V1

for (bedgraph in bedgraph_files) {
  samp <- basename(gsub('.mean.bedGraph.gz$', '', bedgraph))
  tmp_in <- read.table(file = bedgraph, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  gene_mean_depths[tmp_in$V1, samp] <- tmp_in$V4
}


# Save file for future reference.
write.table(x = gene_mean_depths, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_mean_depth.tsv',
            col.names = NA, row.names = TRUE, quote = FALSE, sep = '\t')


# Restrict to core genes.
core_genes <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/core_genes.singletons.above_len.rds")
all_core_genes <- as.character(unlist(core_genes))

core_gene_mean_depths <- gene_mean_depths[all_core_genes, ]

species_mean_depth <- species_presence
species_mean_depth[! is.na(species_mean_depth)] <- 0

for (sp in colnames(species_mean_depth)) {

  for (samp in rownames(species_mean_depth)) {

    # Skip cases where species is missing in sample based on breadth of coverage.
    if (species_presence[samp, sp] == 0) { next }

    sp_core_gene_depth <- core_gene_mean_depths[grep(sp, rownames(core_gene_mean_depths)), samp]

    species_mean_depth[samp, sp] <- mean(sp_core_gene_depth[which(sp_core_gene_depth > 0)])

  }

}

# Remove Bombella_sp as it was never called as present.
species_mean_depth <- species_mean_depth[, -which(colnames(species_mean_depth) %in% c("Bombella_sp"))]

species_mean_depth_rel <- species_mean_depth
species_mean_depth_rel <- data.frame(sweep(x = species_mean_depth, MARGIN = 1, STATS = rowSums(species_mean_depth), FUN = '/')) * 100

write.table(x = species_mean_depth, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_mean_depth.tsv',
            col.names = NA, row.names = TRUE, quote = FALSE, sep = '\t')

write.table(x = species_mean_depth_rel, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_rel_abun.tsv',
            col.names = NA, row.names = TRUE, quote = FALSE, sep = '\t')
