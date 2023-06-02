rm(list = ls(all.names = TRUE))

# Write out number of final gene families per species, for future reference.
trimmed_genes <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.bed",
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

singleton_members <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/singleton_members.txt.gz",
                                stringsAsFactors = FALSE, header = FALSE, sep = "\t")$V1

passing_genes <- intersect(trimmed_genes, singleton_members)

core_genes <- readRDS(file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/core_genes.singletons.above_len.rds")

final_counts <- data.frame(all = rep(NA, length(core_genes)),
                           core = rep(NA, length(core_genes)))
rownames(final_counts) <- names(core_genes)

for (sp in names(core_genes)) {
  final_counts[sp, ] <- c(length(grep(sp, passing_genes)),
                          length(core_genes[[sp]]))
}

final_counts$species <- rownames(final_counts)
final_counts <- final_counts[, c('species', 'all', 'core')]
write.table(x = final_counts, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_gene_families_per_species.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')