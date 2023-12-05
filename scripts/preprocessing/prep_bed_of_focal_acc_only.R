rm(list = ls(all.names = TRUE))

# Subset bedfiles to only be focal accessory genes of interest.
filt_genes <- character()

gene_id_files <- list.files(path = '/scratch/gdouglas/projects/honey_bee/gene_sets/accessory',
                            pattern = '.txt.gz', full.names = TRUE)

species_to_ignore <- c()

for (gene_id_file in gene_id_files) {
  filt_genes <- c(filt_genes,
                  read.table(file = gene_id_file, header = FALSE, stringsAsFactors = FALSE)$V1)
}

old_bed <- read.table(file = "/scratch/gdouglas/projects/honey_bee/ref_annot/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed",
                      header = FALSE, sep = '\t', stringsAsFactors = FALSE)
length(setdiff(filt_genes, old_bed$V1))

head(filt_genes[which(filt_genes %in% old_bed$V1)])

old_bed <- old_bed[which()]