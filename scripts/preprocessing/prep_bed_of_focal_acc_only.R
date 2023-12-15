rm(list = ls(all.names = TRUE))

# Subset bedfiles to only be focal accessory genes of interest.
filt_genes <- character()

species_presence <- read.table(file = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

present_species <- colnames(species_presence)[which(colSums(species_presence) > 3)]

gene_id_files <- list.files(path = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory',
                            pattern = '.txt.gz', full.names = TRUE)

species_to_ignore <- c()

for (gene_id_file in gene_id_files) {
  file_sp <- gsub(".txt.gz$", "", basename(gene_id_file))
  if (file_sp %in% present_species) {
    filt_genes <- c(filt_genes,
                    read.table(file = gene_id_file, header = FALSE, stringsAsFactors = FALSE)$V1)
  }
}

old_bed <- read.table(file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed",
                      header = FALSE, sep = '\t', stringsAsFactors = FALSE)
length(setdiff(filt_genes, old_bed$V1))

new_bed <- old_bed[which(old_bed$V1 %in% filt_genes), ]

write.table(x = new_bed,
            file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.focal.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
