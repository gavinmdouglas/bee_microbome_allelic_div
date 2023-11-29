rm(list = ls(all.names = TRUE))

# Subset BEDfile specifically to accessory genes to analyze.

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

all_acc_genes <- character()

ref_bed <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                      header = FALSE, sep = '\t', stringsAsFactors = FALSE)

for (sp in all_species) {
  sp_acc_genes <- read.table(file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory/', sp, '.txt', sep = ''),
                             header = FALSE, stringsAsFactors = FALSE)$V1
  all_acc_genes <- c(all_acc_genes, sp_acc_genes)
}

length(which(ref_bed$V1 %in% all_acc_genes))

missing_genes <- setdiff(ref_bed$V1, all_acc_genes)