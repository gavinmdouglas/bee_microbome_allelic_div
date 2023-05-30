rm(list = ls(all.names = TRUE))

# Get final number of genes, as well as the num/percentage removed at the trimming and clustering step.
initial_genes <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.bed",
                             header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

trimmed_genes <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.bed",
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1

singleton_members <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/singleton_members.txt.gz",
                                stringsAsFactors = FALSE, header = FALSE, sep = "\t")$V1



final_genes <- intersect(trimmed_genes, initial_genes)
final_genes <- intersect(singleton_members, final_genes)

round(((length(initial_genes) - length(trimmed_genes)) / length(initial_genes)) * 100, 2)
round(((length(initial_genes) - length(singleton_members)) / length(initial_genes)) * 100, 2)

sprintf("%.2f", round((length(final_genes) / length(initial_genes)) * 100, 2))

