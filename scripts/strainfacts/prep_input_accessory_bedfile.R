# Remove core genes and multi-seq cluster members from bedfile.
# Also remove genes associated with Apilactobacillus apinorum,  Bombella apis, and Bombella sp, as
# these species were called as present in only a max of <= 2 samples anyway, so the results won't be
# trustworthy.

core_genes <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/core_genes.singletons.above_len.rds")
all_core_genes <- as.character(unlist(core_genes))

singleton_cluster_members <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/singleton_members.txt.gz",
                                        stringsAsFactors = FALSE, header = FALSE)$V1

bedfile <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.bed",
                      stringsAsFactors = FALSE, header = FALSE, sep = "\t")

bedfile <- bedfile[which(! bedfile$V1 %in% all_core_genes), ]

bedfile <- bedfile[which(bedfile$V1 %in% singleton_cluster_members), ]

bedfile <- bedfile[-grep("Apilactobacillus_apinorum", bedfile$V1), ]
bedfile <- bedfile[-grep("Bombella_apis", bedfile$V1), ]
bedfile <- bedfile[-grep("Bombella_sp", bedfile$V1), ]

write.table(x = bedfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t",
            file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed")
