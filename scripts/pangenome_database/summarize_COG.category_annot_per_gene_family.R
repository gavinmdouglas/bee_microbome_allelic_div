rm(list = ls(all.names = TRUE))

library(stringr)

# Get annotations for gene families
# Note that the COG category classification in the eggNOG output is based on an outdated version that doesn't include mobile elements.

# Read in input
eggNOG_output <- read.table(file = '/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/2022_01_07_func_annotating/eggNOG_mapper/eggnog_mapper_annot/derep_gene_annot.emapper.annotations.tsv',
                            header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

COG2category <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

# Get COG info for clusters based on COG definition of any orthogroups in EggNOG hierarchy definition.
eggNOG_output$all_COG <- sapply(eggNOG_output$eggNOG_OGs, function(x) { paste(grep("^COG", unique(gsub("@.*$", "", strsplit(x, ",")[[1]])), value = TRUE), collapse = ",") })
eggNOG_output$Original_COG_category <- eggNOG_output$COG_category
eggNOG_output$COG_category <- "-"

genes_with_COG_i <- which(eggNOG_output$all_COG != "")

eggNOG_output$COG_category[genes_with_COG_i] <- sapply(genes_with_COG_i,
                                                               function(i) {
                                                                 
                                                                 all_COGs <- eggNOG_output[i, "all_COG"]
                                                                 categories <- as.character()
                                                                 COGs <- strsplit(all_COGs, ",")[[1]]
                                                                 for (COG in COGs) {
                                                                   categories <- c(categories, COG2category[COG, "category"])
                                                                 }
                                                                 
                                                                 return(paste(sort(unique(categories)), collapse = ","))
                                                               })

write.table(x = eggNOG_output,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv',
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")
