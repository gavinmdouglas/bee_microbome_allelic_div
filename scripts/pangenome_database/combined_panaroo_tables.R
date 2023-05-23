# Make a simple combined panaroo table for easier access.

rm(list = ls(all.names = TRUE))

species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1

panaroo_out <- list()

for (sp in species) {
  
  panaroo_out[[sp]] <- read.table(paste("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/", sp, "/gene_presence_absence_roary.csv.gz", sep = ""),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  rownames(panaroo_out[[sp]]) <- paste(sp, panaroo_out[[sp]]$Gene, sep = "_")
  
  panaroo_out[[sp]] <- panaroo_out[[sp]][, 1:5]
  
  panaroo_out[[sp]]$species <- sp
  
}

all_panaroo_out <- do.call(rbind, panaroo_out)

rownames(all_panaroo_out) <- gsub("^.*\\.", "", rownames(all_panaroo_out))

saveRDS(object = all_panaroo_out, file = "/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")
