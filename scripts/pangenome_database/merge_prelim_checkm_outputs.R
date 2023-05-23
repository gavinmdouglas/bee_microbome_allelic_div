# CheckM (with same parameters and version) was run on different genome subsets at different times.
# The first batch of genomes actually included many non-Apis mellifera-associated genomes, so should
# be ignored. This script was used to merge the tables, and to subset them to the genomes of interest only.

rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/honey_bee/ref_genomes/checkm_output/")

orig_checkm_output <- read.table("Carrie_orig_genome_checkm_output.tsv",
                                 header = TRUE, sep ="\t", stringsAsFactors = FALSE)

new_checkm_output <- read.table("2021_12_19_checkm_additional_downloads.tsv",
                                header = TRUE, sep ="\t", stringsAsFactors = FALSE)
new_checkm_output <- new_checkm_output[, -which(colnames(new_checkm_output) == "Marker.lineage.id")]

new_checkm_output_Lactobacillus_sp <- read.table("2022_01_05_checkm_additional_Lactobacillus_sp.txt",
                                header = TRUE, sep ="\t", stringsAsFactors = FALSE)
new_checkm_output_Lactobacillus_sp <- new_checkm_output_Lactobacillus_sp[, -which(colnames(new_checkm_output_Lactobacillus_sp) == "Marker.Lineage.ID")]


colnames(orig_checkm_output) <- colnames(new_checkm_output)
colnames(new_checkm_output_Lactobacillus_sp) <- colnames(new_checkm_output)

checkm_output <- rbind(orig_checkm_output, new_checkm_output)
checkm_output <- rbind(checkm_output, new_checkm_output_Lactobacillus_sp)

checkm_output$accession <- gsub("^GCA_", "GCAXXXXX", checkm_output$Bin)
checkm_output$accession <- gsub("^GCF_", "GCFXXXXX", checkm_output$accession)
checkm_output$accession <- gsub("_.*$", "", checkm_output$accession)
checkm_output$accession <- gsub("XXXXX", "_", checkm_output$accession)


# Remove duplicated id ("GCF_900094785.1")
checkm_output <- checkm_output[-which(duplicated(checkm_output$accession)), ]

# Check that no accessions duplicated between RefSeq / Genbank (i.e., same id except GCF/GCA swapped)
accession_GCA_only <- gsub("GCF_", "GCX_", checkm_output$accession)
which(duplicated(accession_GCA_only))

# Read in accessions (note that GCF and GCA were both changed)
accessions_of_interest <- as.character()
species_names <- gsub(".tsv$", "", list.files("../adding_new_microbiota_genomes/accessions_to_process/"))
accession_filepaths <- list.files("../adding_new_microbiota_genomes/accessions_to_process/", full.names = TRUE)
for (i in 1:length(accession_filepaths)) {
  accessions_of_interest <- c(accessions_of_interest,
                              read.table(accession_filepaths[i],
                                         stringsAsFactors = FALSE,
                                         sep = "\t", header = TRUE)$accession)
}

# Many genome accessions were the same, but differed due to being in GenBank or RefSeq (GCA vs GCF).
# Changed accession labels to all be "GCA.GCF" to avoid this problem.
rownames(checkm_output) <- gsub('GC._', 'GCA.GCF_', checkm_output$accession)
checkm_output_subset <- checkm_output[gsub('GC._', 'GCA.GCF_', accessions_of_interest), ]

write.table(x = checkm_output_subset,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/CheckM_output.tsv',
            col.names = NA, row.names = TRUE, quote = FALSE, sep = '\t')
