rm(list = ls(all.names = TRUE))

# Call species as present across specific samples based on core gene presence profiles.
# Based on both 25% and 90% of core genes being present.

library(ggplot2)

all_species <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz",
                          stringsAsFactors = FALSE)$V1

all_present <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz",
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)

# Breakdown of % core genes (of those above size cut-off) called as present per sample per species.
all_samples <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/all_five_datasets.txt.gz",
                          stringsAsFactors = FALSE)$V1

core_gene_breakdown <- data.frame(matrix(NA, nrow = length(all_species) * length(all_samples), ncol = 4))
colnames(core_gene_breakdown) <- c("Species", "Sample", "Num_core", "Percent_core")

core_gene_breakdown$Species <- rep(x = all_species, each = length(all_samples))
core_gene_breakdown$Sample <- rep(x = all_samples, times = length(all_species))

# Read in core gene sets.
core_genes <- list()
for (sp in all_species) {
  print(sp)
  core_gene_file <- paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/core/', sp, '.txt.gz', sep = '')
  sp_core_genes <- read.table(file = core_gene_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  missing_genes <- setdiff(sp_core_genes, rownames(all_present))
  
  if (length(missing_genes) > 0) { 
   print('Dropping this many genes: ', length(missing_genes))
   print('Out of ', length(sp_core_genes))
  }

  core_genes[[sp]] <- intersect(sp_core_genes, rownames(all_present))
}


for (sp in all_species) {
  
  sp_core <- core_genes[[sp]]
  
  if (length(which(! sp_core %in% rownames(all_present))) > 0) { stop('Error - expected genes missing!') }
  
  for (samp in all_samples) {
    
    row_i <- which(core_gene_breakdown$Species == sp & core_gene_breakdown$Sample == samp)

    core_gene_breakdown[row_i, "Num_core"] <- length(which(all_present[sp_core, samp] > 0))

    core_gene_breakdown[row_i, "Percent_core"] <- (core_gene_breakdown[row_i, "Num_core"] / length(sp_core)) * 100

  }

}


# Get species presence profile based on 25% and 90% core genes present.
species_25per <- data.frame(matrix(0, nrow = length(all_samples), ncol = length(all_species)))
colnames(species_25per) <- all_species
rownames(species_25per) <- all_samples

species_90per <- species_25per

for (sp in all_species) {

  sp_core <- core_genes[[sp]]

  for (samp in all_samples) {
 
    percent_core <- core_gene_breakdown[which(core_gene_breakdown$Species == sp & core_gene_breakdown$Sample == samp), "Percent_core"]

    if (percent_core >= 25) {
      species_25per[samp, sp] <- 1
    }

    if (percent_core >= 90) {
      species_90per[samp, sp] <- 1
    }

  }

}

write.table(x = species_25per, file = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_25percent.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")

write.table(x = species_90per, file = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
