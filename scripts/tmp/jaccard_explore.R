rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

jaccard_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE)
jaccard_summary_clean <- jaccard_summary
jaccard_summary$datatype <- tolower(jaccard_summary$datatype)

strainfacts_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'strainfacts'), -2]
rownames(strainfacts_jaccard) <- paste(strainfacts_jaccard$species, strainfacts_jaccard$dataset)
strainfacts_jaccard <- strainfacts_jaccard[, -c(1:2)]
colnames(strainfacts_jaccard) <- paste('strainfacts', colnames(strainfacts_jaccard), sep = '_')

straingst_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'straingst'), -2]
rownames(straingst_jaccard) <- paste(straingst_jaccard$species, straingst_jaccard$dataset)
straingst_jaccard <- straingst_jaccard[, -c(1:2)]
colnames(straingst_jaccard) <- paste('straingst', colnames(straingst_jaccard), sep = '_')

gene_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'gene'), -2]
rownames(gene_jaccard) <- paste(gene_jaccard$species, gene_jaccard$dataset)
gene_jaccard <- gene_jaccard[, -c(1:2)]
colnames(gene_jaccard) <- paste('gene', colnames(gene_jaccard), sep = '_')

allele_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'allele'), -2]
rownames(allele_jaccard) <- paste(allele_jaccard$species, allele_jaccard$dataset)
allele_jaccard <- allele_jaccard[, -c(1:2)]
colnames(allele_jaccard) <- paste('allele', colnames(allele_jaccard), sep = '_')


intersecting_rows <- intersect(rownames(strainfacts_jaccard), rownames(straingst_jaccard))
intersecting_rows <- intersect(intersecting_rows, rownames(gene_jaccard))
intersecting_rows <- intersect(intersecting_rows, rownames(allele_jaccard))

combined_data <- do.call(cbind, list(strainfacts_jaccard[intersecting_rows, ],
                                     straingst_jaccard[intersecting_rows, ],
                                     gene_jaccard[intersecting_rows, ],
                                     allele_jaccard[intersecting_rows, ]))

combined_data <- combined_data[, grep("mean", colnames(combined_data))]

pairs(combined_data)