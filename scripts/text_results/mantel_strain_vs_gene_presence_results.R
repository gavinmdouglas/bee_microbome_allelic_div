rm(list = ls(all.names = TRUE))

mantel_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_presence_mantel.tsv.gz',
                             header = TRUE, sep = '\t', row.names = 1)


# Only consider species with at least 10 strains across at least 10 samples,
# to help reduce noise.
mantel_summary <- mantel_summary[which(mantel_summary$num_strains >= 10 & mantel_summary$num_samples >= 10), ]

# Mean, median, and SD of Mantel statistics across species.
format(round(mean(mantel_summary$mantel_kendall), 2), nsmall = 2)
format(round(median(mantel_summary$mantel_kendall), 2), nsmall = 2)
format(round(sd(mantel_summary$mantel_kendall), 2), nsmall = 2)

# Number of significant and non-significant (BH < 0.05)
mantel_summary$mantel_BH <- p.adjust(mantel_summary$mantel_p, 'BH')
length(which(mantel_summary$mantel_BH < 0.05))
length(which(mantel_summary$mantel_BH >= 0.05))

# Details on non-significant species.
format(round(mantel_summary[which(mantel_summary$mantel_BH >= 0.05), ], 2), nsmall = 2)
