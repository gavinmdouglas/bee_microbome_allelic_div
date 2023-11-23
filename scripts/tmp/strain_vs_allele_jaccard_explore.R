rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

jaccard_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE)

min_num_samples <- 20
min_num_other_features <- 3
jaccard_summary <- jaccard_summary[which(jaccard_summary$num_samples >= min_num_samples), ]
jaccard_summary <- jaccard_summary[which(jaccard_summary$num_alleles >= min_num_other_features), ]
jaccard_summary <- jaccard_summary[which(jaccard_summary$num_strainfacts_strains >= min_num_other_features), ]
jaccard_summary <- jaccard_summary[which(jaccard_summary$num_straingst_strains >= min_num_other_features), ]

sufficiently_frequent_species <- names(table(jaccard_summary$species))[which(table(jaccard_summary$species) >= 50)]

jaccard_summary <- jaccard_summary[which(jaccard_summary$species %in% sufficiently_frequent_species), ]

# Convert to long-format.
jaccard_long <- data.frame(species=c(jaccard_summary$species, jaccard_summary$species),
                          mean_min=c(jaccard_summary$strainfacts_mean_min, jaccard_summary$straingst_mean_min),
                          OR=c(jaccard_summary$strainfacts_mean_min / jaccard_summary$mean_permuted_strainfacts_vs_allele_jaccard,
                               jaccard_summary$straingst_mean_min / jaccard_summary$mean_permuted_straingst_vs_allele_jaccard),
                          tool=c(rep('StrainFacts', nrow(jaccard_summary)),
                                 rep('StrainGST', nrow(jaccard_summary))))



mean_min_plot <- ggplot(data = jaccard_long, aes(x=mean_min, y=species)) +
  geom_violin(col='grey', fill='grey') +
  geom_boxplot(alpha=0.8, outlier.shape=NA) +
  theme_bw() +
  #xlab("Strain-level inference tool") +
  ylab("Mean minimum Jaccard distance\nbetween strain and allele profiles") +
  facet_wrap(. ~ tool)

jaccard_long$log2OR <- log2(jaccard_long$OR)

OR_plot <- ggplot(data = jaccard_long, aes(y=OR, x=tool)) +
  geom_violin(col='grey', fill='grey') +
                    geom_boxplot(alpha=0.8, outlier.shape=NA) +
                    theme_bw() +
                    xlab("Strain-level inference tool") +
                    ylab("Odd's ratio\n(Observed/permuted min. Jaccard distance)")
