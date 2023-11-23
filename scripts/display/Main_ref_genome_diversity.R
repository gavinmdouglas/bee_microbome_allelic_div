rm(list = ls(all.names = TRUE))

# Breakdown of diversity in reference genomes, based on several metrics.

library(cowplot)
library(ggbeeswarm)
library(ggplot2)
library(reshape2)

source('~/scripts/bee_microbome_allelic_div/scripts/functions.R')

# First plot gene vs species tree distances.
tree_results <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_vs_species_tree_dist.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Remove genes for which neither tree distances nor D values could be computed.
# Only consider genes with at least 5 tips, as otherwise the tree distance is unstable, and is almost always 0 or 1.
tree_results <- tree_results[-which(is.na(tree_results$tree_dist_subset5)), ]

tree_results$type <- 'Core genes'
tree_results[which(tree_results$prevalence < 1), 'type'] <- 'Accessory genes'
tree_results$type <- factor(tree_results$type, levels = c('Core genes', 'Accessory genes'))

tree_results$species <- cleanup_species_names(tree_results$species)

# Sort species by mean values.
mean_tree_dist5_by_species <- aggregate(tree_dist_subset5 ~ species, data = tree_results, FUN = mean)
mean_tree_dist5_by_species <- mean_tree_dist5_by_species[order(mean_tree_dist5_by_species$tree_dist_subset5), ]

tree_results$species <- factor(tree_results$species, levels = mean_tree_dist5_by_species$species)

names(species_plot_colours) <- cleanup_species_names(names(species_plot_colours))

tree_dist_panel <- ggplot(data = tree_results, aes(x = tree_dist_subset5, y = species)) +
  geom_violin(aes(fill=species, col=species)) +
  geom_boxplot(size=0.5, alpha=0.2, outlier.shape = NA) +
  theme_bw() +
  xlab('Normalized distance between gene and species trees') +
  ylab('Species') +
  facet_wrap(type ~ .) +
  theme(axis.text.y=element_text(face="italic")) +
  scale_fill_manual(values = species_plot_colours[levels(tree_results$species)]) +
  scale_colour_manual(values = species_plot_colours[levels(tree_results$species)]) +
  theme(legend.position = "none")

# Then plot Hamming distances for core genome and accessory genes against each other,
# and against Jaccard of accessory gene presence/absence.
genome_hamming <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genome_hamming.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)
genome_hamming[, c('seq1', 'seq2')] <- data.frame(t(apply(genome_hamming[, c('seq1', 'seq2')], 1, function(row) sort(row))))

rownames(genome_hamming) <- paste(genome_hamming$seq1, genome_hamming$seq2, sep = '_')


individual_hamming <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/individual_gene_hamming.tsv.gz',
                                 header = FALSE, sep = '\t', stringsAsFactors = FALSE)

# Exclude core genes and those with multiple copies per genome.
individual_hamming <- individual_hamming[which(! individual_hamming$V1 %in% tree_results[which(tree_results$type == 'Core genes'), 'gene']), ]

multi_copy_genes <- unique(c(individual_hamming$V1[grep('_\\d$', individual_hamming$V2)],
                             individual_hamming$V1[grep('_\\d$', individual_hamming$V3)]))

individual_hamming <- individual_hamming[which(! individual_hamming$V1 %in% multi_copy_genes), ]

individual_hamming$pair <- paste(individual_hamming$V2, individual_hamming$V3, sep = '_')

mean_acc_identity <- aggregate(x = V7 ~ pair, data = individual_hamming, FUN = mean)
rownames(mean_acc_identity) <- mean_acc_identity$pair

sd_acc_identity <- aggregate(x = V7 ~ pair, data = individual_hamming, FUN = sd)
rownames(sd_acc_identity) <- sd_acc_identity$pair

# Read in gene Jaccard vs. phylogenetic distance (although note that I decided to use the core genome Hamming distance-based percenty identity,
# rather than branch length distance, just to make it easier to interpret distance across species.
acc_jaccard <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/pairwise_genome_acc_gene_jaccard.tsv.gz',
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE)
acc_jaccard[, c('genome1', 'genome2')] <- data.frame(t(apply(acc_jaccard[, c('genome1', 'genome2')], 1, function(row) sort(row))))
rownames(acc_jaccard) <- paste(acc_jaccard$genome1, acc_jaccard$genome2, sep = '_')

all_dist <- data.frame(Species = genome_hamming$annot,
                       genome_identity = genome_hamming$percent_identity,
                       acc_jaccard = acc_jaccard[rownames(genome_hamming), 'acc_gene_jaccard'],
                       acc_mean_identity = NA,
                       acc_sd_identity = NA)
rownames(all_dist) <- rownames(genome_hamming)

all_dist[rownames(mean_acc_identity), 'acc_mean_identity'] <- mean_acc_identity$V7
all_dist[rownames(sd_acc_identity), 'acc_sd_identity'] <- sd_acc_identity$V7
all_dist$acc_cv_identity <- (all_dist$acc_sd_identity / all_dist$acc_mean_identity) * 100

all_dist$Species <- cleanup_species_names(all_dist$Species)

genome_identity_vs_gene_content_jaccard <- ggplot(
  data = all_dist, aes(x = genome_identity, y = 1 - acc_jaccard, colour = Species)) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=species_plot_colours[unique(all_dist$Species)]) +
  theme(legend.text = element_text(face='italic')) +
  ylab('Accessory gene content Jaccard similarity') +
  xlab('Core genome percent identity') +
  guides(color = guide_legend(ncol = 2))

genome_identity_vs_mean_acc_identity <- ggplot(
  data = all_dist, aes(x = genome_identity, y = acc_mean_identity, colour = Species)) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=species_plot_colours[unique(all_dist$Species)]) +
  theme(legend.text = element_text(face='italic')) +
  ylab('Mean accessory gene percent identity') +
  xlab('Core genome percent identity') +
  guides(color = guide_legend(ncol = 2))

legend <- get_legend(genome_identity_vs_mean_acc_identity)

top_row <- plot_grid(tree_dist_panel, legend,
                     nrow = 1,
                     labels = c('a', ''))

bottom_row <- plot_grid(genome_identity_vs_gene_content_jaccard + theme(legend.position = "none"),
                        genome_identity_vs_mean_acc_identity + theme(legend.position = "none"),
                        labels = c('b', 'c'), nrow = 1)

combined_plot <- plot_grid(top_row, bottom_row, nrow = 2)

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_ref_genome_diversity.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 600,
       width = 12,
       height = 8)
