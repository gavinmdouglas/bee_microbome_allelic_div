rm(list = ls(all.names = TRUE))

# Compare StrainFacts vs StrainGE strain inferences.

library(ape)
library(cowplot)
library(ggtree)
library(ggplot2)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

# Read in strain results.
strain_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
names(strain_abun)[which(names(strain_abun) == 'Bifidobacterium_coryneforme_indicum')] <- 'Bifidobacterium_coryneforme'

straingst_strain_files <- list.files(path = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/straingst_strain_outfiles',
                                     pattern = '.strains.tsv.gz', full.names = TRUE)

straingst <- list()
for (straingst_strain_file in straingst_strain_files) {
  file_info <- strsplit(sub('.strains.tsv.gz', '', basename(straingst_strain_file)), '_')[[1]]
  SRR <- file_info[1]
  species <- paste(file_info[2], file_info[3], sep = '_')
  straingst_sample_in <- read.table(file = straingst_strain_file, header = TRUE, sep = '\t')
  if (nrow(straingst_sample_in) == 0 | length(which(straingst_sample_in$rapct > 0)) == 0) { next }
  straingst_sample_in <- straingst_sample_in[which(straingst_sample_in$rapct > 0), ]
  if (nrow(straingst_sample_in) == 0) { next }
  if (! species %in% names(straingst)) { straingst[[species]] <- list() }
  straingst[[species]][[SRR]] <- straingst_sample_in$strain
}


# Number of samples with strain calls
straingst_sample_nums <- sapply(straingst, length)

strainfacts_sample_nums <- sapply(strain_abun,
                          function(x) {
                            sum(sapply(x, ncol))
                          })

sample_num_breakdown <- data.frame(Species = names(straingst_sample_nums),
                                   StrainFacts = NA,
                                   StrainGST = straingst_sample_nums)

rownames(sample_num_breakdown) <- sample_num_breakdown$Species

sample_num_breakdown[names(strainfacts_sample_nums), 'StrainFacts'] <- strainfacts_sample_nums

sample_num_breakdown$Species <- cleanup_species_names(names_vec = sample_num_breakdown$Species, shorten = TRUE)

sample_num_breakdown <- sample_num_breakdown[which(! is.na(sample_num_breakdown$StrainFacts)), ]

straingst_vs_strainfacts_sample_breakdown <- ggplot(data = sample_num_breakdown,
                                              aes(x = StrainGST,
                                                  y = StrainFacts)) +
                                              geom_abline(slope = 1, lwd = 1, lty = 2, col = 'red') +
                                              geom_point(size = 2) +
                                              xlim(c(0, 275)) +
                                              ylim(c(0, 275)) +
                                              ggrepel::geom_text_repel(aes(label = Species),
                                                                       col = 'grey20',
                                                                       hjust = -0.1,
                                                                       fontface = 'italic',
                                                                       size = 3) +
                                              theme_bw() +
                                              xlab('StrainGST - Number of samples with strains') +
                                              ylab('StrainFacts - Number of samples with strains') +
                                              theme(plot.title = element_text(hjust = 0.5)) +
                                              ggtitle('Number of samples with strains, across all datasets')


# Number of strains per species (for the Ellegaard 2019 dataset as an example).
Ellegaard2019_samples <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Ellegaard2019_SRRs.txt.gz',
                                    stringsAsFactors = FALSE, header = FALSE)$V1

strainfacts_strain_nums <- sapply(strain_abun,
                          function(x) {
                            if ('Ellegaard2019' %in% names(x)) {
                              return(nrow(x$Ellegaard2019))
                            } else {
                              return(0)
                            }
                          })

straingst_strain_nums <- sapply(straingst,
                                function(x) {
                                  x <- x[Ellegaard2019_samples]
                                  length(unique(unlist(x)))
                                })

strain_num_breakdown <- data.frame(Species = names(straingst_strain_nums),
                                   StrainFacts = NA,
                                   StrainGST = straingst_strain_nums)
rownames(strain_num_breakdown) <- strain_num_breakdown$Species
strain_num_breakdown[names(strainfacts_strain_nums), 'StrainFacts'] <- strainfacts_strain_nums


strain_num_breakdown$Species <- cleanup_species_names(names_vec = strain_num_breakdown$Species, shorten = TRUE)

strain_num_breakdown <- strain_num_breakdown[which(! is.na(strain_num_breakdown$StrainFacts)), ]


straingst_vs_strainfacts_num_strains <- ggplot(data = strain_num_breakdown,
                                         aes(x = StrainGST,
                                             y = StrainFacts)) +
  geom_abline(slope = 1, lwd = 1, lty = 2, col = 'red') +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = Species),
                           col = 'grey20',
                           hjust = -0.1,
                           fontface = 'italic',
                           size = 3) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Number of strains per species in Ellegaard 2019 dataset') +
  xlab('StrainGST - Number of strains') +
  ylab('StrainFacts - Number of strains') +
  xlim(c(0, 15)) +
  ylim(c(0, 15))


# Then summarize tree-based analyses.
# Many genomes were dropped because they are redundant. Identify these before processing trees.
all_species <- read.table('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt', stringsAsFactors = FALSE)$V1

genomes_to_exclude <- list()
for (sp in all_species) {
  sp_genomes <- gsub('\\.fna$', '', list.files(path = paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes/', sp, sep = ''), pattern = '.fna'))
  
  sp_genomes_retained_files <- read.table(paste('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/strainge/kmer_comparison/references_to_keep/', sp, '.txt', sep = ''),
                                          stringsAsFactors = FALSE)$V1
  sp_genomes_retained <- gsub('\\.hdf5$', '', basename(sp_genomes_retained_files))
  
  if (length(setdiff(sp_genomes_retained, sp_genomes)) > 0) {
    stop('Error - Retained genomes not found in orig set.')
  }
  
  genomes_to_exclude[[sp]] <- setdiff(sp_genomes, sp_genomes_retained)
}


# Read in trees.
trees <- list()
tree_paths <- list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strain_trees/',
                         pattern = '.tree',
                         full.names = TRUE)

for (tree_path in tree_paths) {
  tree_id <- gsub('.tree', '', basename(tree_path))
  sp_id <- gsub('\\..*$', '', tree_id)
  trees[[tree_id]] <- ape::read.tree(tree_path)
  
  if (length(genomes_to_exclude[[sp_id]]) > 0) {
    trees[[tree_id]] <- ape::drop.tip(trees[[tree_id]], genomes_to_exclude[[sp_id]])
  }
}

species <- sort(unique(sapply(tree_paths, function(x) { strsplit(basename(x), '\\.')[[1]][1] })))

# Concordance between StrainFacts and StrainGE for a representative example.

dataset <- 'Ellegaard2019'
sp <- 'Snodgrassella_alvi'
  
strainfacts_strains <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_relabun/', dataset, '.', sp, '.tsv.gz', sep = ''),
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
  
# sample(colnames(strainfacts_strains), 1)
sample_id <- 'SRR7287247'
  
tree_id <- paste(sp, dataset, sep = '.')
  
strainfacts_strains <- rownames(strainfacts_strains)[which(strainfacts_strains[, sample_id] > 0)]
strainfacts_strains <- paste(dataset, strainfacts_strains, sep = '.')
  
strainge_strains <- read.table(paste('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/strainge/straingst_out/', sample_id, '_', sp, '.strains.tsv', sep = ''),
                               header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)$strain
  
strainfacts_nodes_to_highlight <- nodeid(tree=trees[[tree_id]],
                                         label = trees[[tree_id]]$tip.label[which(trees[[tree_id]]$tip.label %in% strainfacts_strains)])
  
strainge_nodes_to_highlight <- nodeid(tree=trees[[tree_id]],
                                      label = trees[[tree_id]]$tip.label[which(trees[[tree_id]]$tip.label %in% strainge_strains)])
  
tip_group <- data.frame(taxa = trees[[tree_id]]$tip.label,
                        Presence = 'Absent',
                        Tool = 'StrainGST (ref. genomes)')

tip_group[c(strainge_nodes_to_highlight, strainfacts_nodes_to_highlight), 'Presence'] <- 'Present'
tip_group[grep(dataset, tip_group$taxa), 'Tool'] <- 'StrainFacts (de novo)'

  
example_tree <- ggtree(trees[[tree_id]]) %<+% tip_group +
                     geom_tippoint(aes(colour=Presence, shape=Tool), size=3.5) +
                     scale_colour_manual(values=c(Absent="grey10", Present="#FC4E07")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Example strain tree')


# Observed percentile vs. real data
# Compare how low the branch-length distance is between strainFacts and strainGE strains per sample vs. all other possible comparisons of the same number of strains for each tool.
# Key metric is the "p-value", or the percentile of the observed distance on the background of all (or at least a large number of) possible strain combinations.

strain_dists_summary <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strainfacts_vs_strainge_dist_summary.tsv.gz',
                                   header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Restrict to samples with >= 100 strain combos.
strain_dists_summary_atleast100 <- strain_dists_summary[which(strain_dists_summary$num_background_combos >= 100), ]

wilcox.test(strain_dists_summary_atleast100$background_dist_percentile, mu = 0.5)

percentile_hist <- ggplot(data = strain_dists_summary_atleast100, aes(x = background_dist_percentile)) +
                          geom_histogram(bins = 100) +
                          theme_bw() +
                          ylab('Number of samples') +
                          xlab('Percentile of observed StrainFacts vs. StrainGST\nstrain distance in permuted background') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Permutation test summary')

combined_plot <- plot_grid(straingst_vs_strainfacts_sample_breakdown,
                          straingst_vs_strainfacts_num_strains,
                          example_tree,
                          percentile_hist,
                          labels = c('a', 'b', 'c', 'd'))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_strainfacts_vs_straingst.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 9)
