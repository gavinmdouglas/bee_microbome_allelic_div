rm(list = ls(all.names = TRUE))

library(ggplot2)
library(cowplot)
library(ggbeeswarm)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

# Summarize single nucleotide variants from mapped reads per gene.
datasets <- c('Ellegaard_2019', 'Ellegaard_2020', 'Sun_2022', 'Wu_2021', 'Zhang_2022')

breakdown <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/mean_pi_and_nseg_per_gene.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Read in core and accessory gene IDs, add this information to the table, along with species name.
all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

breakdown$partition <- NA
breakdown$species <- NA

for (sp in all_species) {
  sp_core <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/core/', sp, '.txt', sep = ''),
                        header = FALSE, stringsAsFactors = FALSE)$V1

  sp_accessory <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_sets/accessory/', sp, '.txt', sep = ''),
                             header = FALSE, stringsAsFactors = FALSE)$V1
  
  breakdown[which(breakdown$gene %in% sp_core), 'partition'] <- 'Core'
  breakdown[which(breakdown$gene %in% sp_core), 'species'] <- sp
  
  breakdown[which(breakdown$gene %in% sp_accessory), 'partition'] <- 'Accessory'
  breakdown[which(breakdown$gene %in% sp_accessory), 'species'] <- sp
}

breakdown <- breakdown[which(! is.na(breakdown$species)), ]

# Exclude rare species.
species_tallies <- table(breakdown$species)
rare_species <- names(species_tallies)[which(species_tallies < 1000)]
breakdown <- breakdown[which(! breakdown$species %in% rare_species), ]

breakdown$species <- cleanup_species_names(breakdown$species)

dataset_col <- c("#638ccc",
                 "#c57c3c",
                 "#ab62c0",
                 "#72a555",
                 "#ca5670")

colnames(breakdown)[which(colnames(breakdown) == 'dataset')] <- 'Dataset'
breakdown$Dataset <- sub('_', ' ', breakdown$Dataset)

colnames(breakdown)[which(colnames(breakdown) == 'species')] <- 'Species'
breakdown$Species <- factor(breakdown$Species, levels = rev(unique(breakdown$Species)))

breakdown$percent_seg <- (breakdown$mean_nseg / breakdown$mean_num_sites) * 100

breakdown$partition <- factor(breakdown$partition, levels = c('Core', 'Accessory'))

within_distributions <- ggplot(data = breakdown[which(breakdown$type == 'Within'), ],
                                   aes(y = Species,
                                       x = percent_seg,
                                       fill = Dataset)) +
                              geom_boxplot(outlier.size = 0.25) +
                              facet_wrap(. ~ partition) +
                              scale_fill_manual(values = dataset_col) +
                              theme_bw() +
                              xlab('Percent variable sites') +
                              theme(axis.text.y = element_text(face = 'italic'),
                                    plot.title = element_text(hjust = 0.5)) +
                              ggtitle('Within samples')

between_distributions <- ggplot(data = breakdown[which(breakdown$type == 'Between'), ],
                               aes(y = Species,
                                   x = percent_seg,
                                   fill = Dataset)) +
  geom_boxplot(outlier.size = 0.25) +
  facet_wrap(. ~ partition) +
  scale_fill_manual(values = dataset_col) +
  theme_bw() +
  xlab('Percent variable sites') +
  theme(axis.text.y = element_text(face = 'italic'),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('Between samples')

legend <- get_legend(between_distributions)

combined_plot <-  plot_grid(within_distributions,
                           between_distributions,
                           labels = c('a', 'b'),
                           ncol = 1)
ggsave(plot = combined_plot,
       filename = "/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_percent_variable_sites_SNVs.pdf",
       device = "pdf", width = 7, height = 10, units = "in", dpi = 400)
