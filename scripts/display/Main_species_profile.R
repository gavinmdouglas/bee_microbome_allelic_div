rm(list = ls(all.names = TRUE))

# Main figure
  # Heatmap of species presence/absence with sample metadata information shown as well.
  # Species relative abundance.

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(reshape2)

species_90per <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz',
                            row.names = 1, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

colnames(species_90per) <- gsub('_', ' ', colnames(species_90per))
colnames(species_90per) <- gsub('sp$', 'sp.', colnames(species_90per))
colnames(species_90per) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', colnames(species_90per))

species_90per_char <- species_90per
species_90per_char[species_90per_char == 0] <- 'Absent'
species_90per_char[species_90per_char == '1'] <- 'Present'

Ellegaard_2019_SRRs <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Ellegaard2019_SRRs.txt.gz', stringsAsFactors = FALSE)$V1
Ellegaard_2020_SRRs <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Ellegaard2020_SRRs.txt.gz', stringsAsFactors = FALSE)$V1
Sun_2022_SRRs <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Sun2022_SRRs.txt.gz', stringsAsFactors = FALSE)$V1
Wu_2021_SRRs <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Wu2021_SRRs.txt.gz', stringsAsFactors = FALSE)$V1
Zhang_2022_SRRs <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/Zhang2022_SRRs.txt.gz', stringsAsFactors = FALSE)$V1

dataset_breakdown <- data.frame(Dataset = c(rep(x = 'Ellegaard 2019', length(Ellegaard_2019_SRRs)),
                                            rep(x = 'Ellegaard 2020', length(Ellegaard_2020_SRRs)),
                                            rep(x = 'Sun 2022', length(Sun_2022_SRRs)),
                                            rep(x = 'Wu 2021', length(Wu_2021_SRRs)),
                                            rep(x = 'Zhang 2022', length(Zhang_2022_SRRs))))
rownames(dataset_breakdown) <- c(Ellegaard_2019_SRRs, Ellegaard_2020_SRRs, Sun_2022_SRRs, Wu_2021_SRRs, Zhang_2022_SRRs)

dataset_colours <- c("seagreen1",
                     "#5fa375",
                     "black",
                     "lightsteelblue1",
                     "yellow2")

names(dataset_colours) <- unique(dataset_breakdown$Dataset)
dataset_breakdown_annot <- ComplexHeatmap::rowAnnotation(Dataset=dataset_breakdown$Dataset,
                                                         col = list(Dataset = dataset_colours),
                                                         show_annotation_name = FALSE)

# Plot all datasets together, with only metadata being what dataset the sample is from.
# For this purpose, perform hierarchical clustering of samples and species based on all data.
species_cluster <- stats::hclust(dist(t(species_90per), method = 'binary'), method = 'complete')
sample_cluster <- stats::hclust(dist(species_90per, method = 'binary'), method = 'complete')

species_ordered <- species_cluster$labels[species_cluster$order]
samples_ordered <- sample_cluster$labels[sample_cluster$order]

column_label_col <- c('black', # 'Apilactobacillus apinorum'
                      'black', # 'Apilactobacillus kunkeei'
                      'yellow2', # 'Bartonella apis'
                      '#63bfaf', # 'Bifidobacterium asteroides'
                      '#62d07a', # 'Bifidobacterium cor./indicum'
                      'black', # 'Bombella apis'
                      'black', # 'Bombella sp.'
                      '#d56c24', # 'Bombilactobacillus mellifer'
                      '#bf4d81', # 'Bombilactobacillus mellis'
                      '#d2cb3b', # 'Commensalibacter sp.'          
                      'brown', #'Frischella perrara'
                      'royalblue1', # 'Gilliamella apicola'
                      'skyblue1', # 'Gilliamella apis'
                      'mediumblue', # 'Gilliamella sp.'
                      'grey40', # 'Lactobacillus apis'            
                      'grey50', #'Lactobacillus helsingborgensis'
                      'grey60', # 'Lactobacillus kimbladii'
                      'grey70', # 'Lactobacillus kullabergensis'
                      'grey80', # 'Lactobacillus melliventris'
                      'black', # 'Serratia marcescens'           
                      '#9f4ad4'# 'Snodgrassella alvi'
)

species_presence <- Heatmap(as.matrix(species_90per_char),
                            col = c('grey95', 'cornflowerblue'),
                            heatmap_legend_param = list(title = 'Species'),
                            column_labels = ComplexHeatmap::gt_render(paste0("*", colnames(species_90per_char), "*")),
                            cluster_rows = sample_cluster,
                            cluster_columns = species_cluster,
                            show_row_names = FALSE,
                            show_column_names = TRUE,
                            right_annotation = dataset_breakdown_annot,
                            column_names_rot = 45,
                            row_title_rot = 0,
                            row_title = 'Sample',
                            row_title_gp = gpar(fontsize = 12),
                            column_names_gp = gpar(fontsize = 12, col = column_label_col))


# Also plot species' relative abundances across samples.
species_mean_depth_rel <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_rel_abun.tsv.gz',
                                     header = TRUE, row.names = 1, sep = '\t')

# Collapse all species with mean relative abundances < 1% into a category called 'Other'.
rare_species <- colnames(species_mean_depth_rel)[which(colMeans(species_mean_depth_rel) < 1)]
species_mean_depth_rel$Other <- rowSums(species_mean_depth_rel[, rare_species])
species_mean_depth_rel <- species_mean_depth_rel[, -which(colnames(species_mean_depth_rel) %in% rare_species)]


# Stacked barchart.
species_relabun_long <- reshape2::melt(as.matrix(species_mean_depth_rel))
species_relabun_long <- species_relabun_long[-which(species_relabun_long$value == 0), ]

colnames(species_relabun_long) <- c('Sample', 'Species', 'Relabun')

species_relabun_long$Sample <- factor(species_relabun_long$Sample, levels = rev(samples_ordered))

species_relabun_long$Species <- gsub('_', ' ', species_relabun_long$Species)
species_relabun_long$Species <- gsub('sp$', 'sp.', species_relabun_long$Species)
species_relabun_long$Species <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', species_relabun_long$Species)
species_relabun_long$Species <- factor(species_relabun_long$Species, levels = unique(species_relabun_long$Species))

plot_colours <- c("yellow1", # Bartonella apis

                  "#63bfaf", # Bifiobacteria
                  "#62d07a",

                  "#d56c24", # Bombilactobacilli
                  "#bf4d81",

                  "#d2cb3b", # Commensalibacter sp.

                  "brown", # Frischella perrara

                  "royalblue1", # Three Gilliamella
                  "skyblue1",
                  "mediumblue",
                  
                  "grey40", # Five Lactobacilli
                  "grey50",
                  "grey60",
                  "grey70",
                  "grey80", 
                  
                  "#9f4ad4", #Snodgrassella
                  
                  'black') # Other

italicized_labels <- c()

stacked_species_relabun <- ggplot(data = species_relabun_long, aes(y = Sample, x = Relabun, fill = Species)) +
                                  geom_bar(position='stack', stat='identity') +
                                  scale_fill_manual(values = plot_colours) +
                                  xlab('Species relative abundance (%)') +
                                  ylab('Sample') +
                                  coord_cartesian(expand = FALSE) +
                                  theme_bw() +
                                  theme(axis.text.y=element_blank(),
                                        axis.ticks = element_blank(),
                                        axis.title.y = element_text(angle = 0, vjust = 0.5),
                                        legend.text = element_text(face = 'italic'))

relabun_panel <- cowplot::plot_grid(NULL, stacked_species_relabun, NULL,
                                    nrow = 3, rel_heights = c(0.46, 8, 1.55))

overall_species_profiles <- cowplot::plot_grid(NULL,
                                               grid.grabExpr(draw(species_presence)),
                                               relabun_panel,
                                               rel_widths = c(0.03, 1, 1),
                                               nrow = 1,
                                               ncol = 3,
                                               labels = c('', 'a', 'b'))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_overall_species_profile.pdf',
       plot = overall_species_profiles,
       device = 'pdf',
       dpi = 600,
       width = 18,
       height = 9)
