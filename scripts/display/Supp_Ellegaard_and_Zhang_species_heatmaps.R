rm(list = ls(all.names = TRUE))

# Heatmap of species presence/absence overlaid with detailed metadata from Ellegaard and Zhang datasets.

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)

species_90per <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz',
                            row.names = 1, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

colnames(species_90per) <- gsub('_', ' ', colnames(species_90per))
colnames(species_90per) <- gsub('sp$', 'sp.', colnames(species_90per))
colnames(species_90per) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', colnames(species_90per))

species_90per_char <- species_90per
species_90per_char[species_90per_char == 0] <- 'Absent'
species_90per_char[species_90per_char == '1'] <- 'Present'

Ellegaard2019_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Ellegaard2019_metadata.tsv.gz',
                                     header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
Ellegaard2020_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Ellegaard2020_metadata.tsv.gz',
                                     header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
Ellegaard_combined_metadata <- rbind(Ellegaard2019_metadata, Ellegaard2020_metadata)
Ellegaard_combined_metadata$Year <- as.factor(Ellegaard_combined_metadata$Year)

Ellegaard_Country_colours <- c('#7A942E','#D52B1E')
names(Ellegaard_Country_colours) <- c('Japan', 'Switzerland')

Ellegaard_Apiary_colours <- c('#67bda7', '#72c365', '#cf5930', '#b1504f')
names(Ellegaard_Apiary_colours) <- c("Ai", "Iu", "Les Droites", "Grammont")

Ellegaard_Year_colours <- c('#E6E6FA', '#c15498', '#7e8ec3')
names(Ellegaard_Year_colours) <- c('2015', '2016', '2020')

Ellegaard_Age_colours <- c('#c6a93c', '#ada1a1', '#4a3548')
names(Ellegaard_Age_colours) <- c('Young', 'Middle-aged', 'Old')

Ellegaard_column_annot <- ComplexHeatmap::HeatmapAnnotation(Country = Ellegaard_combined_metadata$Country,
                                                            Apiary = factor(Ellegaard_combined_metadata$Apiary,
                                                                            levels = c("Ai", "Iu", "Les Droites", "Grammont")),
                                                            Year = Ellegaard_combined_metadata$Year,
                                                            Age = factor(Ellegaard_combined_metadata$Age,
                                                                         levels = c('Young', 'Middle-aged', 'Old')),
                                                            col = list(Country = Ellegaard_Country_colours,
                                                                       Apiary = Ellegaard_Apiary_colours,
                                                                       Year = Ellegaard_Year_colours,
                                                                       Age = Ellegaard_Age_colours))

Zhang_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Zhang2022_metadata.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
Zhang_metadata$Day <- gsub('Day ', '', Zhang_metadata$Day)
Zhang_metadata$Day <- factor(Zhang_metadata$Day, levels = c('7', '11', '19'))
Zhang_metadata$Treatment <- factor(Zhang_metadata$Treatment,
                                   levels = c('Control', 'Antibiotic'))

Zhang_Treatment_colours <- c('grey75', 'black')
names(Zhang_Treatment_colours) <- c('Control', 'Antibiotic')

Zhang_Day_colours <- c('#7aa457', '#9e6ebd', '#cb6751')
names(Zhang_Day_colours) <- c('7', '11', '19')

Zhang_column_annot <- ComplexHeatmap::HeatmapAnnotation(Day=Zhang_metadata$Day,
                                                        Treatment=Zhang_metadata$Treatment,
                                                        col = list(Treatment = Zhang_Treatment_colours,
                                                                   Day = Zhang_Day_colours))

## Start by performing clustering within each of these datasets (Ellegaard datasets being considered together).
species_90per_Ellegaard <- species_90per[rownames(Ellegaard_combined_metadata), ]
species_90per_Ellegaard_cluster <- stats::hclust(dist(t(species_90per_Ellegaard), method = 'binary'), method = 'complete')
sample_90per_Ellegaard_cluster <- stats::hclust(dist(species_90per_Ellegaard, method = 'binary'), method = 'complete')

species_90per_Ellegaard_order <- species_90per_Ellegaard_cluster$labels[species_90per_Ellegaard_cluster$order]
samples_90per_Ellegaard_order <- sample_90per_Ellegaard_cluster$labels[sample_90per_Ellegaard_cluster$order]

species_Ellegaard_heatmap <- ComplexHeatmap::Heatmap(t(species_90per_char[rownames(Ellegaard_combined_metadata), ]),
                                                     col = c('grey95', 'cornflowerblue'),
                                                     heatmap_legend_param = list(title = 'Species'),
                                                     row_labels = ComplexHeatmap::gt_render(paste0("*", colnames(species_90per_char), "*")),
                                                     cluster_columns = sample_90per_Ellegaard_cluster,
                                                     cluster_rows = species_90per_Ellegaard_cluster,
                                                     show_row_names = TRUE,
                                                     show_column_names = TRUE,
                                                     top_annotation = Ellegaard_column_annot)

species_90per_Zhang <- species_90per[rownames(Zhang_metadata), ]
species_90per_Zhang_cluster <- stats::hclust(dist(t(species_90per_Zhang), method = 'binary'), method = 'complete')
sample_90per_Zhang_cluster <- stats::hclust(dist(species_90per_Zhang, method = 'binary'), method = 'complete')

species_90per_Zhang_order <- species_90per_Zhang_cluster$labels[species_90per_Zhang_cluster$order]
samples_90per_Zhang_order <- sample_90per_Zhang_cluster$labels[sample_90per_Zhang_cluster$order]

species_Zhang_heatmap <- Heatmap(t(species_90per_char[rownames(Zhang_metadata), ]),
                                 col = c('grey95', 'cornflowerblue'),
                                 heatmap_legend_param = list(title = 'Species'),
                                 row_labels = ComplexHeatmap::gt_render(paste0("*", colnames(species_90per_char), "*")),
                                 cluster_columns = sample_90per_Zhang_cluster,
                                 cluster_rows = species_90per_Zhang_cluster,
                                 show_row_names = TRUE,
                                 show_column_names = TRUE,
                                 top_annotation = Zhang_column_annot)

Ellegaard_and_Zhang_heatmaps <- cowplot::plot_grid(grid.grabExpr(draw(species_Ellegaard_heatmap, column_title = 'Ellegaard 2019 and 2020 datasets')),
                                                   grid.grabExpr(draw(species_Zhang_heatmap, column_title = 'Zhang 2022 dataset')),
                                                   nrow = 2,
                                                   ncol = 1,
                                                   labels = c('a', 'b'))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_Ellegaard_and_Zhang_species_heatmap.pdf',
       plot = Ellegaard_and_Zhang_heatmaps,
       device = 'pdf',
       dpi = 600,
       width = 15,
       height = 12)



# Also run PERMANOVA on the Jaccard distance of these matrices, along with the highlighted metadata.
# May want to summarize alongside.
library(vegan)

Ellegaard_sp_jaccard <- species_90per[rownames(Ellegaard_combined_metadata), ]
Ellegaard_sp_adonis <- adonis2(Ellegaard_sp_jaccard ~ Ellegaard_combined_metadata$Country + Ellegaard_combined_metadata$Apiary + Ellegaard_combined_metadata$Year + Ellegaard_combined_metadata$Age)

Zhang_sp_jaccard <- species_90per[rownames(Zhang_metadata), ]
Zhang_sp_adonis <- adonis2(Zhang_sp_jaccard ~ Zhang_metadata$Treatment + Zhang_metadata$Day)


# Print the results
print(result)
