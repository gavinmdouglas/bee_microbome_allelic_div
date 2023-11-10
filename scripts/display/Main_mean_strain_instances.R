rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

# Compare strains, genes, and alleles, based on basic characterizations and also to display relative variation across data types.
strain_diversity <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_diversity.tsv.gz',
                                     header = TRUE, sep = '\t')
strain_diversity$species <- cleanup_species_names(names_vec = strain_diversity$species, shorten=FALSE)
strain_diversity$dataset <- sub('2', ' 2', strain_diversity$dataset)

strainfacts_richness <- reshape2::dcast(data = strain_diversity,
                                        formula = species ~ dataset,
                                        fill = NA,
                                        value.var = 'strainfacts_mean_richness')
rownames(strainfacts_richness) <- strainfacts_richness$species
strainfacts_richness <- strainfacts_richness[, which(colnames(strainfacts_richness) != 'species')]

straingst_richness <- reshape2::dcast(data = strain_diversity,
                                        formula = species ~ dataset,
                                        fill = NA,
                                        value.var = 'straingst_mean_richness')
rownames(straingst_richness) <- straingst_richness$species
straingst_richness <- straingst_richness[, which(colnames(straingst_richness) != 'species')]

combined_richness <- cbind(strainfacts_richness, straingst_richness)

combined_richness_scaled <- scale(combined_richness, center = TRUE, scale = TRUE)

combined_richness_rounded <- format(round(combined_richness, 1), nsmall = 1)
combined_richness_rounded[combined_richness_rounded == ' NA'] <- ''

heatmap_mean_num_strains <- Heatmap(matrix = combined_richness_scaled,
                           
                           na_col = 'grey70',
                           
                           col = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'firebrick')),
                           
                           heatmap_legend_param = list(title = 'Standard\nscore\n(by\ncolumn)'),
                           
                           show_heatmap_legend = TRUE,
                           
                           column_labels = c('Ellegaard 2019',
                                             'Ellegaard 2020',
                                             'Sun 2022',
                                             'Wu 2021',
                                             'Zhang 2022',
                                             'Ellegaard 2019',
                                             'Ellegaard 2020',
                                             'Sun 2022',
                                             'Wu 2021',
                                             'Zhang 2022'),
                           
                           cluster_rows = TRUE,
                           cluster_columns = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_row_slices = FALSE,
                           
                           row_labels = gt_render(paste('*', rownames(combined_richness_rounded), '*', sep = '')),
                           row_names_side = 'left',
                           row_dend_side = 'right',
                           column_names_rot = 45,
                           
                           column_split = c(rep('StrainFacts', times = 5),
                                            rep('StrainGST', times = 5)),
                           
                           column_gap = unit(5, 'mm'),
                           
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             if(! is.na(combined_richness_rounded[i, j] > 0))
                               grid.text(combined_richness_rounded[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                           })


heatmap_mean_num_strains_parsed <- plot_grid(grid.grabExpr(draw(heatmap_mean_num_strains)))
                                                 #padding = unit(c(2, 15, 2, 2), "mm"))))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_mean_num_strains.pdf',
       plot = heatmap_mean_num_strains_parsed,
       device = 'pdf',
       dpi = 600,
       width = 8,
       height = 7)
