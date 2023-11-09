rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

# Compare strains, genes, and alleles, based on basic characterizations and also to display relative variation across data types.
strain_basic_diversity <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_diversity.tsv.gz',
                                     header = TRUE, sep = '\t')
rownames(strain_basic_diversity) <- paste(strain_basic_diversity$species, strain_basic_diversity$dataset)

strain_and_gene_count_info <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_count_cor_summary.tsv.gz',
                                         header = TRUE, sep = '\t')
strain_and_gene_count_info <- strain_and_gene_count_info[, -which(colnames(strain_and_gene_count_info) %in% c('cor_kendall', 'cor_p'))]
rownames(strain_and_gene_count_info) <- paste(strain_and_gene_count_info$species, strain_and_gene_count_info$dataset)
strain_and_gene_count_info <- strain_and_gene_count_info[, -c(1:2)]

allele_count_info <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/accessory_allele_per_sample_summary.tsv.gz',
                                header = TRUE, sep = '\t')
rownames(allele_count_info) <- paste(allele_count_info$species, allele_count_info$dataset)
allele_count_info <- allele_count_info[, -c(1:2)]
colnames(allele_count_info) <- paste('allele', colnames(allele_count_info), sep = '_')

intersecting_rows <- intersect(rownames(strain_basic_diversity), rownames(strain_and_gene_count_info))
intersecting_rows <- intersect(intersecting_rows, rownames(allele_count_info))

combined_data <- do.call(cbind, list(strain_basic_diversity[intersecting_rows, ],
                                     strain_and_gene_count_info[intersecting_rows, ],
                                     allele_count_info[intersecting_rows, ]))
combined_data <- combined_data[order(combined_data$dataset, combined_data$species), ]

rownames(combined_data) <- gsub('_', ' ', rownames(combined_data))
rownames(combined_data) <- gsub(' sp ', ' sp. ', rownames(combined_data))
rownames(combined_data) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', rownames(combined_data))

column_ordering <- c('num_intersecting_samples', # Strain profile-specific
                     'strainfacts_total_strains',
                     'straingst_total_strains',
                     
                     'strainfacts_mean_richness', # Mean instances per sample
                     'straingst_mean_richness',
                     'mean_acc.genes_per_sample',
                     'allele_mean_per_sample')

combined_clean <- combined_data[, column_ordering]
combined_clean[is.na(combined_clean)] <- NA

combined_num_instances <- combined_clean[, 4:ncol(combined_clean)]

combined_strain_specific <- combined_clean[, 1:3]

species_vec <- sapply(strsplit(rownames(combined_num_instances), ' '), function(x) { return(paste(x[1], x[2], sep = ' '))})
dataset_vec <- sapply(strsplit(rownames(combined_num_instances), ' '), function(x) { return(x[3])})
dataset_vec <- sub('2', ' 2', dataset_vec)


# Scale by dataset and column.
combined_num_instances_scaled <- matrix(NA, nrow = nrow(combined_num_instances), ncol(combined_num_instances))
combined_strain_specific_scale <- matrix(NA, nrow = nrow(combined_strain_specific), ncol(combined_strain_specific))

for (d in unique(dataset_vec)) {
  d_i <- which(dataset_vec == d)
  combined_num_instances_scaled[d_i, ] <- scale(combined_num_instances[d_i, ], center = TRUE, scale = TRUE)
  combined_strain_specific_scale[d_i, ] <- scale(combined_strain_specific[d_i, ], center = TRUE, scale = TRUE)
}

combined_num_instances <- format(round(combined_num_instances, 2), nsmall = 2)
combined_num_instances$mean_acc.genes_per_sample <- format(round(combined_data$mean_acc.genes_per_sample, 0), nsmall = 0)


# Get cluster in advance, to use same species order for both plots.
# combined_num_instances_scaled_hclust <- hclust(dist(combined_num_instances_scaled, method = 'euclidean'), method = 'complete')

summary_heatmap <- Heatmap(matrix = combined_num_instances_scaled,
                           
                           na_col = 'grey70',
                           col = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'firebrick')),
                           heatmap_legend_param = list(title = 'Standard\nscore\n(by\ncolumn\nand\ndataset)'),
                           show_heatmap_legend = FALSE,
  
                           column_labels = c( 'StrainFacts',
                                              'StrainGST',
                                              'Accessory genes',
                                              'Accessory gene alleles'),
                           
                           cluster_rows = FALSE,
                           cluster_columns = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_row_slices = FALSE,
                           
                           row_split = dataset_vec,
                           
                           row_labels = gt_render(paste('*', species_vec, '*', sep = '')),
                           row_names_side = 'left',
                           row_dend_side = 'right',
                           column_names_rot = 45,
                           
                           column_split = rep('Mean instances\nper sample', times = 4),
                           
                           row_gap = unit(5, 'mm'),
                           
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             if(! is.na(combined_num_instances[i, j] > 0))
                               grid.text(combined_num_instances[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                           })


strain_specific_info <- Heatmap(matrix = combined_strain_specific_scale,
                                
                                col = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'firebrick')),
                                heatmap_legend_param = list(title = 'Standard\nscore\n(by\ncolumn\nand\ndataset)'),
                                show_heatmap_legend = TRUE,
                                
                                column_labels = c('No. samples (with strains)',
                                                  'No. StrainFacts strains',
                                                  'No. StrainGST strains'),

                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                cluster_column_slices = FALSE,
                                cluster_row_slices = FALSE,
                                
                                row_split = dataset_vec,
                                
                                row_labels = gt_render(paste('*', species_vec, '*', sep = '')),
                                show_row_dend = FALSE,
                                row_names_side = 'left',
                                column_names_rot = 45,
                                
                                column_split = rep('Strain\nprofile-specific', times = 3),
                                
                                row_gap = unit(5, 'mm'),
                                
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if(! is.na(combined_strain_specific[i, j] > 0))
                                    grid.text(combined_strain_specific[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                })


combined_heatmap <- plot_grid(grid.grabExpr(draw(strain_specific_info +
                                                 summary_heatmap,
                                                 main_heatmap = 1,
                                                 ht_gap = unit(c(5), "mm"),
                                                 padding = unit(c(2, 15, 2, 2), "mm"))))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_mgs_datatype_compare.pdf',
       plot = combined_heatmap,
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 12)
