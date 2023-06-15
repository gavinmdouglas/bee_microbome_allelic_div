rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

# Compare strains, genes, and alleles, based on basic characterizations and also to display relative variation across data types.
strain_basic_diversity <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_basic_diversity.tsv.gz',
                                     header = TRUE, sep = '\t', row.names = 1)

strain_and_gene_count_info <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_count_cor_summary.tsv.gz',
                                         header = TRUE, sep = '\t', row.names = 1)
strain_and_gene_count_info <- strain_and_gene_count_info[, -which(colnames(strain_and_gene_count_info) %in% c('cor_kendall', 'cor_p'))]

allele_count_info <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/accessory_allele_per_sample_summary.tsv.gz',
                                header = TRUE, sep = '\t', row.names = 1)
colnames(allele_count_info) <- paste('allele', colnames(allele_count_info), sep = '_')

jaccard_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE)
jaccard_summary_clean <- jaccard_summary
jaccard_summary$datatype <- tolower(jaccard_summary$datatype)

strain_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'strain'), -2]
rownames(strain_jaccard) <- strain_jaccard$species
strain_jaccard <- strain_jaccard[, -1]
colnames(strain_jaccard) <- paste('strain', colnames(strain_jaccard), sep = '_')

gene_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'gene'), -2]
rownames(gene_jaccard) <- gene_jaccard$species
gene_jaccard <- gene_jaccard[, -1]
colnames(gene_jaccard) <- paste('gene', colnames(gene_jaccard), sep = '_')

allele_jaccard <- jaccard_summary[which(jaccard_summary$datatype == 'allele'), -2]
rownames(allele_jaccard) <- allele_jaccard$species
allele_jaccard <- allele_jaccard[, -1]
colnames(allele_jaccard) <- paste('allele', colnames(allele_jaccard), sep = '_')

combined_data <- do.call(cbind, list(strain_basic_diversity,
                                     strain_and_gene_count_info,
                                     allele_count_info,
                                     strain_jaccard,
                                     gene_jaccard,
                                     allele_jaccard))

rownames(combined_data) <- gsub('_', ' ', rownames(combined_data))
rownames(combined_data) <- gsub(' sp$', ' sp.', rownames(combined_data))
rownames(combined_data) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', rownames(combined_data))

column_ordering <- c('total_samples_w_strains', # Strain profile-specific
                     'total_strains',
                     
                     'mean_strains_per_sample', # Mean instances per sample
                     'mean_acc.genes_per_sample',
                     'allele_mean_per_sample',
                     
                     'strain_mean_jaccard', # Mean Jaccard distance
                     'gene_mean_jaccard',
                     'allele_mean_jaccard')

combined_data <- combined_data[, column_ordering]
combined_data[is.na(combined_data)] <- NA

# Set CV values for species with low strains/samples to be NA.
combined_data[c('Serratia marcescens', 'Apilactobacillus kunkeei'), grep('^cv_', colnames(combined_data))] <- NA

combined_clean <- combined_data[, 3:ncol(combined_data)]
combined_clean <- format(round(combined_clean, 2), nsmall = 2)
combined_clean$mean_acc.genes_per_sample <- format(round(combined_data$mean_acc.genes_per_sample, 0), nsmall = 0)

combined_scaled <- as.matrix(scale(combined_data[, 3:ncol(combined_data)], center = TRUE, scale = TRUE))

# Get cluster in advance, to use same species order for both plots.
combined_scaled_hclust <- hclust(dist(combined_scaled, method = 'euclidean'), method = 'complete')

summary_heatmap <- Heatmap(matrix = combined_scaled,
                           
                           na_col = 'grey70',
                           col = circlize::colorRamp2(c(-4, 0, 4), c('blue', 'white', 'firebrick')),
                           heatmap_legend_param = list(title = 'Standard\nscore'),
                           show_heatmap_legend = TRUE,
                           
                           column_labels = c( 'Strains',
                                              'Accessory genes',
                                              'Accessory gene alleles',
                                              
                                              'Strains',
                                              'Accessory genes',
                                              'Accessory gene alleles'),
                           
                           cluster_rows = stats::as.dendrogram(combined_scaled_hclust),
                           cluster_columns = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_row_slices = FALSE,
                           
                           row_labels = gt_render(paste('*', rownames(combined_scaled), '*', sep = '')),
                           row_names_side = 'left',
                           row_dend_side = 'right',
                           column_names_rot = 45,
                           
                           column_split = factor(c(rep('Mean instances\nper sample', times = 3),
                                                   rep('Mean Jaccard\ndistances', times = 3)),
                                                 levels = c('Mean instances\nper sample',
                                                            'Mean Jaccard\ndistances')),
                           
                           column_gap = unit(5, 'mm'),
                           
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             if(! is.na(combined_clean[i, j] > 0))
                               grid.text(combined_clean[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                           })

strain_specific_info <- Heatmap(matrix = as.matrix(combined_data[, c(1, 2)]),
                           
                           col = circlize::colorRamp2(breaks = c(-10000, 10000),
                                                      colors = c('white', 'white')),
                           show_heatmap_legend = FALSE,
                           
                           column_labels = c('No. samples (with strains)',
                                             'No. strains'),
                           
                           cluster_rows = stats::as.dendrogram(combined_scaled_hclust),
                           cluster_columns = FALSE,
                           cluster_column_slices = FALSE,
                           cluster_row_slices = FALSE,
                           
                           row_labels = gt_render(paste('*', rownames(combined_scaled), '*', sep = '')),
                           show_row_dend = TRUE,
                           row_names_side = 'left',
                           column_names_rot = 45,
                           
                           column_split = rep('Strain\nprofile-specific', times = 2),
                           
                           column_gap = unit(5, 'mm'),
                           
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             if(! is.na(combined_data[, c(1, 2)][i, j] > 0))
                               grid.text(combined_data[, c(1, 2)][i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                           })

jaccard_summary_clean <- jaccard_summary_clean[which(jaccard_summary_clean$datatype != 'Species'), ]
jaccard_summary_clean$datatype <- factor(jaccard_summary_clean$datatype, levels = c('Strain', 'Gene', 'Allele'))

jaccard_boxplot <- ggplot(data = jaccard_summary_clean, aes(x = datatype, y = mean_jaccard)) +
                          geom_boxplot(outlier.shape = NA, fill = 'grey90') +
                          ggbeeswarm::geom_beeswarm(cex = 3) +
                          ggpubr::geom_pwc(label = 'P = {p}', method = 'wilcox.test') +
                          theme_bw() +
                          xlab('Data type') +
                          ylab('Mean pairwise\nJaccard distance\n(per species)') +
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0))


options(ggrepel.max.overlaps = Inf)

ppcor_gene_vs_allele_jaccard <- ppcor::pcor.test(x = combined_data$gene_mean_jaccard,
                                                 y = combined_data$allele_mean_jaccard,
                                                 z = combined_data$total_strains,
                                                 method = 'kendall')
ppcor_gene_vs_allele_jaccard_estimate <- format(round(ppcor_gene_vs_allele_jaccard$estimate, digits=4), nsmall = 4)
ppcor_gene_vs_allele_jaccard_p <- format(round(ppcor_gene_vs_allele_jaccard$p.value, digits=4), nsmall = 4)

ppcor_gene_vs_allele_jaccard_string <- paste0('Partial ~ tau ~', "'= '*",
                                              ppcor_gene_vs_allele_jaccard_estimate,
                                              "*','", '~ italic(P) ~', "'= '*",
                                              ppcor_gene_vs_allele_jaccard_p)

allele_vs_gene_scatterplot <- ggplot(data = combined_data, aes(x = gene_mean_jaccard,
                                                               y = allele_mean_jaccard,
                                                               colour = total_strains)) +
                                    geom_point() +
                                    ggrepel::geom_text_repel(aes(label = rownames(combined_data)),
                                                             col = 'grey75',
                                                             hjust = -0.1,
                                                             fontface = 'italic',
                                                             size = 1.5)  +
                                    annotate(geom = "text",
                                             y = 0.25,
                                             x = 0.40,
                                             label = ppcor_gene_vs_allele_jaccard_string,
                                             parse = TRUE) +
                                    theme_bw() +
                                    xlab('Gene mean pairwise Jaccard distances') +
                                    ylab('Allele mean\npairwise\nJaccard\ndistances') +
                                    scale_colour_gradient(limits = c(0, 50),
                                                          low = 'grey70',
                                                          high = 'black', 
                                                          name = 'No. strains') +
                                    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0),
                                          legend.position = c(0.80, 0.35),
                                          legend.background = element_rect(fill = 'white',
                                                                           colour = 'black'))
options(ggrepel.max.overlaps = NULL)


heatmap_panel <- plot_grid(grid.grabExpr(draw(strain_specific_info +
                                                summary_heatmap,
                                              main_heatmap = 2,
                                              ht_gap = unit(c(5), "mm"),
                                              padding = unit(c(2, 15, 2, 2), "mm"))),
                           labels = 'a')

bottom_row <- plot_grid(jaccard_boxplot, allele_vs_gene_scatterplot,
                        labels = c('b', 'c'), ncol = 2)

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_mgs_datatype_compare.pdf',
       plot = plot_grid(heatmap_panel,
                        bottom_row,
                        nrow = 2,
                        rel_heights = c(2, 1)),
       device = 'pdf',
       dpi = 600,
       width = 10,
       height = 12)
