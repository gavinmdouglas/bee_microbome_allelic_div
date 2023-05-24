rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)

# Plot number of genomes and genes per species.
species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                      stringsAsFactors = FALSE)$V1

genome_counts <- data.frame(Count = rep(NA, times = length(species)))
rownames(genome_counts) <- species

for (sp in species) {
  genome_counts[sp, 'Count'] <- length(read.table(file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_genome_sets/', sp, '.txt.gz', sep = ''),
                                                  stringsAsFactors = FALSE)$V1)
}
rownames(genome_counts) <- gsub('_', ' ', rownames(genome_counts))
rownames(genome_counts) <- gsub(' sp$', ' sp.', rownames(genome_counts))
rownames(genome_counts) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium coryneforme/indicum', rownames(genome_counts))


gene_counts <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/orthologs_per_species.tsv.gz',
                          header = FALSE, sep = '\t', row.names = 1)
colnames(gene_counts) <- 'Count'
rownames(gene_counts) <- gsub('_', ' ', rownames(gene_counts))
rownames(gene_counts) <- gsub(' sp$', ' sp.', rownames(gene_counts))
rownames(gene_counts) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium coryneforme/indicum', rownames(gene_counts))

genome_count_heatmap_raw <- Heatmap(matrix = as.matrix(genome_counts),

                          col = circlize::colorRamp2(c(0, 35), c("grey95", "firebrick2")),
                          
                          heatmap_legend_param = list(title = 'Number of\ngenomes'),
                         
                          cluster_rows = TRUE,
                          show_row_dend = FALSE,
                          show_column_names = FALSE,
                          cluster_columns = FALSE,
                          row_names_side = 'left',
                          row_labels = gt_render(paste('*', rownames(genome_counts), '*', sep = '')),
                  
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(! is.na(genome_counts[i, j] > 0))
                              grid.text(genome_counts[i, j], x, y, gp = gpar(fontsize = 10))
                          })


gene_count_heatmap_raw <- Heatmap(matrix = as.matrix(gene_counts),
                                    
                                    col = circlize::colorRamp2(c(0, 6500), c("grey95", "#1976D2")),
                                    
                                    heatmap_legend_param = list(title = 'Number of\northologs'),
                                    
                                    cluster_rows = TRUE,
                                    show_row_dend = FALSE,
                                    show_column_names = FALSE,
                                    cluster_columns = FALSE,
                                    row_names_side = 'left',
                                    row_labels = gt_render(paste('*', rownames(gene_counts), '*', sep = '')),
                                    
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if(! is.na(gene_counts[i, j] > 0))
                                        grid.text(gene_counts[i, j], x, y, gp = gpar(fontsize = 10))
                                    })

count_heatmap <- plot_grid(grid.grabExpr(draw(genome_count_heatmap_raw + gene_count_heatmap_raw,
                                              padding = unit(c(2, 15, 2, 2), "mm"))))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_number_of_ref_genomes.pdf',
       plot = count_heatmap,
       device = 'pdf',
       dpi = 600,
       width = 4.7,
       height = 8)
