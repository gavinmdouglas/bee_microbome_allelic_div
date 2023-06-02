rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(cowplot)
library(ggplot2)

# Make summary plot with this information:
  # Species names
  # Higher-level taxonomic lineages
  # Gram stain
  # Core phylotype
  # Number of genomes
  # Number of genes
  # Percentage core genes

# All of this information will be plotted as individual heatmaps, to make them easier to visually parse.

# Read in species information.
species_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/db_species_summary_info.txt.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

species <- rownames(species_info)

rownames(species_info) <- gsub('_', ' ', rownames(species_info))
rownames(species_info) <- gsub(' sp$', ' sp.', rownames(species_info))
rownames(species_info) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium coryneforme/indicum', rownames(species_info))

# Sort by higher lineage, and then by species name.
species_info <- species_info[order(species_info$Phylum_Class_Order_Family, rownames(species_info)), ]

# Prep higher lineage heatmap.
higher_lineage <- species_info[, 'Phylum_Class_Order_Family', drop = FALSE]

higher_lineage_dummy_colour <- rep('grey95', times = length(unique(higher_lineage$Phylum_Class_Order_Family)))
names(higher_lineage_dummy_colour) <- unique(higher_lineage$Phylum_Class_Order_Family)

higher_lineage_heatmap <- Heatmap(matrix = as.matrix(higher_lineage),
                                    
                                    col = higher_lineage_dummy_colour,
                                    
                                    show_heatmap_legend = FALSE,  
                                    
                                    column_labels = 'Phylum; Class; Order; Family',
                                    cluster_rows = FALSE,
                                    show_column_names = TRUE,
                                    column_names_rot = 45,
                                    cluster_columns = FALSE,
                                    row_names_side = 'left',
                                    row_labels = gt_render(paste('*', rownames(higher_lineage), '*', sep = '')),
                                    
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if(! is.na(higher_lineage[i, j] > 0))
                                        grid.text(higher_lineage[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                    },
                                  width = 2.2)


# Prep Gram_stain heatmap.
Gram_stain <- species_info[, 'Gram_stain', drop = FALSE]
Gram_stain$Gram_stain[which(Gram_stain$Gram_stain == 'Negative')] <- '-'
Gram_stain$Gram_stain[which(Gram_stain$Gram_stain == 'Positive')] <- '+'

Gram_stain_colours <- c('#6a1890', '#e130b2')
names(Gram_stain_colours) <- c('+', '-')

Gram_stain_heatmap <- Heatmap(matrix = as.matrix(Gram_stain),
                                      
                                      col = Gram_stain_colours,
                                      
                                      show_heatmap_legend = FALSE,  
                                      
                                      column_labels = 'Gram stain',
                                      cluster_rows = FALSE,
                                      show_column_names = TRUE,
                                      column_names_rot = 45,
                                      cluster_columns = FALSE,
                                      row_names_side = 'left',
                                      row_labels = gt_render(paste('*', rownames(Gram_stain), '*', sep = '')),
                                      
                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                        if(! is.na(Gram_stain[i, j] > 0))
                                          grid.text(Gram_stain[i, j], x, y, gp = gpar(fontsize = 10, col = 'white'), just = 'centre')
                                      },
                              width = 0.1)

# Prep Phylotype heatmap.
Phylotype <- species_info[, 'Core_phylotype', drop = FALSE]

Phylotype_colours <- c('#63bfaf', '#d56c24', 'yellow3', 'skyblue1', '#C9A0DC', 'grey95')
names(Phylotype_colours) <- c('Bifidobacterium', 'Firm-4', 'Firm-5', 'Gilliamella', 'Snodgrassella', 'Non-core phylotype')

Phylotype_heatmap <- Heatmap(matrix = as.matrix(Phylotype),
                                  
                                  col = Phylotype_colours,
                                  
                                  show_heatmap_legend = FALSE,  
                                  
                                  column_labels = 'Core phylotype',
                                  cluster_rows = FALSE,
                                  show_column_names = TRUE,
                                  column_names_rot = 45,
                                  cluster_columns = FALSE,
                                  row_names_side = 'left',
                                  row_labels = gt_render(paste('*', rownames(Phylotype), '*', sep = '')),
                                  
                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                    if(! is.na(Phylotype[i, j] > 0))
                                      grid.text(Phylotype[i, j], x, y, gp = gpar(fontsize = 10, col = 'black'), just = 'centre')
                                  },
                             width = 0.7)

# Prep genome counts heatmap.
genome_counts <- data.frame(Count = rep(NA, times = nrow(species_info)))
rownames(genome_counts) <- species
for (sp in species) {
  genome_counts[sp, 'Count'] <- length(read.table(file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_genome_sets/', sp, '.txt.gz', sep = ''),
                                                  stringsAsFactors = FALSE)$V1)
}

rownames(genome_counts) <- gsub('_', ' ', rownames(genome_counts))
rownames(genome_counts) <- gsub(' sp$', ' sp.', rownames(genome_counts))
rownames(genome_counts) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium coryneforme/indicum', rownames(genome_counts))
genome_counts <- genome_counts[rownames(species_info), , drop = FALSE]


genome_count_heatmap <- Heatmap(matrix = as.matrix(genome_counts),

                          col = circlize::colorRamp2(c(0, 35), c("grey95", "firebrick2")),
                          
                          show_heatmap_legend = FALSE,  
                          
                          column_labels = 'Number of genomes',
                          cluster_rows = FALSE,
                          show_column_names = TRUE,
                          column_names_rot = 45,
                          cluster_columns = FALSE,
                          row_names_side = 'left',
                          row_labels = gt_render(paste('*', rownames(genome_counts), '*', sep = '')),
                  
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            if(! is.na(genome_counts[i, j] > 0))
                              grid.text(genome_counts[i, j], x, y, gp = gpar(fontsize = 10))
                          },
                          width = 0.3)


# Final gene count, and percentage core genes.
gene_counts <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_gene_families_per_species.tsv.gz',
                          header = TRUE, sep = '\t', row.names = 1)

# Convert core counts to percentages.
gene_counts$core <- 100 * (gene_counts$core / gene_counts$all)

colnames(gene_counts) <- c('No. gene families', '% Core gene families')
rownames(gene_counts) <- gsub('_', ' ', rownames(gene_counts))
rownames(gene_counts) <- gsub(' sp$', ' sp.', rownames(gene_counts))
rownames(gene_counts) <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium coryneforme/indicum', rownames(gene_counts))

gene_counts <- gene_counts[rownames(species_info), ]

gene_counts_only <- gene_counts[, 'No. gene families', drop = FALSE]
gene_core_percent <- gene_counts[, '% Core gene families', drop = FALSE]
gene_core_percent$`% Core gene families` <- round(gene_core_percent$`% Core gene families`, 2)

gene_count_heatmap <- Heatmap(matrix = as.matrix(gene_counts_only),
                                    
                                    col = circlize::colorRamp2(c(0, 5200), c("grey95", "#1976D2")),
                                    
                                    show_heatmap_legend = FALSE,  
                                  
                                    column_labels = 'Number of gene families',
                                    cluster_rows = FALSE,
                                    show_column_names = TRUE,
                                    column_names_rot = 45,
                                    cluster_columns = FALSE,
                                    row_names_side = 'left',
                                    row_labels = gt_render(paste('*', rownames(gene_counts_only), '*', sep = '')),
                                    
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if(! is.na(gene_counts_only[i, j] > 0))
                                        grid.text(gene_counts_only[i, j], x, y, gp = gpar(fontsize = 10))
                                    },
                              width = 0.3)

core_percent_heatmap <- Heatmap(matrix = as.matrix(gene_core_percent),
                              
                              col = circlize::colorRamp2(c(0, 100), c("grey95", "seagreen3")),
                              
                              show_heatmap_legend = FALSE,  
                              
                              column_labels = 'Percent core gene families',
                              cluster_rows = FALSE,
                              show_column_names = TRUE,
                              column_names_rot = 45,
                              cluster_columns = FALSE,
                              row_names_side = 'left',
                              row_labels = gt_render(paste('*', rownames(gene_core_percent), '*', sep = '')),
                              
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(! is.na(gene_core_percent[i, j] > 0))
                                  grid.text(gene_core_percent[i, j], x, y, gp = gpar(fontsize = 10))
                              },
                              width = 0.3)

combined_heatmap <- plot_grid(grid.grabExpr(draw(higher_lineage_heatmap +
                                                 Gram_stain_heatmap +
                                                 Phylotype_heatmap +
                                                 genome_count_heatmap +
                                                 gene_count_heatmap +
                                                 core_percent_heatmap,
                                                 padding = unit(c(2, 15, 2, 2), "mm"))))

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_ref_db_species_summary.pdf',
       plot = combined_heatmap,
       device = 'pdf',
       dpi = 600,
       width = 12,
       height = 8)
