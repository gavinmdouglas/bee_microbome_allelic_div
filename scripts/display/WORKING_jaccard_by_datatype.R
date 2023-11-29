rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

jaccard <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

jaccard_species_level <- jaccard[which(jaccard$species == 'All'), ]
jaccard <- jaccard[-which(jaccard$species == 'All'), ]

datasets <- sort(unique(jaccard$dataset))
datatypes <- c("StrainFacts", "StrainGST", "Gene", "Allele")
species <- sort(unique(jaccard$species))

jaccard_tab <- data.frame(matrix(NA, nrow = length(species), ncol = 20))
rownames(jaccard_tab) <- species

full_colnames <- character()
col_i <- 1
for (d_type in datatypes) {
  for (d_set in datasets) {
    col_name <- paste(d_type, d_set)
    full_colnames <- c(full_colnames, col_name)
    
    for (sp in species) {
      tab_subset <- jaccard[which(jaccard$species == sp & jaccard$datatype == d_type & jaccard$dataset == d_set), ]
      
      if (nrow(tab_subset) == 0) { 
        next
      } else {
        jaccard_tab[sp, col_i] <- tab_subset$mean_jaccard
      }
    }
    col_i <- col_i + 1
  }
}

colnames(jaccard_tab) <- full_colnames

datasets_clean <- c('Ellegaard 2019',
                    'Ellegaard 2020',
                    'Sun 2022',
                    'Wu 2021',
                    'Zhang 2022')

jaccard_tab_rounded <- jaccard_tab
jaccard_tab_rounded <- format(round(jaccard_tab_rounded, 1), nsmall = 1)
jaccard_tab_rounded[jaccard_tab_rounded == ' NA'] <- ''

column_split_labels <- c(rep('StrainFacts', times = 5),
                         rep('StrainGST', times = 5),
                         rep('Genes', times = 5),
                         rep('Alleles', times = 5))
column_split_labels <- factor(column_split_labels, c('StrainFacts', 'StrainGST', 'Genes', 'Alleles'))

overall_jaccard <- Heatmap(matrix = as.matrix(jaccard_tab),
                                    
                                    na_col = 'grey70',
                                    
                                    col = circlize::colorRamp2(c(0, 1), c('white', 'firebrick')),
                                    
                                    heatmap_legend_param = list(title = 'Mean\nJaccard'),
                                    
                                    show_heatmap_legend = TRUE,
                                    
                                    column_labels = rep(datasets_clean, 4),
                                    
                                    cluster_rows = TRUE,
                                    cluster_columns = FALSE,
                                    cluster_column_slices = FALSE,
                                    cluster_row_slices = FALSE,
                                    
                                    row_labels = gt_render(paste('*', rownames(jaccard_tab), '*', sep = '')),
                                    row_names_side = 'left',
                                    row_dend_side = 'right',
                                    column_names_rot = 45,
                                    
                                    column_split = column_split_labels,
                           
                                    
                                    column_gap = unit(5, 'mm'),
                                    
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if(! is.na(jaccard_tab_rounded[i, j] > 0))
                                        grid.text(jaccard_tab_rounded[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                                    })
