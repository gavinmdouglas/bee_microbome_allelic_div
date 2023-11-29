rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)
library(poolr)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

results <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/within_sample_allele_permutation_results.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

results$species <- cleanup_species_names(results$species)

datasets <- sort(unique(results$dataset))
all_species <- sort(unique(results$species))

higher_tab <- data.frame(matrix(NA, nrow = length(all_species), ncol = length(datasets)))
rownames(higher_tab) <- all_species
colnames(higher_tab) <- datasets

lower_tab <- higher_tab

total_tab <- data.frame(matrix(NA, nrow = length(all_species), ncol = length(datasets)))
rownames(total_tab) <- all_species
colnames(total_tab) <- datasets

results <- results[which(results$num_alleles >= 5), ]
results <- results[which(results$num_samples >= 5), ]
results <- results[which(results$mean_alleles_per_sample >= 1.5), ]

results$p_lower_BH <- p.adjust(results$p_lower, 'BH')
results$p_higher_BH <- p.adjust(results$p_higher, 'BH')

for (sp in all_species) {
 
   for (d in datasets) {

     results_subset <- results[which(results$species == sp & results$dataset == d), ]
     
     if (nrow(results_subset) == 0) { next }

     lower_tab[sp, d] <- length(which(results_subset$p_lower_BH < 0.05))
     higher_tab[sp, d] <- length(which(results_subset$p_higher_BH < 0.05))

     total_tab[sp, d] <- nrow(results_subset)

   }
}

percent_lower <- (lower_tab / total_tab) * 100
percent_lower <- apply(percent_lower, 2, function(x) { round(x, 2)})

percent_higher <- (higher_tab / total_tab) * 100
percent_higher <- apply(percent_higher, 2, function(x) { round(x, 2)})

Heatmap(matrix = as.matrix(percent_higher),
        col = circlize::colorRamp2(c(0, 100), c('cornflowerblue', 'red')),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(percent_higher), '*', sep = '')),
        
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(higher_tab[i, j] > 0))
            grid.text(higher_tab[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })

Heatmap(matrix = as.matrix(percent_lower),
        col = circlize::colorRamp2(c(0, 100), c('cornflowerblue', 'red')),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(percent_lower), '*', sep = '')),
        
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(lower_tab[i, j] > 0))
            grid.text(lower_tab[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })
