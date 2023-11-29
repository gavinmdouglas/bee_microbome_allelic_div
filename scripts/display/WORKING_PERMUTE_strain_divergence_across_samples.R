rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)
library(poolr)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

results <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/within_sample_strain_permutation_results.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

results$species <- cleanup_species_names(results$species)

datasets <- sort(unique(results$dataset))
all_species <- sort(unique(results$species))

sig_tab <- data.frame(matrix(NA, nrow = length(all_species), ncol = length(datasets)))
rownames(sig_tab) <- all_species
colnames(sig_tab) <- datasets

OR_tab <- data.frame(matrix(NA, nrow = length(all_species), ncol = length(datasets)))
rownames(OR_tab) <- all_species
colnames(OR_tab) <- datasets


results <- results[which(results$num_strains >= 3), ]
results <- results[which(results$num_samples >= 3), ]
results <- results[which(results$mean_strains_per_sample >= 1.5), ]

for (sp in all_species) {
 
   for (d in datasets) {

     row_i <- which(results$species == sp & results$dataset == d)
     
     if (length(row_i) == 0) {
       sig_tab[sp, d] <- 'Insufficient data'
       next
     }

     if (results[row_i, 'p_lower'] < 0.05) {
       sig_tab[sp, d] <- 'Lower'
     } else if (results[row_i, 'p_higher'] < 0.05) {
       sig_tab[sp, d] <- 'Higher'
     } else {
       sig_tab[sp, d] <- 'Not significant'
     }
     
     OR_tab[sp, d] <- results[row_i, 'obs_mean'] / results[row_i, 'permuted_mean']
     
   }
}

OR_tab <- apply(OR_tab, 2, function(x) { format(round(x, 4), nsmall = 4) })

OR_tab[OR_tab == "    NA"] <- ''

Heatmap(matrix = as.matrix(sig_tab),
        col = c('red', 'grey', 'grey95'),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(sig_tab), '*', sep = '')),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(OR_tab[i, j] > 0))
            grid.text(OR_tab[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })


sig_tab_combined <- data.frame(matrix(NA, nrow = length(all_species), ncol = 1))
rownames(sig_tab_combined) <- all_species
colnames(sig_tab_combined) <- 'combined'

for (sp in all_species) {
    
    tab_subset <- results[which(results$species == sp), ]
    
    if (nrow(tab_subset) == 0) { next }
    
    combined_p <- poolr::fisher(p = tab_subset$p_higher)$p
    
    if (combined_p < 0.05) {
      sig_tab_combined[sp, 1] <- 'Significant'
    } else {
      sig_tab_combined[sp, 1] <- 'Non-significant'
    }
}

sig_tab_combined[is.na(sig_tab_combined)] <- 'Insufficient data'

Heatmap(matrix = as.matrix(sig_tab_combined),
        col = c('grey', 'grey95', 'red'),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(sig_tab_combined), '*', sep = '')))
