rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

within_div <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/within_sample_strain_divergence.tsv.gz',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)

between_div <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/between_sample_strain_divergence.tsv.gz',
                          header = TRUE, sep = '\t', stringsAsFactors = FALSE)

within_div_clean <- within_div
within_div_clean <- within_div_clean[, -1]
colnames(within_div_clean)[1] <- 'mean_div'
within_div_clean$Comparison <- 'Within'

between_div_clean <- between_div
between_div_clean <- between_div_clean[, -c(1, 2)]
colnames(between_div_clean)[1] <- 'mean_div'
between_div_clean$Comparison <- 'Between'

combined_strain_div <- rbind(within_div_clean, between_div_clean)

combined_strain_div$species <- cleanup_species_names(combined_strain_div$species)

ggplot(data = combined_strain_div, aes(x = mean_div,
                                       y = species,
                                       fill = Comparison)) +
  geom_boxplot() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.y = element_text(face = "italic")) +
  ylab('Species') +
  xlab('Percent identity')

combined_strain_div_no.outlier <- combined_strain_div[which(combined_strain_div$mean_div > 90), ]

ggplot(
  data = combined_strain_div_no.outlier,
  aes(x = mean_div,
      y = species,
      fill = Comparison)) +
  geom_boxplot() +
  scale_y_discrete(limits = rev) +
  theme(axis.text.y = element_text(face = "italic")) +
  ylab('Species') +
  xlab('Percent identity')

datasets <- sort(unique(combined_strain_div$dataset))
all_species <- sort(unique(combined_strain_div$species))



sig_tab <- data.frame(matrix(NA, nrow = length(all_species), ncol = length(datasets)))
rownames(sig_tab) <- all_species
colnames(sig_tab) <- datasets

for (sp in all_species) {
 
   for (d in datasets) {

      tab_subset <- combined_strain_div[which(combined_strain_div$species == sp & combined_strain_div$dataset == d), ]
      
      within_vec <- tab_subset[which(tab_subset$Comparison == 'Within'), 'mean_div']
      between_vec <- tab_subset[which(tab_subset$Comparison == 'Between'), 'mean_div']
      
      if (length(within_vec) < 5 || length(between_vec) < 5) {
        sig_tab[sp, d] <- 'Insufficient data'
      } else {
        wilcox_out <- wilcox.test(within_vec, between_vec, exact = FALSE)
        
        wilcox_p <- wilcox_out$p.value
        
        if (is.na(wilcox_p)) {
          wilcox_p <- 1 
        }
        
        if (wilcox_p >= 0.05) {
          sig_tab[sp, d] <- 'Non-significant'
        } else if (mean(within_vec) > mean(between_vec)) {
          sig_tab[sp, d] <- 'Higher within'
        } else if (mean(within_vec) < mean(between_vec)) {
          sig_tab[sp, d] <- 'Higher between'
        } else {
          stop('Error!') 
        }
      }
   }
}

Heatmap(matrix = as.matrix(sig_tab),
        col = c('red', 'blue', 'grey', 'grey95'),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(sig_tab), '*', sep = '')))


sig_tab_combined <- data.frame(matrix(NA, nrow = length(all_species), ncol = 1))
rownames(sig_tab_combined) <- all_species
colnames(sig_tab_combined) <- 'combined'

for (sp in all_species) {
    
    tab_subset <- combined_strain_div[which(combined_strain_div$species == sp), ]
    
    within_vec <- tab_subset[which(tab_subset$Comparison == 'Within'), 'mean_div']
    between_vec <- tab_subset[which(tab_subset$Comparison == 'Between'), 'mean_div']
    
    if (length(within_vec) < 5 || length(between_vec) < 5) {
      sig_tab_combined[sp, 1] <- 'Insufficient data'
    } else {
      wilcox_out <- wilcox.test(within_vec, between_vec, exact = FALSE)
      
      wilcox_p <- wilcox_out$p.value
      
      if (is.na(wilcox_p)) {
        wilcox_p <- 1 
      }
      
      if (wilcox_p >= 0.05) {
        sig_tab_combined[sp, 1] <- 'Non-significant'
      } else if (mean(within_vec) > mean(between_vec)) {
        sig_tab_combined[sp, 1] <- 'Higher within'
      } else if (mean(within_vec) < mean(between_vec)) {
        sig_tab_combined[sp, 1] <- 'Higher between'
      } else {
        stop('Error!') 
      }
    }
}

Heatmap(matrix = as.matrix(sig_tab_combined),
        col = c('blue', 'grey', 'grey95'),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(sig_tab_combined), '*', sep = '')))
