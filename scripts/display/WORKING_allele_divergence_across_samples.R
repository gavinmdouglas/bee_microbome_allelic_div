rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

within_div <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/within_sample_allele_divergence.tsv.gz',
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)

between_div <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/between_sample_allele_divergence.tsv.gz',
                          header = TRUE, sep = '\t', stringsAsFactors = FALSE)

within_div_clean <- within_div
within_div_clean <- within_div_clean[, -1]
colnames(within_div_clean)[1] <- 'mean_div'
#within_div_clean <- aggregate(x = mean_div ~ dataset + gene + species, data = within_div_clean, FUN = mean)
within_div_clean$Comparison <- 'Within'

between_div_clean <- between_div
between_div_clean <- between_div_clean[, -c(1, 2)]
colnames(between_div_clean)[1] <- 'mean_div'
#between_div_clean <- aggregate(x = mean_div ~ dataset + gene + species, data = between_div_clean, FUN = mean)
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

sig_tab_higher_within <- data.frame(matrix(0, nrow = length(all_species), ncol = length(datasets)))
rownames(sig_tab_higher_within) <- all_species
colnames(sig_tab_higher_within) <- datasets

sig_tab_higher_between <- sig_tab_higher_within
total_tab <- sig_tab_higher_within

for (sp in all_species) {
 
   for (d in datasets) {

      tab_subset <- combined_strain_div[which(combined_strain_div$species == sp & combined_strain_div$dataset == d), ]
      
      for (gene in unique(tab_subset$gene)) {
      
        tab_subset_gene <- tab_subset[which(tab_subset$gene == gene), , drop = FALSE]
        
        within_vec <- tab_subset_gene[which(tab_subset_gene$Comparison == 'Within'), 'mean_div']
        between_vec <- tab_subset_gene[which(tab_subset_gene$Comparison == 'Between'), 'mean_div']
        
        if (length(within_vec) < 10 || length(between_vec) < 10) {
          next
        } else {
          
          total_tab[sp, d] <- total_tab[sp, d] + 1
          
          wilcox_out <- wilcox.test(within_vec, between_vec, exact = FALSE)
          
          wilcox_p <- wilcox_out$p.value
          
          if (is.na(wilcox_p)) {
            wilcox_p <- 1 
          }
          
          if (wilcox_p < 0.05) {
            
            if (mean(within_vec) < mean(between_vec)) {
              sig_tab_higher_between[sp, d] <- sig_tab_higher_between[sp, d] + 1
            } else if(mean(within_vec) > mean(between_vec)) {
              sig_tab_higher_within[sp, d] <- sig_tab_higher_within[sp, d] + 1
            }
          }
        }
      }
   }
}

# Only consider species/datasets with > 50 testable genes.
sig_tab_higher_within[total_tab < 50] <- NA
sig_tab_higher_between[total_tab < 50] <- NA
total_tab[total_tab < 50] <- NA

percent_higher_within <- (sig_tab_higher_within / total_tab) * 100
percent_higher_within <- apply(percent_higher_within, 2, function(x) { round(x, 2)})

percent_higher_between <- (sig_tab_higher_between / total_tab) * 100
percent_higher_between <- apply(percent_higher_between, 2, function(x) { round(x, 2)})

Heatmap(matrix = as.matrix(percent_higher_within),
        col = circlize::colorRamp2(c(0, 30), c('cornflowerblue', 'red')),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(percent_higher_within), '*', sep = '')),

        cluster_rows = FALSE,
        cluster_columns = FALSE,

        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(percent_higher_within[i, j] > 0))
            grid.text(percent_higher_within[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })

Heatmap(matrix = as.matrix(percent_higher_between),
        col = circlize::colorRamp2(c(0, 30), c('cornflowerblue', 'red')),
        heatmap_legend_param = list(title = ''),
        row_labels = gt_render(paste('*', rownames(percent_higher_between), '*', sep = '')),
        
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(percent_higher_between[i, j] > 0))
            grid.text(percent_higher_between[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })
