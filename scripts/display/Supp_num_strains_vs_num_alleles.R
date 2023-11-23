rm(list = ls(all.names = TRUE))

library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)

source('~/scripts/bee_microbome_allelic_div/scripts/functions.R')

breakdown <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/num_strain_vs_num_allele_breakdown.tsv.gz',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
breakdown <- breakdown[-which(is.na(breakdown$cor_present_p)), ]
breakdown <- breakdown[which(breakdown$sample_prev >= 20), ]
breakdown <- breakdown[which(breakdown$num_alleles >= 5), ]
breakdown <- breakdown[which(breakdown$prop_samples_w_one_strain < 0.5), ]
breakdown <- breakdown[which(breakdown$mean_number_of_strains >= 2), ]

species_tally <- table(breakdown$Species)
rare_species <- names(species_tally)[which(species_tally < 50)]
breakdown <- breakdown[which(! breakdown$Species %in% rare_species), ]

breakdown$Species <- cleanup_species_names(breakdown$Species)

# First get boxplots of tau per species.
tau_boxplots <- ggplot(
  data = breakdown, aes(x = tau_present, y = Species)) +
  geom_violin(col='grey80') +
  geom_boxplot(fill = 'grey50', outlier.shape = NA, alpha = 0.2) +
  theme_bw() +
  theme(axis.text.y = element_text(face = 'italic')) +
  scale_y_discrete(limits=rev) +
  xlab(expression("Kendall's "*tau*" per accessory gene"))

# Then get breakdown of % and numbers of genes significant based on multinomial test.
unique_species <- sort(unique(breakdown$Species))
multinomial_sig_breakdown <- data.frame(Species = unique_species,
                                        sig = NA,
                                        nonsig = NA)
rownames(multinomial_sig_breakdown) <- unique_species
for (sp in rownames(multinomial_sig_breakdown)) {
  multinomial_sig_breakdown[sp, 'sig'] <- length(which(breakdown[which(breakdown$Species == sp), 'multinomial_p'] < 0.05))
  multinomial_sig_breakdown[sp, 'nonsig'] <- length(which(breakdown[which(breakdown$Species == sp), 'multinomial_p'] >= 0.05))
}

multinomial_sig_breakdown <- multinomial_sig_breakdown[, c(2, 3)]

multinomial_sig_breakdown_percent <- (multinomial_sig_breakdown / rowSums(multinomial_sig_breakdown)) * 100

multinomial_sig_breakdown_char <- multinomial_sig_breakdown_percent

multinomial_sig_breakdown_char$sig <- paste(as.character(multinomial_sig_breakdown$sig),
                                            ' (',
                                            gsub(' ', '', format(round(multinomial_sig_breakdown_percent$sig, digits=1), nsmall = 1)),
                                            '%)', sep = '')

multinomial_sig_breakdown_char$nonsig <- paste(as.character(multinomial_sig_breakdown$nonsig),
                                            ' (',
                                            gsub(' ', '', format(round(multinomial_sig_breakdown_percent$nonsig, digits=1), nsmall = 1)),
                                            '%)', sep = '')

multinomial_sig_heatmap <- Heatmap(
  matrix = as.matrix(multinomial_sig_breakdown_percent),
  name = "Percent\nsig.",
  row_names_side = "left",
  row_labels = gt_render(paste("*", rownames(multinomial_sig_breakdown_percent), "*", sep = "")),
  column_labels = c('Significant', 'Non-significant'),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  col = colorRamp2(c(0, 100), c("white", "red")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(multinomial_sig_breakdown_char[i, j], x, y, gp = gpar(fontsize = 10, fontface="bold"))
  }
)


combined_plot <- plot_grid(
  tau_boxplots,
  plot_grid(grid.grabExpr(draw(multinomial_sig_heatmap))),
  labels = c('a', 'b'),
  nrow = 1
)

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_num_strains_vs_num_alleles_overview.pdf',
       plot = combined_plot,
       device = 'pdf',
       dpi = 600,
       width = 12,
       height = 8)
