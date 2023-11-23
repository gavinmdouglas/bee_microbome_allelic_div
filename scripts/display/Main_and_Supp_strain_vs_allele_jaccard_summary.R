rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(glmmTMB)

model_outputs <- readRDS(file = '/data1/gdouglas/projects/honey_bee/large_files_to_backup/statistics/strain_vs_allele_jaccard_models.rds')

AIC_table <- data.frame(matrix(NA, nrow = length(model_outputs), ncol = 2))
colnames(AIC_table) <- c('predictors', 'AIC')

for (i in 1:length(model_outputs)) {
  AIC_table[i, 'predictors'] <- names(model_outputs)[i]
  AIC_table[i, 'AIC'] <- summary(model_outputs[[i]])$AICtab[1]
}

AIC_table$predictors <- factor(AIC_table$predictors, levels = rev(unique(AIC_table$predictors)))

AIC_table$min_by_AIC <- AIC_table$AIC - min(AIC_table$AIC)
AIC_table$min_by_AIC_norm <- AIC_table$min_by_AIC / max(AIC_table$min_by_AIC)
AIC_table$AIC_round <- format(round(AIC_table$AIC, digits=0), nsmall = 0)

AIC_table$hjust_setting <- -0.05
AIC_table$text_colour <- "black"
AIC_table[which(AIC_table$min_by_AIC_norm >= 0.8), "hjust_setting"] <- 1.075
AIC_table[which(AIC_table$min_by_AIC_norm >= 0.8), "text_colour"] <- "white"
norm_AIC_barplot <- ggplot(data = AIC_table, aes(x = min_by_AIC_norm, y = predictors)) +
                            geom_bar(stat = "identity", fill = "grey20") +
                            geom_text(aes(label=AIC_round), colour = AIC_table$text_colour, hjust = AIC_table$hjust_setting, size = 3) +
                            xlab("Normalized AIC: (AIC - min[AIC]) / max(AIC - min(AIC))") +
                            ylab("") +
                            theme_bw() +
                            theme(axis.text.y=element_blank()) +
                            coord_cartesian(expand = FALSE) +
                            theme(plot.margin = unit(c(2, 0, 2, 5), 'mm'))

dummy <- data.frame(matrix('Excluded', nrow = nrow(AIC_table), ncol = 8))
colnames(dummy) <- c('species', 'COG_category', 'dataset', 'num_samples', 'num_alleles', 'num_strainfacts_strains', 'allelic_jaccard', 'mean_strainfacts_jaccard')
rownames(dummy) <- AIC_table$predictors
for (i in 1:nrow(AIC_table)) {
  for (var_part in colnames(dummy)) {
     if (length(grep(var_part, rownames(dummy)[i])) == 1) {
       dummy[i, var_part] <- 'Included'
     }
  }
}

model_info_heatmap <- Heatmap(as.matrix(dummy),
                              col = c('grey95', 'grey20'),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              column_split = c('Fixed effects', 'Fixed effects',
                                               'Random effects',
                                               'Random effects', 'Random effects',
                                               'Random effects', 'Random effects', 'Random effects'),
                              column_gap = unit(5, "mm"),
                              show_row_names = FALSE,
                              heatmap_legend_param = list(title = 'Model variables'),
                              column_names_rot = 45,
                              column_labels = c('Species', 'COG category', 'Dataset', 'No. samples',
                                                'No. alleles', 'No. strains', 'Jaccard (alleles)',
                                                'Jaccard (strain)'))

norm_AIC_barplot_w_whitespace <- plot_grid(NULL, norm_AIC_barplot, NULL, nrow = 3, rel_heights = c(3, 85, 7.5))

AIC_barplot_panel <- plot_grid(grid.grabExpr(draw(model_info_heatmap,
                                                  heatmap_legend_side = 'left',
                                                  padding = unit(c(2, 2, 2, 0), 'mm'))),
                               norm_AIC_barplot_w_whitespace,
                               NULL,
                               ncol = 3,
                               rel_widths = c(1, 2, 0.08))


# Significant coefficients of best fit model.
best_model <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/model_outputs/strain_vs_allele_jaccard_model_best.rds')
best_fit_summary <- summary(best_model)

# Then plot coefficients for all species, in same order as earlier plots.
species_coefficients <- data.frame(best_fit_summary$coefficients$cond[grep('^species', rownames(best_fit_summary$coefficients$cond)), ])
species_coefficients$variable <- gsub('^species', '', rownames(species_coefficients))
species_coefficients$variable <- gsub('_', ' ', species_coefficients$variable)
species_coefficients$variable <- gsub(' sp$', ' sp.', species_coefficients$variable)
species_coefficients$variable <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', species_coefficients$variable)

species_coefficients$P_category <- 'P >= 0.05'
species_coefficients$P_category[which(species_coefficients$Pr...z.. < 0.05)] <- 'P < 0.05'

species_coefficients <- species_coefficients[order(species_coefficients$Estimate), ]

species_coefficients$variable <- factor(species_coefficients$variable,
                                        levels = species_coefficients$variable)

species_coefficients_barplot <- ggplot(data = species_coefficients, aes(x = Estimate, y = variable, fill = P_category)) +
                                    geom_bar(stat="identity") +
                                    scale_fill_manual(values = c('firebrick', 'grey85')) +
                                    theme_bw() +
                                    ylab('Species') +
                                    xlab('Intercept coefficient') +
                                    theme(axis.text.y = element_text(face = 'italic'),
                                          plot.title = element_text(hjust = 0.5)) +
                                    geom_vline(xintercept = 0, linetype="dotted", 
                                               color = "black") +
                                    geom_errorbar(aes(xmin = Estimate - Std..Error,
                                                      xmax = Estimate + Std..Error),
                                                  width = 0.2, color = "black") +
                                    guides(fill = guide_legend('Coefficient\ntest result')) +
                                    ggtitle('Normalized minimum Jaccard distances\nbetween allele and strain profiles by species')
                                  

# Output main plot.
ggsave(plot = species_coefficients_barplot,
       filename = "/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Main_strain_vs_allele_min_jaccard_by_species.pdf",
       device = "pdf", width = 7, height = 7, units = "in", dpi = 400)

# And supplementary plot:
ggsave(plot = AIC_barplot_panel,
       filename = "/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_strain_vs_allele_model_AICs.pdf",
       device = "pdf", width = 12, height = 8.5, units = "in", dpi = 400)
