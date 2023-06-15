rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

min_jaccard_summary <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv.gz',
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)
min_jaccard_summary$species <- gsub('_', ' ', min_jaccard_summary$species)
min_jaccard_summary$species <- gsub(' sp$', ' sp.', min_jaccard_summary$species)
min_jaccard_summary$species <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', min_jaccard_summary$species)

# Remove all genes found in fewer than 20 samples.
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_samples >= 20), ]

# Ignore Commensalibacter_sp due to low sample size
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$species != 'Commensalibacter sp.'), ]

mean_by_species <- aggregate(x = mean_min ~ species, FUN = mean, data = min_jaccard_summary)
species_order <- mean_by_species$species[order(mean_by_species$mean_min, decreasing = FALSE)]

min_jaccard_summary$species <- factor(min_jaccard_summary$species, levels = species_order)

min_jaccard_summary$permuted_mean_min_P_sig <- 'P >= 0.05'
min_jaccard_summary$permuted_mean_min_P_sig[which(min_jaccard_summary$permuted_mean_min_P < 0.05)] <- 'P < 0.05'

min_jaccard_summary$log2_obs_vs_permuted_mean_min <- log2((min_jaccard_summary$mean_min + 0.1) / (min_jaccard_summary$permuted_mean_min + 0.1))

# Put log2 ratio and raw value into same panel ,as there is a lot of repeated labels for both.
jaccard_faceted_tmp1 <- min_jaccard_summary[, c('species', 'mean_min', 'permuted_mean_min_P_sig')]
colnames(jaccard_faceted_tmp1)[2] <- 'value'
jaccard_faceted_tmp1$variable_type <- "'Mean minimum Jaccard distance\nacross alleles (per gene)'"

jaccard_faceted_tmp2 <- min_jaccard_summary[, c('species', 'log2_obs_vs_permuted_mean_min', 'permuted_mean_min_P_sig')]
colnames(jaccard_faceted_tmp2)[2] <- 'value'
jaccard_faceted_tmp2$variable_type <- "'log'[2]*'([Obs. Jaccard] / [Permuted Jaccard])'"

jaccard_faceted <- rbind(jaccard_faceted_tmp1, jaccard_faceted_tmp2)

jaccard_faceted$variable_type <- factor(jaccard_faceted$variable_type,
                                        levels = unique(jaccard_faceted$variable_type))

combined_mean_min_distributions <- ggplot(data = jaccard_faceted, aes(y = species, x = value)) +
                                      geom_quasirandom(aes(col = permuted_mean_min_P_sig), size = 0.5, ) +
                                      scale_colour_manual(values = c('firebrick', 'grey85')) +
                                      geom_boxplot(outlier.shape = NA, fill = 'white', alpha = 0.1) +
                                      theme_bw() +
                                      xlab('') +
                                      ylab('Species') +
                                      theme() +
                                      guides(col = guide_legend('Permutation\ntest result')) +
                                      facet_wrap(. ~ variable_type,
                                                 labeller = label_parsed,
                                                 scales = 'free_x',
                                                 strip.position = 'bottom') +
                                      theme(axis.text.y = element_text(face = 'italic'),
                                            strip.background = element_blank(),
                                            strip.placement = "outside",
                                            strip.text = element_text(size = 11))

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

dummy <- data.frame(matrix('Excluded', nrow = nrow(AIC_table), ncol = 7))
colnames(dummy) <- c('species', 'COG_category', 'num_samples', 'num_alleles', 'num_strains', 'allelic_jaccard', 'strain_jaccard')
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
                                               'Random effects', 'Random effects',
                                               'Random effects', 'Random effects', 'Random effects'),
                              column_gap = unit(5, "mm"),
                              show_row_names = FALSE,
                              heatmap_legend_param = list(title = 'Model variables'),
                              column_names_rot = 45,
                              column_labels = c('Species', 'COG category', 'No. samples',
                                                'No. alleles', 'No. strains', 'Jaccard (alleles)',
                                                'Jaccard (strain)'))

norm_AIC_barplot_w_whitespace <- plot_grid(NULL, norm_AIC_barplot, NULL, nrow = 3, rel_heights = c(3, 85, 7.5))

AIC_barplot_panel <- plot_grid(grid.grabExpr(draw(model_info_heatmap,
                                                  heatmap_legend_side = 'left',
                                                  padding = unit(c(2, 2, 2, 0), 'mm'))),
                               norm_AIC_barplot_w_whitespace,
                               NULL,
                               ncol = 3,
                               rel_widths = c(1, 2, 0.37))


# Significant coefficients of species + random effect model (best fit).
species_full_summary <- summary(model_outputs$`species + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_sample_strain_jaccard)`)

# First confirm that the only significant coefficients are species:
species_full_summary$coefficients$cond[which(species_full_summary$coefficients$cond[, 4] < 0.05), ]


# Then plot coefficients for all species, in same order as earlier plots.
species_coefficients <- data.frame(species_full_summary$coefficients$cond[grep('^species', rownames(species_full_summary$coefficients$cond)), ])
species_coefficients$variable <- gsub('^species', '', rownames(species_coefficients))
species_coefficients$variable <- gsub('_', ' ', species_coefficients$variable)
species_coefficients$variable <- gsub(' sp$', ' sp.', species_coefficients$variable)
species_coefficients$variable <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', species_coefficients$variable)

species_coefficients$P_category <- 'P >= 0.05'
species_coefficients$P_category[which(species_coefficients$Pr...z.. < 0.05)] <- 'P < 0.05'

species_coefficients$variable <- factor(species_coefficients$variable,
                                        levels = species_order)

species_coefficients_barplot <- ggplot(data = species_coefficients, aes(x = Estimate, y = variable, fill = P_category)) +
                                    geom_bar(stat="identity") +
                                    scale_fill_manual(values = c('firebrick', 'grey85')) +
                                    theme_bw() +
                                    ylab('Species') +
                                    xlab('Intercept coefficient') +
                                    theme(axis.text.y = element_text(face = 'italic')) +
                                    geom_vline(xintercept = 0, linetype="dotted", 
                                               color = "black") +
                                    geom_errorbar(aes(xmin = Estimate - Std..Error,
                                                      xmax = Estimate + Std..Error),
                                                  width = 0.2, color = "black") +
                                    guides(fill = guide_legend('Coefficient\ntest result'))



# Output supplementary plots.
ggsave(plot = combined_mean_min_distributions,
       filename = "/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_strain_vs_allele_distributions.pdf",
       device = "pdf", width = 9, height = 6, units = "in", dpi = 400)


strain_vs_allele_jaccard_models <- plot_grid(AIC_barplot_panel,
                                            species_coefficients_barplot,
                                            nrow = 2,
                                            rel_heights = c(2, 1),
                                            labels = c('a', 'b'))

ggsave(plot = strain_vs_allele_jaccard_models,
       filename = "/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_strain_vs_allele_model_comparisons.pdf",
       device = "pdf", width = 12, height = 13.1, units = "in", dpi = 400)
