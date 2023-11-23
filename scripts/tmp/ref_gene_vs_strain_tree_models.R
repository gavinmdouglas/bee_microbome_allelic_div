rm(list = ls(all.names = TRUE))

library(ggplot2)

tree_results <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_vs_species_dist_and_D.tsv.gz',
                           header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Remove genes for which neither tree distances nor D values could be computed.
tree_results <- tree_results[-which(is.na(tree_results$tree_dist) & is.na(tree_results$D)), ]

# Only consider genes with at least 5 tips, as otherwise the tree distance is always just 1.
tree_results <- tree_results[which(tree_results$num_tips > 4), ]


tree_results_core <- tree_results[which(tree_results$prevalence == 1), ]
tree_results_acc <- tree_results[which(tree_results$prevalence < 1), ]

ggplot(data = tree_results_core, aes(x = tree_dist, y = species)) + geom_boxplot()
ggplot(data = tree_results_acc, aes(x = tree_dist, y = species)) + geom_boxplot()

ggplot(data = tree_results, aes(x = D, y = species)) + geom_boxplot()

# Run orderedNorm transformation on D values to make them normally distributed
# (although they should not be interpreted as standard scores/ranks).
tree_results$D_orderedNorm <- bestNormalize::orderNorm(tree_results$D)$x.t

tree_results$num_tips_factor <- factor(as.character(tree_results$num_tips))

D_model_core <- glmmTMB::glmmTMB(formula = D ~ species + num_tips,
                            data = tree_results_core,
                            family = 'gaussian',
                            control = glmmTMBControl(optimizer = nlminb,
                                                     parallel = 10,
                                                     profile = TRUE,
                                                     optCtrl = list(iter.max = 100000,
                                                                    eval.max = 100000)))
tree_results$tree_results_binary <- 0
tree_results$tree_results_binary[which(tree_results$tree_dist > 0.5)] <- 1



tree_results$nodes_to_tips <- tree_results$num_internal_nodes / tree_results$num_tips

D_model <- lm(formula = D_orderedNorm ~ species + num_tips, data = tree_results)
dist_model <- lm(formula = tree_dist ~ species + num_tips, data = tree_results)


ggplot(data = tree_results, aes(x = num_tips, y = species)) + geom_boxplot()

tmp <- bestNormalize::bestNormalize(tree_results$tree_dist / tree_results$num_tips)

tree_results$dist_by_num_tips_binary <- 0
tree_results$dist_by_num_tips_binary[which(tmp$x.t > 0.5)] <- 1

dist_model <- glmmTMB::glmmTMB(formula = dist_by_num_tips_binary ~ species,
                            data = tree_results,
                            family = 'binomial',
                            control = glmmTMBControl(optimizer = nlminb,
                                                     parallel = 10,
                                                     profile = TRUE,
                                                     optCtrl = list(iter.max = 100000,
                                                                    eval.max = 100000)))

tree_results_atleast5 <- tree_results[which(tree_results$num_tips >= 5), ]
ggplot(data = tree_results_atleast5, aes(x = tree_dist, y = species)) + geom_boxplot()

boxplot(data = tree_results, nodes_to_tips ~ tree_results_binary)

tree_results_Ntips <- tree_results[which(tree_results$num_tips == 6), ]
ggplot(data = tree_results_Ntips, aes(x = tree_dist, y = species)) + geom_boxplot()

tree_results$below_ten_tips <- 'Yes'
tree_results[which(tree_results$num_tips >= 10), 'below_ten_tips'] <- 'No'
ggplot(data = tree_results, aes(x = tree_dist, y = species, fill = below_ten_tips)) + geom_boxplot()

summary(lm(formula = tree_dist ~ species, data = tree_results_Ntips))

tree_results_Ntips$tree_dist_norm <- bestNormalize::bestNormalize(tree_results_Ntips$tree_dist)$x.t

summary(lm(formula = tree_dist ~ species, data = tree_results_Ntips))


dist_model <- glmmTMB::glmmTMB(formula = dist_by_num_tips_binary ~ species,
                               data = tree_results,
                               family = 'binomial',
                               control = glmmTMBControl(optimizer = nlminb,
                                                        parallel = 10,
                                                        profile = TRUE,
                                                        optCtrl = list(iter.max = 100000,
                                                                       eval.max = 100000)))
