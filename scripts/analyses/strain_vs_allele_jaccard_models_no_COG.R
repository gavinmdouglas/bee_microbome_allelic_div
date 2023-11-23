rm(list = ls(all.names = TRUE))

# Re-run best fitting model (i.e., all random effects, but not COG categories as fixed effect).
# Do so because it is no longer necessary to duplicate rows for genes with multiple COG category annotations.
# Can also include *all* genes, regardless of whether they are COG-annotated.
# Also, can run GLMM fit with more iterations (100,000 rather than 1,000).

library(bestNormalize)
library(glmmTMB)

min_jaccard_summary <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv.gz',
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Remove rare features.
min_num_samples <- 10
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_samples >= min_num_samples), ]

min_num_features <- 2
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_alleles >= min_num_features), ]
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_strainfacts_strains >= min_num_features), ]
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_straingst_strains >= min_num_features), ]

# Transformed to ~normal with ordered quantile normalizing transformation.
ORQ_info <- list()
min_jaccard_ORQ <- min_jaccard_summary

for (variable in colnames(min_jaccard_summary)) {
  if (inherits(x = min_jaccard_summary[, variable], 'character')) {
    next
  }
  ORQ_info[[variable]] <- bestNormalize::orderNorm(min_jaccard_summary[, variable])
  min_jaccard_ORQ[, variable] <- ORQ_info[[variable]]$x.t
}

# Remove rare species.
sufficiently_frequent_species <- names(table(min_jaccard_ORQ$species))[which(table(min_jaccard_ORQ$species) >= 100)]
min_jaccard_ORQ <- min_jaccard_ORQ[which(min_jaccard_ORQ$species %in% sufficiently_frequent_species), ]

# Chose 'Lactobacillus_melliventris' as the species reference, as this species tended to have intermediate mean min values and had the highest sample size of all core species.
species_order <- sort(unique(min_jaccard_ORQ$species))
species_order <- c('Lactobacillus_melliventris', species_order[-which(species_order == 'Lactobacillus_melliventris')])
min_jaccard_ORQ$species <- factor(min_jaccard_ORQ$species, levels = species_order)

model_formula <- as.formula('strainfacts_mean_min ~ species + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_strainfacts_jaccard)')
   
model_out <- glmmTMB::glmmTMB(formula = model_formula,
                              data = min_jaccard_ORQ,
                              family = 'gaussian',
                              control = glmmTMBControl(optimizer = nlminb,
                                                       parallel = 10,
                                                       profile = TRUE,
                                                       optCtrl = list(iter.max = 100000,
                                                                      eval.max = 100000)))

saveRDS(object = model_out,
        file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/model_outputs/strain_vs_allele_jaccard_model_best.rds')
