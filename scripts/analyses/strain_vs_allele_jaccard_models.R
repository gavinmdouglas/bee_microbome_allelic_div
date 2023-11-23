rm(list = ls(all.names = TRUE))

# Fit linear mixed-effects model of transformed mean min. Jaccard distance between alleles and strains
# vs. fixed effects:
#   species,
#   COG-categories,
#   and interaction
# and random effects: 
#   species' strains mean pairwise Jaccard (based on samples),
#   gene's alleles mean pairwise Jaccard (based on samples),
#   Num of samples
#   Num of alleles
#   Num of strains
#   (and interactions as appropriate)

library(bestNormalize)
library(glmmTMB)

source('~/scripts/bee_microbome_allelic_div/scripts/functions.R')

min_jaccard_summary <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv.gz',
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Remove rare features.
min_num_samples <- 10
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_samples >= min_num_samples), ]

min_num_features <- 2
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_alleles >= min_num_features), ]
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_strainfacts_strains >= min_num_features), ]
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_straingst_strains >= min_num_features), ]

## Used bestNormalize to compare transformation approaches for all quantitative variables.
## Which turned out to be ordered quantile normalizing transformation by far in each case.
# mean_min_bn <- bestNormalize::bestNormalize(min_jaccard_summary$mean_min)
# num_samples_bn <- bestNormalize::bestNormalize(min_jaccard_summary$num_samples)
# num_alleles_bn <- bestNormalize::bestNormalize(min_jaccard_summary$num_alleles)
# num_strains_bn <- bestNormalize::bestNormalize(min_jaccard_summary$num_strains)
# mean_sample_allelic_jaccard_bn <- bestNormalize::bestNormalize(min_jaccard_summary$mean_sample_allelic_jaccard)
# mean_sample_strain_jaccard_bn <- bestNormalize::bestNormalize(min_jaccard_summary$mean_sample_strain_jaccard)

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

# Add in COG category information.
COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

min_jaccard_ORQ <- min_jaccard_ORQ[which(min_jaccard_ORQ$gene %in% rownames(COG_category_info)), ]
min_jaccard_ORQ$COG_category <- COG_category_info[min_jaccard_ORQ$gene, 'COG_category']

min_jaccard_ORQ <- split_multi_category_rows(in_df = min_jaccard_ORQ,
                                             category_col = 'COG_category',
                                             delimiter = ',',
                                             num_cores = 10)

# Remove rare COG categories (< 100 instances).
COG_category_tallies <- table(min_jaccard_ORQ$COG_category)
min_jaccard_ORQ <- min_jaccard_ORQ[which(! min_jaccard_ORQ$COG_category %in% names(COG_category_tallies)[which(COG_category_tallies < 100)]), ]

# Remove unannotated genes:
min_jaccard_ORQ <- min_jaccard_ORQ[which(min_jaccard_ORQ$COG_category != '-'), ]

# Use COG Category E as reference as it is abundant and centred around 0.
COG_category_order <- sort(unique(min_jaccard_ORQ$COG_category))
COG_category_order <- c('E', COG_category_order[-which(COG_category_order == 'E')])
min_jaccard_ORQ$COG_category <- factor(min_jaccard_ORQ$COG_category, levels = COG_category_order)

# Remove rare species.
sufficiently_frequent_species <- names(table(min_jaccard_ORQ$species))[which(table(min_jaccard_ORQ$species) >= 100)]
min_jaccard_ORQ <- min_jaccard_ORQ[which(min_jaccard_ORQ$species %in% sufficiently_frequent_species), ]

# Chose 'Lactobacillus_melliventris' as the species reference, as this species tended to have intermediate mean min values and had the highest sample size of all core species.
species_order <- sort(unique(min_jaccard_ORQ$species))
species_order <- c('Lactobacillus_melliventris', species_order[-which(species_order == 'Lactobacillus_melliventris')])
min_jaccard_ORQ$species <- factor(min_jaccard_ORQ$species, levels = species_order)

raw_predictor_strings <- c('species',
                           'species + (1 | dataset)',
                           'species + (1 | dataset) + (1 | num_samples)',
                           'species + (1 | dataset) + (1 | num_samples) + (1 | num_alleles)',
                           'species + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains)',
                           'species + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard)',
                           'species + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_strainfacts_jaccard)',

                           'COG_category',
                           'COG_category + (1 | dataset)',
                           'COG_category + (1 | dataset) + (1 | num_samples)',
                           'COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles)',
                           'COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains)',
                           'COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard)',
                           'COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_strainfacts_jaccard)',
                           
                           'species + COG_category',
                           'species + COG_category + (1 | dataset)',
                           'species + COG_category + (1 | dataset) + (1 | num_samples)',
                           'species + COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles)',
                           'species + COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains)',
                           'species + COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard)',
                           'species + COG_category + (1 | dataset) + (1 | num_samples) + (1 | num_alleles) + (1 | num_strainfacts_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_strainfacts_jaccard)'
                           )

model_outputs <- list()

for (raw_predictors in raw_predictor_strings) {
 
  print(raw_predictors)
  
  model_formula <- as.formula(paste('strainfacts_mean_min ~ ', raw_predictors, sep = ''))
   
  model_outputs[[raw_predictors]] <- glmmTMB::glmmTMB(formula = model_formula,
                                                      data = min_jaccard_ORQ,
                                                      family = 'gaussian',
                                                      control = glmmTMBControl(optimizer = nlminb,
                                                                               parallel = 10,
                                                                               profile = TRUE,
                                                                               optCtrl = list(iter.max = 1000,
                                                                                              eval.max = 1000)))
  
}

saveRDS(object = model_outputs,
        file = '/data1/gdouglas/projects/honey_bee/large_files_to_backup/statistics/strain_vs_allele_jaccard_models.rds')
