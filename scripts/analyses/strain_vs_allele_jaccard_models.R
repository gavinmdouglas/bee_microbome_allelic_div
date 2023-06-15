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

split_multi_category_rows <- function(in_df, category_col, delimiter = ",", num_cores = 1) {
  
  multi_category_row_i <- grep(delimiter, in_df[, category_col])
  
  in_df_nonmulti <- in_df[-multi_category_row_i, , drop = FALSE]
  in_df_multi <- in_df[multi_category_row_i, , drop = FALSE]
  
  in_df_multi_merged_raw <- parallel::mclapply(1:nrow(in_df_multi),
                                               function(i) {
                                                 row_subset <- in_df_multi[i, , drop = FALSE]
                                                 category_split <- base::strsplit(x = row_subset[, category_col], split = ",")[[1]]
                                                 split_row <- row_subset[rep(1, length(category_split)), , drop = FALSE]
                                                 split_row[, category_col] <- category_split
                                                 return(split_row)
                                               },
                                               mc.cores = num_cores)
  
  in_df_multi_merged <- do.call(rbind, in_df_multi_merged_raw)
  
  return(rbind(in_df_nonmulti, in_df_multi_merged))
}

min_jaccard_summary <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_allele_min_jaccard_breakdown.tsv.gz',
                                  header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)

# Remove all genes found in fewer than 20 samples.
min_jaccard_summary <- min_jaccard_summary[which(min_jaccard_summary$num_samples >= 20), ]

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

min_jaccard_ORQ$COG_category <- NA
intersecting_genes <- intersect(rownames(min_jaccard_ORQ), rownames(COG_category_info))
min_jaccard_ORQ[intersecting_genes, 'COG_category'] <- COG_category_info[intersecting_genes, 'COG_category']
min_jaccard_ORQ[which(min_jaccard_ORQ$COG_category == '-'), 'COG_category'] <- NA

# Ignore rows that could not be annotated with COG categories.
min_jaccard_ORQ <- min_jaccard_ORQ[which(! is.na(min_jaccard_ORQ$COG_category)), ]

# Ignore Commensalibacter_sp due to low sample size.
min_jaccard_ORQ <- min_jaccard_ORQ[which(min_jaccard_ORQ$species != 'Commensalibacter_sp'), ]

min_jaccard_ORQ <- split_multi_category_rows(in_df = min_jaccard_ORQ,
                                             category_col = 'COG_category',
                                             delimiter = ',',
                                             num_cores = 10)

# Remove rare COG categories (< 100 instances).
COG_category_tallies <- table(min_jaccard_ORQ$COG_category)

min_jaccard_ORQ <- min_jaccard_ORQ[which(! min_jaccard_ORQ$COG_category %in% names(COG_category_tallies)[which(COG_category_tallies < 100)]), ]

# Use COG Category E as reference as it is abundant and centred around 0.
COG_category_order <- sort(unique(min_jaccard_ORQ$COG_category))
COG_category_order <- c('E', COG_category_order[-which(COG_category_order == 'E')])
min_jaccard_ORQ$COG_category <- factor(min_jaccard_ORQ$COG_category, levels = COG_category_order)

# Chose 'Lactobacillus_melliventris' as the species reference, as this species tended to have intermediate mean min values and had the highest sample size of all core species.
species_order <- sort(unique(min_jaccard_ORQ$species))
species_order <- c('Lactobacillus_melliventris', species_order[-which(species_order == 'Lactobacillus_melliventris')])
min_jaccard_ORQ$species <- factor(min_jaccard_ORQ$species, levels = species_order)

raw_predictor_strings <- c('species',
                           'species + (1 | num_samples)',
                           'species + (1 | num_samples) + (1 | num_alleles)',
                           'species + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains)',
                           'species + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard)',
                           'species + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_sample_strain_jaccard)',

                           'COG_category',
                           'COG_category + (1 | num_samples)',
                           'COG_category + (1 | num_samples) + (1 | num_alleles)',
                           'COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains)',
                           'COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard)',
                           'COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_sample_strain_jaccard)',
                           
                           'species + COG_category',
                           'species + COG_category + (1 | num_samples)',
                           'species + COG_category + (1 | num_samples) + (1 | num_alleles)',
                           'species + COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains)',
                           'species + COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard)',
                           'species + COG_category + (1 | num_samples) + (1 | num_alleles) + (1 | num_strains) + (1 | mean_sample_allelic_jaccard) + (1 | mean_sample_strain_jaccard)'
                           )

model_outputs <- list()

for (raw_predictors in raw_predictor_strings) {
 
  model_formula <- as.formula(paste('mean_min ~ ', raw_predictors, sep = ''))
   
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
