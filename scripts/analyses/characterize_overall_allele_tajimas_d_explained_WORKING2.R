rm(list = ls(all.names = TRUE))

library(bestNormalize)
library(glmmTMB)
library(performance)

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

allele_popgen_by_species <- readRDS(file = '/home/gdouglas/tmp/allele_popgen_by_species.rds')


# Add in COG category information.
COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

allele_popgen_by_species$COG_category <- '-'
allele_popgen_by_species$COG <- '-'
intersecting_genes_i <- which(allele_popgen_by_species$gene %in% rownames(COG_category_info))
allele_popgen_by_species[intersecting_genes_i, 'COG_category'] <- COG_category_info[allele_popgen_by_species$gene[intersecting_genes_i], 'COG_category']
allele_popgen_by_species[intersecting_genes_i, 'COG'] <- COG_category_info[allele_popgen_by_species$gene[intersecting_genes_i], 'all_COG']
allele_popgen_by_species[which(allele_popgen_by_species$COG == ''), 'COG'] <- '-'

allele_popgen_by_species <- allele_popgen_by_species[which(allele_popgen_by_species$mean_num_reads >= 50), ]
allele_popgen_by_species <- allele_popgen_by_species[which(allele_popgen_by_species$num_seqs >= 3), ]
allele_popgen_by_species <- allele_popgen_by_species[which(! is.na(allele_popgen_by_species$tajimas_d)), ]

datasets <- c('Ellegaard2020', 'Wu2021', 'Ellegaard2019', 'Sun2022', 'Zhang2022')
allele_popgen_by_species$dataset <- NA
for (dataset_id in datasets) {
  dataset_samples <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/', dataset_id, '_SRRs.txt.gz', sep = ''),
                                header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1
  allele_popgen_by_species[which(allele_popgen_by_species$sample %in% dataset_samples), 'dataset'] <- dataset_id
}

all_species <- gsub('.tsv.gz', '', list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                                              pattern = '.tsv.gz'))
glmm_COG_category_full_out <- list()
glmm_COG_full_out <- list()
glmm_COG_category_summary_out <- list()
glmm_COG_summary_out <- list()
glmm_COG_category_R2_out <- list()
glmm_COG_R2_out <- list()

for (sp in all_species) {
  tab_subset <- allele_popgen_by_species[which(allele_popgen_by_species$species == sp), ]
  
  # Remove samples with fewer than 50 unique genes.
  samples_by_genes <- tab_subset[, c('dataset', 'sample', 'gene')]
  samples_by_genes <- samples_by_genes[which(! duplicated(samples_by_genes)), ]
  genes_per_sample <- table(samples_by_genes$sample)
  samples_w_sufficient_genes <- names(genes_per_sample)[which(genes_per_sample >= 50)]
  tab_subset <- tab_subset[which(tab_subset$sample %in% samples_w_sufficient_genes), ]
  
  if (nrow(tab_subset) < 100) { next }
  
  tab_subset_COG_category <- split_multi_category_rows(in_df = tab_subset,
                                                       category_col = 'COG_category',
                                                       delimiter = ',',
                                                       num_cores = 40)
  COG_category_tallies <- table(tab_subset_COG_category$COG_category)
  higher_freq_COG_categories <- names(COG_category_tallies)[which(COG_category_tallies >= 20)]
  tab_subset_COG_category <- tab_subset_COG_category[which(tab_subset_COG_category$COG_category %in% higher_freq_COG_categories), , drop = FALSE]
  
  all_COG_category_levels <- unique(tab_subset_COG_category$COG_category)
  tab_subset_COG_category$COG_category <- factor(tab_subset_COG_category$COG_category,
                                                 levels = c('K', all_COG_category_levels[-which(all_COG_category_levels == 'K')]))
  
  tab_subset_COG <- split_multi_category_rows(in_df = tab_subset,
                                              category_col = 'COG',
                                              delimiter = ',',
                                              num_cores = 40)
  COG_tallies <- table(tab_subset_COG$COG)
  higher_freq_COGs <- names(COG_tallies)[which(COG_tallies >= 5)]
  tab_subset_COG <- tab_subset_COG[which(tab_subset_COG$COG %in% higher_freq_COGs), , drop = FALSE]
  
  tab_subset_COG_category$tajimas_d_norm <- bestNormalize::orderNorm(tab_subset_COG_category$tajimas_d)$x
  tab_subset_COG_category$num_seqs_norm <- bestNormalize::orderNorm(tab_subset_COG_category$num_seqs)$x
  tab_subset_COG_category$mean_num_reads_norm <- bestNormalize::orderNorm(tab_subset_COG_category$mean_num_reads)$x
  
  tab_subset_COG$tajimas_d_norm <- bestNormalize::orderNorm(tab_subset_COG$tajimas_d)$x
  tab_subset_COG$num_seqs_norm <- bestNormalize::orderNorm(tab_subset_COG$num_seqs)$x
  tab_subset_COG$mean_num_reads_norm <- bestNormalize::orderNorm(tab_subset_COG$mean_num_reads)$x
  
  num_datasets <- length(unique(tab_subset_COG_category$dataset))
  num_samples <- length(unique(tab_subset_COG_category$sample))
  num_genes <- length(unique(tab_subset_COG_category$gene))
  num_unique_num_seqs <- length(unique(tab_subset_COG_category$num_seqs))
  num_unique_mean_num_reads <- length(unique(tab_subset_COG_category$mean_num_reads))

  random_effect_levels <- as.character()

  if (num_datasets > 1) { random_effect_levels <- c(random_effect_levels, '(1 | dataset)') }
  if (num_samples > 1) { random_effect_levels <- c(random_effect_levels, '(1 | sample)') }
  if (num_genes > 1) { random_effect_levels <- c(random_effect_levels, '(1 | gene)') }
  if (num_unique_num_seqs > 1) { random_effect_levels <- c(random_effect_levels, '(1 | num_seqs_norm)') }
  if (num_unique_mean_num_reads > 1) { random_effect_levels <- c(random_effect_levels, '(1 | mean_num_reads_norm)') }
  
  random_effect_formula_section <- paste(random_effect_levels, collapse = ' + ')
  
  formula_COG_category <- as.formula(paste('tajimas_d_norm ~ COG_category', random_effect_formula_section, sep = ' + '))
  glmm_COG_category_full_out[[sp]] <- glmmTMB::glmmTMB(formula = formula_COG_category,
                                                       data = tab_subset_COG_category,
                                                       family = 'gaussian',
                                                       control = glmmTMBControl(optimizer = nlminb,
                                                                                parallel = 40,
                                                                                profile = TRUE,
                                                                                optCtrl = list(iter.max = 1000,
                                                                                               eval.max = 1000)))
  glmm_COG_category_summary_out[[sp]] <- summary(glmm_COG_category_full_out[[sp]])
  glmm_COG_category_R2_out[[sp]] <- performance::r2(glmm_COG_category_full_out[[sp]])

  formula_COG <- as.formula(paste('tajimas_d_norm ~ COG', random_effect_formula_section, sep = ' + '))
  glmm_COG_full_out[[sp]] <- glmmTMB::glmmTMB(formula = formula_COG,
                                              data = tab_subset_COG,
                                              family = 'gaussian',
                                              control = glmmTMBControl(optimizer = nlminb,
                                                                       parallel = 40,
                                                                       profile = TRUE,
                                                                       optCtrl = list(iter.max = 1000,
                                                                                      eval.max = 1000)))
  glmm_COG_summary_out[[sp]] <- summary(glmm_COG_full_out[[sp]])
  glmm_COG_R2_out[[sp]] <- performance::r2(glmm_COG_full_out[[sp]])

}

saveRDS(object = glmm_COG_category_full_out, file = '/data1/gdouglas/tmp/glmm_COG_category_full_by_species.rds')
saveRDS(object = glmm_COG_category_summary_out, file = '/data1/gdouglas/tmp/glmm_COG_category_summary_by_species.rds')
saveRDS(object = glmm_COG_category_R2_out, file = '/data1/gdouglas/tmp/glmm_COG_category_R2_by_species.rds')

saveRDS(object = glmm_COG_full_out, file = '/data1/gdouglas/tmp/glmm_COG_full_by_species.rds')
saveRDS(object = glmm_COG_summary_out, file = '/data1/gdouglas/tmp/glmm_COG_summary_by_species.rds')
saveRDS(object = glmm_COG_R2_out, file = '/data1/gdouglas/tmp/glmm_COG_R2_by_species.rds')
