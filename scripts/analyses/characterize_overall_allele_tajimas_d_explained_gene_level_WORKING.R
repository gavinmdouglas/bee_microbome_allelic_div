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

which.median <- function(x) which.min(abs(x - median(x)))

allele_popgen_by_species <- readRDS(file = '/home/gdouglas/tmp/allele_popgen_by_species.rds')

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
glmm_gene_full_out <- list()
glmm_gene_summary_out <- list()
glmm_gene_R2_out <- list()

for (sp in all_species) {
  tab_subset <- allele_popgen_by_species[which(allele_popgen_by_species$species == sp), ]
  
  # Remove samples with fewer than 50 unique genes.
  samples_by_genes <- tab_subset[, c('dataset', 'sample', 'gene')]
  samples_by_genes <- samples_by_genes[which(! duplicated(samples_by_genes)), ]
  genes_per_sample <- table(samples_by_genes$sample)
  samples_w_sufficient_genes <- names(genes_per_sample)[which(genes_per_sample >= 50)]
  tab_subset <- tab_subset[which(tab_subset$sample %in% samples_w_sufficient_genes), ]
  
  # Remove genes that did not pass cut-offs in at least three samples.
  samples_per_gene <- table(tab_subset$gene)
  genes_w_sufficient_samples <- names(samples_per_gene)[which(samples_per_gene >= 3)]
  tab_subset <- tab_subset[which(tab_subset$gene %in% genes_w_sufficient_samples), ]

  if (nrow(tab_subset) < 100) { next }
  
  tab_subset$tajimas_d_norm <- bestNormalize::orderNorm(tab_subset$tajimas_d)$x

  num_datasets <- length(unique(tab_subset$dataset))
  num_samples <- length(unique(tab_subset$sample))

  random_effect_levels <- as.character()

  if (num_datasets > 1) { random_effect_levels <- c(random_effect_levels, '(1 | dataset)') }
  if (num_samples > 1) { random_effect_levels <- c(random_effect_levels, '(1 | sample)') }

  random_effect_formula_section <- paste(random_effect_levels, collapse = ' + ')

  # Choose the gene with the median normalized tajima's d norm value as the reference category.
  mean_D_by_gene <- aggregate(x = tajimas_d_norm ~ gene, data = tab_subset, FUN=mean)
  ref_gene_level <- mean_D_by_gene[which.median(mean_D_by_gene$tajimas_d_norm), 'gene']
  all_gene_levels <- unique(tab_subset$gene)
  tab_subset$gene <- factor(tab_subset$gene,
                            levels = c(ref_gene_level, all_gene_levels[-which(all_gene_levels == ref_gene_level)]))
  
  formula_gene <- as.formula(paste('tajimas_d_norm ~ gene', random_effect_formula_section, sep = ' + '))
  glmm_gene_full_out[[sp]] <- glmmTMB::glmmTMB(formula = formula_gene,
                                              data = tab_subset,
                                              se = FALSE,
                                              family = 'gaussian',
                                              control = glmmTMBControl(optimizer = nlminb,
                                                                       parallel = 40,
                                                                       profile = TRUE,
                                                                       optCtrl = list(iter.max = 1000,
                                                                                      eval.max = 1000)))
  glmm_gene_summary_out[[sp]] <- summary(glmm_gene_full_out[[sp]])
  glmm_gene_R2_out[[sp]] <- performance::r2(glmm_gene_full_out[[sp]])

}

saveRDS(object = glmm_gene_full_out, file = '/data1/gdouglas/tmp/glmm_gene_full_by_species.rds')
saveRDS(object = glmm_gene_summary_out, file = '/data1/gdouglas/tmp/glmm_gene_summary_by_species.rds')
saveRDS(object = glmm_gene_R2_out, file = '/data1/gdouglas/tmp/glmm_gene_R2_by_species.rds')
