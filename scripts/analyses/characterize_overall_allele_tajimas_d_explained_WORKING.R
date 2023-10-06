rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)

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

allele_popgen_by_species <- split_multi_category_rows(in_df = allele_popgen_by_species,
                                                      category_col = 'COG_category',
                                                      delimiter = ',',
                                                      num_cores = 40)

all_COG_category_levels <- unique(allele_popgen_by_species$COG_category)
allele_popgen_by_species$COG_category <- factor(allele_popgen_by_species$COG_category,
                                                levels = c('K', all_COG_category_levels[-which(all_COG_category_levels == 'K')]))
allele_popgen_by_species <- allele_popgen_by_species[which(allele_popgen_by_species$mean_num_reads >= 50), ]
allele_popgen_by_species <- allele_popgen_by_species[which(allele_popgen_by_species$num_seqs >= 3), ]

allele_popgen_by_species$tajimas_d_norm <- bestNormalize::orderNorm(allele_popgen_by_species$tajimas_d)$x

glmm_out_w_interact <- glmmTMB::glmmTMB(formula = tajimas_d_norm ~ COG_category + COG_category:species + (1 | species) + (1 | sample) + (1 | gene) + (1 | num_seqs) + (1 | mean_num_reads),
                                        data = allele_popgen_by_species)
glmm_out_w_interact_r2 <- performance::r2(glmm_out_w_interact)
glmm_out_w_interact_summary <- summary(glmm_out_w_interact)

saveRDS(object = glmm_out_w_interact, file = '/data1/gdouglas/tmp/glmm_out_w_interact.rds')
saveRDS(object = glmm_out_w_interact_r2, file = '/data1/gdouglas/tmp/glmm_out_w_interact_r2.rds')
saveRDS(object = glmm_out_w_interact_summary, file = '/data1/gdouglas/tmp/glmm_out_w_interact_summary.rds')


glmm_out_no_interact <- glmmTMB::glmmTMB(formula = tajimas_d_norm ~ COG_category + (1 | species) + (1 | sample) + (1 | gene) + (1 | num_seqs) + (1 | mean_num_reads),
                                        data = allele_popgen_by_species)
glmm_out_no_interact_r2 <- performance::r2(glmm_out_no_interact)
glmm_out_no_interact_summary <- summary(glmm_out_no_interact)

saveRDS(object = glmm_out_no_interact, file = '/data1/gdouglas/tmp/glmm_out_no_interact.rds')
saveRDS(object = glmm_out_no_interact_r2, file = '/data1/gdouglas/tmp/glmm_out_no_interact_r2.rds')
saveRDS(object = glmm_out_no_interact_summary, file = '/data1/gdouglas/tmp/glmm_out_no_interact_summary.rds')
