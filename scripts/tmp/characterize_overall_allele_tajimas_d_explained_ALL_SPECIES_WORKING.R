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

# Remove samples with fewer than 500 unique genes.
samples_by_genes <- allele_popgen_by_species[, c('dataset', 'sample', 'gene')]
samples_by_genes <- samples_by_genes[which(! duplicated(samples_by_genes)), ]
genes_per_sample <- table(samples_by_genes$sample)
samples_w_sufficient_genes <- names(genes_per_sample)[which(genes_per_sample >= 500)]
allele_popgen_by_species <- allele_popgen_by_species[which(allele_popgen_by_species$sample %in% samples_w_sufficient_genes), ]

allele_popgen_by_species_COG_category <- split_multi_category_rows(in_df = allele_popgen_by_species,
                                                     category_col = 'COG_category',
                                                     delimiter = ',',
                                                     num_cores = 40)
COG_category_tallies <- table(allele_popgen_by_species_COG_category$COG_category)
higher_freq_COG_categories <- names(COG_category_tallies)[which(COG_category_tallies >= 100)]
allele_popgen_by_species_COG_category <- allele_popgen_by_species_COG_category[which(allele_popgen_by_species_COG_category$COG_category %in% higher_freq_COG_categories), , drop = FALSE]

all_COG_category_levels <- unique(allele_popgen_by_species_COG_category$COG_category)
allele_popgen_by_species_COG_category$COG_category <- factor(allele_popgen_by_species_COG_category$COG_category,
                                               levels = c('K', all_COG_category_levels[-which(all_COG_category_levels == 'K')]))

all_species_levels <- unique(allele_popgen_by_species_COG_category$species)
allele_popgen_by_species_COG_category$species <- factor(allele_popgen_by_species_COG_category$species,
                                                        levels = c('Lactobacillus_melliventris', all_species_levels[-which(all_species_levels == 'Lactobacillus_melliventris')]))

allele_popgen_by_species_COG_category$tajimas_d_norm <- bestNormalize::orderNorm(allele_popgen_by_species_COG_category$tajimas_d)$x
allele_popgen_by_species_COG_category$num_seqs_norm <- bestNormalize::orderNorm(allele_popgen_by_species_COG_category$num_seqs)$x
allele_popgen_by_species_COG_category$mean_num_reads_norm <- bestNormalize::orderNorm(allele_popgen_by_species_COG_category$mean_num_reads)$x

random_effect_formula_section <- '(1 | dataset) + (1 | sample) + (1 | gene) + (1 | num_seqs_norm) + (1 | mean_num_reads_norm)'

formula_COG_category <- as.formula(paste('tajimas_d_norm ~ species + COG_category', random_effect_formula_section, sep = ' + '))
glmm_COG_category_full_out <- glmmTMB::glmmTMB(formula = formula_COG_category,
                                                     data = allele_popgen_by_species_COG_category,
                                                     family = 'gaussian',
                                                     control = glmmTMBControl(optimizer = nlminb,
                                                                              parallel = 40,
                                                                              profile = TRUE,
                                                                              optCtrl = list(iter.max = 1000,
                                                                                             eval.max = 1000)))
glmm_COG_category_summary_out <- summary(glmm_COG_category_full_out)
glmm_COG_category_R2_out <- performance::r2(glmm_COG_category_full_out)

saveRDS(object = glmm_COG_category_full_out, file = '/data1/gdouglas/tmp/glmm_COG_category_full_ALL_species.rds')
saveRDS(object = glmm_COG_category_summary_out, file = '/data1/gdouglas/tmp/glmm_COG_category_summary_ALL_species.rds')
saveRDS(object = glmm_COG_category_R2_out, file = '/data1/gdouglas/tmp/glmm_COG_category_R2_ALL_species.rds')
