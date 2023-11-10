rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)

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

species <- gsub('.tsv.gz', '',
                list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                           pattern = '*.tsv.gz'))

overall_mean <- data.frame(species = species, pi = NA, dnds = NA)
rownames(overall_mean) <- species

raw_per_sample_mean <- list()

for (sp in species) {

  relabun <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                              sp, '.tsv.gz', sep = ''), header = TRUE)
  
  dnds <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/inferred_core_genome_data/pairwise_dnds/',
                           sp, '.tsv.gz', sep = ''), header = TRUE)
  
  diff <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/inferred_core_genome_data/pairwise_num_diff/',
                           sp, '.strains.tsv.gz', sep = ''), header = TRUE)
  
  overall_mean[sp, 'pi'] <- mean(diff$num_diff / diff$num_comparable_sites)
  overall_mean[sp, 'dnds'] <- mean(dnds$dnds)
  
  tmp <- data.frame(species = sp,
                    sample = colnames(relabun),
                    pi = NA,
                    dnds = NA,
                    num_strains = NA)
  rownames(tmp) <- colnames(relabun)

  for (samp in colnames(relabun)) {

    present_strains <- rownames(relabun)[which(relabun[, samp] > 0)]

    if (length(present_strains) == 1) { 
      tmp[samp, c('pi', 'dnds', 'num_strains')] <- c(NA, NA, NA)
      next
    }

    dnds_subset <- dnds[which(dnds$seq1 %in% present_strains & dnds$seq2 %in% present_strains), ]
    diff_subset <- diff[which(diff$seq1 %in% present_strains & diff$seq2 %in% present_strains), ]
    
    tmp[samp, 'dnds'] <- mean(dnds_subset$dnds)
    tmp[samp, 'pi'] <- mean(diff_subset$num_diff / diff_subset$num_comparable_sites)
    tmp[samp, 'num_strains'] <- length(present_strains)
    
  }
  
  raw_per_sample_mean[[sp]] <- tmp

}

per_sample_mean <- do.call(rbind, raw_per_sample_mean)

ggplot(data = per_sample_mean, aes(x = pi, y = species)) +
  geom_quasirandom() +
  geom_boxplot(alpha = 0.5) +
  geom_quasirandom(data = overall_mean, colour = 'red') +
  theme_bw()


# Read in allele data too.
allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

to_ignore <- c('Bartonella_apis_group_69',
               'Gilliamella_apicola_group_712',
               'Snodgrassella_alvi_group_337',
               'Bartonella_apis_metZ_2')

raw_sp_gene_data <- list()
raw_overall_gene_data <- list()

for (sp in names(allele_relabun)) {

  raw_overall_gene_data[[sp]] <- data.frame(species = sp,
                                            gene = names(allele_relabun[[sp]]),
                                            mean_pi = NA,
                                            mean_dnds = NA)
  rownames(raw_overall_gene_data[[sp]]) <- raw_overall_gene_data[[sp]]$gene
  
  for (g in names(allele_relabun[[sp]])) {
    
    if (g %in% to_ignore) { next }
    
    multi_strain_samples <- sort(unique(allele_relabun[[sp]][[g]]$sample[which(duplicated(allele_relabun[[sp]][[g]]$sample))]))
    
    if (length(multi_strain_samples) == 0) { next }
    
    gene_subset <- allele_relabun[[sp]][[g]][which(allele_relabun[[sp]][[g]]$sample %in% multi_strain_samples), ]
    
    gene_data <- data.frame(species = sp,
                            gene = g,
                            sample = multi_strain_samples,
                            mean_pi = NA,
                            mean_dnds = NA,
                            num_alleles = NA)
    rownames(gene_data) <- multi_strain_samples
    
    diff_filename <- paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/allele_analyses/pairwise_num_diff/',
                           g,
                           '.tsv.gz', sep = '')
    
    dnds_filename <- paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/allele_analyses/pairwise_dnds/',
                           g,
                           '.tsv.gz', sep = '')
    
    dnds <- read.table(dnds_filename, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
    
    diff <- read.table(diff_filename, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
    
    raw_overall_gene_data[[sp]][g, 'mean_dnds'] <- mean(dnds$dnds)
    raw_overall_gene_data[[sp]][g, 'mean_pi'] <- mean(diff$num_diff / diff$num_comparable_sites)
    
    for (samp in rownames(gene_data)) {
      
      present_alleles <- paste('allele', gene_subset[which(gene_subset$sample == samp), 'strain'], sep = '_')
      
      dnds_subset <- dnds[which(dnds$seq1 %in% present_alleles & dnds$seq2 %in% present_alleles), ]
      diff_subset <- diff[which(diff$seq1 %in% present_alleles & diff$seq2 %in% present_alleles), ]
      
      gene_data[samp, 'mean_dnds'] <- mean(dnds_subset$dnds)
      gene_data[samp, 'mean_pi'] <- mean(diff_subset$num_diff / diff_subset$num_comparable_sites)
      gene_data[samp, 'num_alleles'] <- length(present_alleles)
      
    }
    
    raw_sp_gene_data[[g]] <- gene_data
    
  }

}

sp_allele_data <- do.call(rbind, raw_sp_gene_data)
overall_allele_data <- do.call(rbind, raw_overall_gene_data)

sp_allele_data_averaged <- aggregate(x =  mean_pi ~ gene,
                                     data = sp_allele_data,
                                     FUN = mean)
intersecting_genes <- intersect(sp_allele_data_averaged$gene,
                                overall_allele_data$gene)

rownames(sp_allele_data_averaged) <- sp_allele_data_averaged$gene
rownames(overall_allele_data) <- overall_allele_data$gene
sp_allele_data_averaged <- sp_allele_data_averaged[intersecting_genes, ]
overall_allele_data <- overall_allele_data[intersecting_genes, ]

sp_allele_data_averaged$mean_pi_ratio <- sp_allele_data_averaged$mean_pi / overall_allele_data$mean_pi
sp_allele_data_averaged$species <- overall_allele_data[sp_allele_data_averaged$gene, 'species']

ggplot(data = sp_allele_data, aes(y = species, x = mean_pi)) +
  geom_quasirandom(size = 0.1, col = 'grey50') +
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  theme_bw()



# Add in COG category information.
COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")

sp_allele_data$COG_category <- NA
intersecting_genes  <- intersect(sp_allele_data$gene, rownames(COG_category_info))
sp_allele_data <- sp_allele_data[which(sp_allele_data$gene %in% intersecting_genes), ]
sp_allele_data$COG_category <- COG_category_info[sp_allele_data$gene, 'COG_category']
sp_allele_data$COG <- COG_category_info[sp_allele_data$gene, 'all_COG']

sp_allele_data <- split_multi_category_rows(in_df = sp_allele_data,
                                            category_col = 'COG',
                                            delimiter = ',',
                                            num_cores = 40)

sp_allele_data <- sp_allele_data[which(sp_allele_data$COG_category != '-'), ]

test <- sp_allele_data[which(sp_allele_data$species == 'Lactobacillus_melliventris'), ]


test$mean_pi_norm <- log(bestNormalize::orderNorm(test$mean_pi)$x)
test <- test[-which(test$mean_pi_norm == -Inf), ]
tmp <- lm(mean_pi_norm ~ COG, data = test)


test$mean_dnds_norm <- log(bestNormalize::orderNorm(test$mean_dnds)$x)
test <- test[-which(test$mean_dnds_norm == -Inf), ]
tmp2 <- lm(mean_dnds_norm ~ COG_category, data = test)


test <- test[which(test$mean_pi > 0.01), ]

ggplot(data = test, aes(y = COG_category, x = mean_dnds)) +
  geom_quasirandom(size = 0.1, col = 'grey50') +
  geom_boxplot(outlier.shape = NA, alpha=0.8) +
  theme_bw()
  


ggplot(data = sp_allele_data, aes(y = species, x = mean_pi)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(data = overall_mean, aes(x = pi, y = species))


# Ignore rows that could not be annotated with COG categories.
min_jaccard_ORQ <- min_jaccard_ORQ[which(! is.na(min_jaccard_ORQ$COG_category)), ]

sp_allele_data


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

