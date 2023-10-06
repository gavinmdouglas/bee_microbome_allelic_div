rm(list = ls(all.names = TRUE))

# Compute Tajima's D based on output of strain & allele sequences and presence/absence.
# Compute based on just each unique sequences (without rel. abun.), and also based on
# relative abundance along with the mean number of reads mapped to the gene overall
# per sample.
depth_folder <- '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/mean_depth_per_site/'
mean_depth_tables <- list()
for (depth_file in list.files(depth_folder, full.names = TRUE, pattern = '.mean.bedGraph.gz')) {
  depth_sample_id <- gsub('.mean.bedGraph.gz', '', basename(depth_file))
  mean_depth_tables[[depth_sample_id]] <- read.table(file = depth_file,
                                                     header = FALSE,
                                                     sep = '\t',
                                                     stringsAsFactors = FALSE,
                                                     row.names = 1)[, 3, drop = FALSE]
}

core_genes <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/core_genes.singletons.above_len.rds')

# Return n, S, Watterson's theta, nucleotide diff for a list of input sequences.
strain_tajimas_d_from_relabun <- function(species_id, diff_df, relabun_df) {
  
  raw_popgen_out <- list()
  
  for (sample_id in colnames(relabun_df)) {

    mean_rounded_depth <- ceiling(mean(mean_depth_tables[[sample_id]][core_genes[[species_id]], 1]))
    if (mean_rounded_depth < 10) { next }

    sample_seqs <- rownames(relabun_df)[which(relabun_df[, sample_id] > 0)]
    if (length(sample_seqs) == 1) { next }
    diff_df_subset <- diff_df[which(diff_df$seq1 %in% sample_seqs & diff_df$seq2 %in% sample_seqs), ]
    seq_relabun <- relabun_df[sample_seqs, sample_id, drop = TRUE] / 100
    names(seq_relabun) <- sample_seqs
    seq_weighted <- round(seq_relabun * mean_rounded_depth)
    if (length(which(seq_weighted == 0)) > 0) {
      seq_weighted[which(seq_weighted == 0)] <- 1 
    }

    diff_df_subset$subsamples <- NA
    all_segregating_sites <- as.character()
    for (i in 1:nrow(diff_df_subset)) {
      diff_df_subset[i, 'subsamples'] <- seq_weighted[diff_df_subset[i, 'seq1']] * seq_weighted[diff_df_subset[i, 'seq2']]
      all_segregating_sites <- unique(c(all_segregating_sites, strsplit(x = diff_df_subset[i, 'diff_positions'], split = ',')[[1]]))
    }
    
    multi_count_seqs <- seq_weighted[which(seq_weighted > 1)]
    multi_count_seqs_tallies <- data.frame(seq1 = names(multi_count_seqs),
                                           seq2 = names(multi_count_seqs),
                                           num_comparable_sites = 1,
                                           num_diff = 0,
                                           diff_positions = '',
                                           subsamples = as.integer(sapply(multi_count_seqs, choose, 2)))
    diff_df_subset <- rbind(diff_df_subset, multi_count_seqs_tallies)
    
    num_seg_sites <- length(all_segregating_sites)
    mean_num_diff <- sum((diff_df_subset$subsamples / sum(diff_df_subset$subsamples)) * diff_df_subset$num_diff)
    mean_num_diff_per_site <- mean_num_diff / mean(diff_df_subset$num_comparable_sites)
    prop_segregating_sites <- num_seg_sites / mean(diff_df_subset$num_comparable_sites)

    a1 <- sum(sapply(1:mean_rounded_depth, function(x) { 1 / x} ))
    a2 <- sum(sapply(1:mean_rounded_depth, function(x) { 1 / (x ** 2)} ))
    b1 <- (mean_rounded_depth + 1) / (3 * (mean_rounded_depth - 1))
    b2 <- (2 * (mean_rounded_depth ** 2 + mean_rounded_depth + 3)) / (9 * mean_rounded_depth * (mean_rounded_depth - 1))
    c1 <- b1 - (1 / a1)
    c2 <- b2 - ((mean_rounded_depth + 2) / (a1 * mean_rounded_depth)) + (a2 / a1 ** 2)
    e1 <- c1 / a1
    e2 <- c2 / (a1 ** 2 + a2)
    
    expected_sd = sqrt(e1 * num_seg_sites + e2 * num_seg_sites * (num_seg_sites - 1))
    wattersons_theta = num_seg_sites / a1
    
    if (expected_sd == 0) {
      tajimas_d = NA
    } else {
      tajimas_d <- (mean_num_diff - wattersons_theta) / expected_sd
    }
    
    raw_popgen_out[[sample_id]] <- data.frame(species = species_id,
                                              sample = sample_id,
                                              num_seqs = length(sample_seqs),
                                              mean_num_reads = mean_rounded_depth,
                                              pi = mean_num_diff,
                                              pi_per_site = mean_num_diff_per_site,
                                              num_seg_sites = num_seg_sites,
                                              prop_seg_sites = prop_segregating_sites,
                                              wattersons_theta  = wattersons_theta,
                                              tajimas_d = tajimas_d)
  }
  return(do.call(rbind, raw_popgen_out))
}


allele_tajimas_d_from_relabun <- function(species_id, gene_id, diff_df, relabun_df) {
  
  diff_df$diff_positions <- as.character(diff_df$diff_positions)
  
  raw_popgen_out <- list()
  
  for (sample_id in unique(relabun_df$sample)) {
    
    mean_rounded_depth <- ceiling(mean_depth_tables[[sample_id]][gene_id, 1])
    if (mean_rounded_depth < 10) { next }
    
    relabun_df_subset <- relabun_df[which(relabun_df$sample == sample_id), , drop = FALSE]
    
    sample_seqs <- paste('allele', relabun_df_subset$strain, sep = '_')
    if (length(sample_seqs) == 1) { next }
    diff_df_subset <- diff_df[which(diff_df$seq1 %in% sample_seqs & diff_df$seq2 %in% sample_seqs), ]
    seq_relabun <- relabun_df_subset$community / 100
    names(seq_relabun) <- sample_seqs
    seq_weighted <- round(seq_relabun * mean_rounded_depth)
    if (length(which(seq_weighted == 0)) > 0) {
      seq_weighted[which(seq_weighted == 0)] <- 1 
    }
    
    diff_df_subset$subsamples <- NA
    all_segregating_sites <- as.character()
    for (i in 1:nrow(diff_df_subset)) {
      diff_df_subset[i, 'subsamples'] <- seq_weighted[diff_df_subset[i, 'seq1']] * seq_weighted[diff_df_subset[i, 'seq2']]
      
      if (! is.na(diff_df_subset[i, 'diff_positions'])) {
        all_segregating_sites <- unique(c(all_segregating_sites, strsplit(x = diff_df_subset[i, 'diff_positions'], split = ',')[[1]]))
      }
    }
    
    multi_count_seqs <- seq_weighted[which(seq_weighted > 1)]
    multi_count_seqs_tallies <- data.frame(seq1 = names(multi_count_seqs),
                                           seq2 = names(multi_count_seqs),
                                           num_comparable_sites = 1,
                                           num_diff = 0,
                                           diff_positions = '',
                                           subsamples = as.integer(sapply(multi_count_seqs, choose, 2)))
    diff_df_subset <- rbind(diff_df_subset, multi_count_seqs_tallies)
    
    num_seg_sites <- length(all_segregating_sites)
    mean_num_diff <- sum((diff_df_subset$subsamples / sum(diff_df_subset$subsamples)) * diff_df_subset$num_diff)
    mean_num_diff_per_site <- mean_num_diff / mean(diff_df_subset$num_comparable_sites)
    prop_segregating_sites <- num_seg_sites / mean(diff_df_subset$num_comparable_sites)
    
    a1 <- sum(sapply(1:mean_rounded_depth, function(x) { 1 / x} ))
    a2 <- sum(sapply(1:mean_rounded_depth, function(x) { 1 / (x ** 2)} ))
    b1 <- (mean_rounded_depth + 1) / (3 * (mean_rounded_depth - 1))
    b2 <- (2 * (mean_rounded_depth ** 2 + mean_rounded_depth + 3)) / (9 * mean_rounded_depth * (mean_rounded_depth - 1))
    c1 <- b1 - (1 / a1)
    c2 <- b2 - ((mean_rounded_depth + 2) / (a1 * mean_rounded_depth)) + (a2 / a1 ** 2)
    e1 <- c1 / a1
    e2 <- c2 / (a1 ** 2 + a2)
    
    expected_sd = sqrt(e1 * num_seg_sites + e2 * num_seg_sites * (num_seg_sites - 1))
    wattersons_theta = num_seg_sites / a1
    
    if (expected_sd == 0) {
      tajimas_d = NA
    } else {
      tajimas_d <- (mean_num_diff - wattersons_theta) / expected_sd
    }
    
    raw_popgen_out[[sample_id]] <- data.frame(species = species_id,
                                              sample = sample_id,
                                              gene = gene_id,
                                              num_seqs = length(sample_seqs),
                                              mean_num_reads = mean_rounded_depth,
                                              pi = mean_num_diff,
                                              pi_per_site = mean_num_diff_per_site,
                                              num_seg_sites = num_seg_sites,
                                              prop_seg_sites = prop_segregating_sites,
                                              wattersons_theta  = wattersons_theta,
                                              tajimas_d = tajimas_d)
  }
  return(do.call(rbind, raw_popgen_out))
}
all_species <- gsub('.tsv.gz', '', list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                                              pattern = '.tsv.gz'))
raw_sp_popgen_out <- list()
for (sp in all_species) {

  # Compute stats for strains.
  relabun <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_rel_abun/',
                              sp, '.tsv.gz', sep = ''), header = TRUE, sep = '\t', row.names = 1)
  diff <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/inferred_core_genome_data/pairwise_num_diff/',
                           sp, '.strains.tsv.gz', sep = ''), header = TRUE, sep = '\t')
  raw_sp_popgen_out[[sp]] <- strain_tajimas_d_from_relabun(species_id = sp,
                                                           diff_df = diff,
                                                           relabun_df = relabun)
}


strain_popgen_by_species <- do.call(rbind, raw_sp_popgen_out)
rownames(strain_popgen_by_species) <- NULL
saveRDS(object = strain_popgen_by_species, file = '/home/gdouglas/tmp/strain_popgen_by_species.rds')


raw_gene_popgen_out <- list()
allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')
for (sp in all_species) {
  
  raw_gene_popgen_out[[sp]] <- parallel::mclapply(X = names(allele_relabun[[sp]]),
                                            function(gene_id) {
                                              
                                              gene_diff_filename <- paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/allele_analyses/pairwise_num_diff/',
                                                                          gene_id, '.tsv.gz', sep = '')
                                              if (file.exists(gene_diff_filename)) {
                                                gene_diff <- read.table(gene_diff_filename, header = TRUE, sep = '\t')
                                                allele_tajimas_d_from_relabun(species_id = sp,
                                                                              gene_id = gene_id,
                                                                              diff_df = gene_diff,
                                                                              relabun_df = allele_relabun[[sp]][[gene_id]])
                                              }
                                            },
                                            mc.cores = 40)
}

raw_gene_popgen_out_combined <- list()
for (sp in all_species) {
  raw_gene_popgen_out_combined[[sp]] <- do.call(rbind, raw_gene_popgen_out[[sp]])
}

allele_popgen_by_species <- do.call(rbind, raw_gene_popgen_out_combined)
rownames(allele_popgen_by_species) <- NULL


saveRDS(object = allele_popgen_by_species, file = '/home/gdouglas/tmp/allele_popgen_by_species.rds')
