raw_to_short <- list(
  'Apilactobacillus_apinorum' = 'A. apinorum',
  'Apilactobacillus_kunkeei' = 'A. kunkeei',
  'Bartonella_apis' = 'Bart. apis',
  'Bifidobacterium_asteroides' = 'B. asteroides',
  'Bifidobacterium_coryneforme_indicum' = 'B. cor./indicum',
  'Bombella_apis' = 'Bom. apis',
  'Bombella_sp' = 'B. sp.',
  'Bombilactobacillus_mellifer' = 'B. mellifer', 
  'Bombilactobacillus_mellis' = 'B. mellis',
  'Commensalibacter_sp' = 'C. sp.',
  'Frischella_perrara' = 'F. perrara',
  'Gilliamella_apicola' = 'G. apicola',
  'Gilliamella_apis' = 'G. apis',
  'Gilliamella_sp' = 'G. sp.',
  'Lactobacillus_apis' = 'L. apis',
  'Lactobacillus_helsingborgensis' = 'L. helsingborgensis',
  'Lactobacillus_kimbladii' = 'L. kimbladii',
  'Lactobacillus_kullabergensis' = 'L. kullabergensis',
  'Lactobacillus_melliventris' = 'L. melliventris',
  'Serratia_marcescens' = 'S. marcescens',
  'Snodgrassella_alvi' = 'S. alvi'
)

cleanup_species_names <- function(names_vec, shorten=FALSE) {
  
  if (shorten) {
    for (raw in names(raw_to_short)) {
      if (raw %in% names_vec) {
        names_vec[which(names_vec == raw)] <- raw_to_short[[raw]]
      }
    }
    
    if ('Bifidobacterium_coryneforme' %in% names_vec) {
      names_vec[which(names_vec == 'Bifidobacterium_coryneforme')] <- raw_to_short[['Bifidobacterium_coryneforme_indicum']]
    }
    
  } else {
    
    names_vec <- gsub('_', ' ', names_vec)
    if (length(grep('Bifidobacterium coryneforme', names_vec)) > 0) {
      names_vec[grep('Bifidobacterium coryneforme', names_vec)] <- 'Bifidobacterium cor./indicum'
    }
    
    names_vec <- gsub('sp$', 'sp.', names_vec)
    
  }
  
  return(names_vec)
  
}


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

species_plot_colours <- c(
  'pink1', # 'Apilactobacillus apinorum'
  'pink2', # 'Apilactobacillus kunkeei'
  'yellow2', # 'Bartonella apis'
  '#63bfaf', # 'Bifidobacterium asteroides'
  '#62d07a', # 'Bifidobacterium cor./indicum'
  'cyan', # 'Bombella apis'
  'cyan2', # 'Bombella sp.'
  '#d56c24', # 'Bombilactobacillus mellifer'
  '#bf4d81', # 'Bombilactobacillus mellis'
  '#d2cb3b', # 'Commensalibacter sp.'          
  'brown', #'Frischella perrara'
  'royalblue1', # 'Gilliamella apicola'
  'skyblue1', # 'Gilliamella apis'
  'mediumblue', # 'Gilliamella sp.'
  'grey40', # 'Lactobacillus apis'            
  'grey50', #'Lactobacillus helsingborgensis'
  'grey60', # 'Lactobacillus kimbladii'
  'grey70', # 'Lactobacillus kullabergensis'
  'grey80', # 'Lactobacillus melliventris'
  'black', # 'Serratia marcescens'           
  '#9f4ad4' # 'Snodgrassella alvi'
)

names(species_plot_colours) <- c('Apilactobacillus_apinorum',
                                 'Apilactobacillus_kunkeei',
                                 'Bartonella_apis',
                                 'Bifidobacterium_asteroides',
                                 'Bifidobacterium_coryneforme_indicum',
                                 'Bombella_apis',
                                 'Bombella_sp',
                                 'Bombilactobacillus_mellifer',
                                 'Bombilactobacillus_mellis',
                                 'Commensalibacter_sp',
                                 'Frischella_perrara',
                                 'Gilliamella_apicola',
                                 'Gilliamella_apis',
                                 'Gilliamella_sp',
                                 'Lactobacillus_apis',
                                 'Lactobacillus_helsingborgensis',
                                 'Lactobacillus_kimbladii',
                                 'Lactobacillus_kullabergensis',
                                 'Lactobacillus_melliventris',
                                 'Serratia_marcescens',
                                 'Snodgrassella_alvi')

identify_enriched_categories <- function(genes,
                                         background,
                                         gene_to_category_map,
                                         min_category_count = 0,
                                         to_ignore = character()) {
  

  
  enrichments_out <- data.frame(matrix(NA, nrow = length(gene_to_category_map), ncol = 8))
  rownames(enrichments_out) <- names(gene_to_category_map)
  colnames(enrichments_out) <- c("category", "genes_num_category", "genes_num_other",
                                 "background_num_category", "background_num_other", "OR", "p", "fdr")
  
  enrichments_out[names(gene_to_category_map), "category"] <- names(gene_to_category_map)
  
  for (category in rownames(enrichments_out)) {
    
    if (category %in% to_ignore) { next }
    
    genes_num_category <- length(which(genes %in% gene_to_category_map[[category]]))
    genes_num_other <- length(genes) - genes_num_category
    
    background_num_category <- length(which(background %in% gene_to_category_map[[category]]))
    background_num_other <- length(background) - background_num_category
    
    count_table <- matrix(c(genes_num_category, genes_num_other, background_num_category, background_num_other), nrow = 2, ncol = 2)
    
    if (min(c(genes_num_category + background_num_category, genes_num_other + background_num_other)) < min_category_count) {
      next
    }
    
    fisher_out <- fisher.test(count_table)
    
    enrichments_out[category, c("genes_num_category",
                                "genes_num_other",
                                "background_num_category",
                                "background_num_other", "p")] <- c(genes_num_category,
                                                                   genes_num_other,
                                                                   background_num_category,
                                                                   background_num_other,
                                                                   fisher_out$p.value)
    if (genes_num_other > 0) {
      ratio_numer <- genes_num_category / genes_num_other
    } else {
      ratio_numer <- genes_num_category / 1 
    }
    
    if (background_num_other == 0) {
      ratio_denom <- 1
    } else if(background_num_category == 0) {
      ratio_denom <- 1 / background_num_other
    } else {
      ratio_denom <- background_num_category / background_num_other
    }
    
    enrichments_out[category, "OR"] <- ratio_numer / ratio_denom
  }
  
  if (length(which(rowSums(is.na(enrichments_out)) > 1)) > 0) {
    enrichments_out <- enrichments_out[-which(rowSums(is.na(enrichments_out)) > 1), ]
  }
  
  enrichments_out$fdr <- p.adjust(enrichments_out$p, "fdr")
  
  rownames(enrichments_out) <- NULL
  
  return(enrichments_out)
  
}

pairwise_matrix_summed_dist <- function(dist_mats) {
  
  if (length(dist_mats) <= 1) {
    return(NULL)
  }
  
  sample_subset <- colnames(as.matrix(dist_mats[[1]]))
  for (i in 2:length(dist_mats)) {
    sample_subset <- intersect(sample_subset, colnames(as.matrix(dist_mats[[i]])))
  }
  
  all_mats <- list()
  for (i in 1:length(dist_mats)) {
    all_mats[[i]] <- as.matrix(dist_mats[[i]])[sample_subset, sample_subset]
  }
  
  summed_dist <- numeric()
  rep_i <- character()
  rep_j <- character()
  for (i in 1:(length(all_mats) - 1)) {
    for (j in (i + 1):length(dist_mats)) {
      summed_dist <- c(summed_dist, sum(abs(all_mats[[i]] - all_mats[[j]])))
      rep_i <- c(rep_i, names(dist_mats)[i])
      rep_j <- c(rep_j, names(dist_mats)[j])
    }
  }
  
  return(data.frame(rep_i=rep_i, rep_j=rep_j, summed_dist=summed_dist))
  
}

read_acc_comm_matrix <- function(file_in, sample_map, freq_cutoff=0.1) {
  raw_tab <- read.table(file_in, header = TRUE, sep = "\t")
  raw_tab <- raw_tab[which(raw_tab$community >= freq_cutoff), ]
  if (nrow(raw_tab) == 0) { return(NULL) }
  raw_tab$sample <- sample_map[as.character(raw_tab$sample), 'V2']
  wide_tab <- reshape2::dcast(data = raw_tab, formula = sample ~ strain, value.var = 'community', fill = 0)
  rownames(wide_tab) <- wide_tab$sample
  wide_tab <- wide_tab[, -1, drop = FALSE]
  wide_tab <- sweep(x = wide_tab, MARGIN = 1, STATS = rowSums(wide_tab), FUN = '/') * 100
  return(wide_tab)
}

identify_centroid_comm_replicate <- function(comm_file_struc) {
  
  gene_comm_files <- list.files(comm_file_path,
                                pattern = basename(comm_file_struc),
                                full.names = TRUE)
  
  # Return NA if there are not at least five replicate files.
  if (length(gene_comm_files) < 5) { return(NA) }
  
  comms <- list()
  
  base_file_split <- strsplit(gsub('_metagenotype.comm.tsv.gz', '', basename(gene_comm_files[1])), '\\.')[[1]]
  dataset <- base_file_split[1]
  gene <- base_file_split[4]
  
  for (gene_comm_file in gene_comm_files) {
    rep <- strsplit(gsub('_metagenotype.comm.tsv.gz', '', basename(gene_comm_file)), '\\.')[[1]][3]
    
    sample_input <- read.table(paste(prepped_path, '/', dataset, '/', rep, '/samples/', gene, '_samples.tsv', sep = ''),
                               sep = "\t", row.names = 1, header = FALSE)
    
    comms[[gene_comm_file]] <- read_acc_comm_matrix(file_in = gene_comm_file, sample_map = sample_input)
  }
  
  comms_dist <- lapply(comms,
                       function(x) {
                         vegan::vegdist(x=x, method="bray")
                       })
  
  pairwise_summary <- pairwise_matrix_summed_dist(comms_dist)
  
  pairwise_summary_long <- data.frame(rep_id = c(pairwise_summary$rep_i, pairwise_summary$rep_j),
                                      summed_dist = c(pairwise_summary$summed_dist, pairwise_summary$summed_dist))
  summed_dist_median <- aggregate(x = summed_dist ~ rep_id, FUN = median, data = pairwise_summary_long)
  
  ## I originally was computing Mantel test statistics, but the summed distance is easy to interpret (and robust across all matrix sizes).
  # 
  # pairwise_mantel_summary <- pairwise_matrix_mantel(comms_dist)
  # pairwise_mantel_summary_long <- data.frame(rep_id = c(pairwise_mantel_summary$rep_i, pairwise_mantel_summary$rep_j),
  #                                            mantel_r = c(pairwise_mantel_summary$mantel_r, pairwise_mantel_summary$mantel_r))
  # mantel_r_median <- aggregate(x = mantel_r ~ rep_id, FUN = median, data = pairwise_mantel_summary_long)
  # mantel_r_median[order(mantel_r_median$mantel_r, decreasing = TRUE), ]
  
  return(summed_dist_median[order(summed_dist_median$summed_dist, decreasing = FALSE), ][1, 'rep_id'])
}

preprocess_subsampled_comm <- function(comm_file) {
  
  comm_file_split <- strsplit(comm_file, '\\.')[[1]]
  
  dataset <- basename(comm_file_split[1])
  subsample <- comm_file_split[2]
  rep <- comm_file_split[3]
  gene <- sub('_metagenotype$', '', comm_file_split[4])
  
  if (subsample == "subsample50") {
    gene_samples_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input_subsampled/',
                               dataset,
                               '/',
                               rep,
                               '/samples/',
                               gene,
                               '_samples.tsv',
                               sep = '')
  } else if (subsample == "subsample20") {
    gene_samples_file <- paste('/scratch/gdouglas/projects/honey_bee/strainfacts_working/prepped_input/prepped_accessory_input_subsampled20/',
                               dataset,
                               '/',
                               rep,
                               '/samples/',
                               gene,
                               '_samples.tsv',
                               sep = '')
  }
  
  gene_samples <- read.table(gene_samples_file, sep = "\t", row.names = 1, header = FALSE)
  
  allele_freq <- read.table(comm_file, header = TRUE, sep = "\t")
  
  if (length(unique(allele_freq$sample)) != nrow(gene_samples)) {
    stop('Mismatch in expected number of samples and those in sample mapfile!') 
  }
  
  allele_freq$sample <- as.character(allele_freq$sample)
  
  if (length(which(! allele_freq$sample %in% rownames(gene_samples))) > 0) {
    stop('Error - sample missing in mapfile.') 
  }
  
  allele_freq$sample <- gene_samples[allele_freq$sample, 1]
  
  orig_samples <- sort(unique(allele_freq$sample))
  
  # Only keep strain calls that were at least >= 10%
  allele_freq <- allele_freq[which(allele_freq$community >= 0.1), ]
  
  if (! identical(orig_samples, sort(unique(allele_freq$sample)))) {
    message("WARNING - samples dropped after excluding low abun strains.")
  }
  
  # Total-sum scale by sample remaining rel. abun.
  sample_sums <- aggregate(x = community ~ sample,
                           data = allele_freq,
                           FUN = sum)
  rownames(sample_sums) <- sample_sums$sample
  allele_freq$community <- (allele_freq$community / sample_sums[allele_freq$sample, "community"]) * 100
  
  return(allele_freq)
  
}

