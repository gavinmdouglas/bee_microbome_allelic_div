rm(list = ls(all.names = TRUE))

# Summarize single nucleotide variants from mapped reads per gene.
# Create simplified tables of the mean within/between values per gene.

create_mean_tab <- function(in_tab) {
  in_tab_mean_pi <- aggregate(x = pi ~ gene, data = in_tab, FUN = mean)
  in_tab_mean_nseg <- aggregate(x = num_seg ~ gene, data = in_tab, FUN = mean)
  in_tab_mean_num_sites <- aggregate(x = num_sites ~ gene, data = in_tab, FUN = mean)
  
  rownames(in_tab_mean_pi) <- in_tab_mean_pi$gene
  rownames(in_tab_mean_nseg) <- in_tab_mean_nseg$gene
  rownames(in_tab_mean_num_sites) <- in_tab_mean_num_sites$gene
  
  summary_df <- data.frame(gene = in_tab_mean_pi$gene,
                           mean_pi = in_tab_mean_pi$pi)
  rownames(summary_df) <- summary_df$gene
  
  summary_df[rownames(in_tab_mean_nseg), 'mean_nseg'] <- in_tab_mean_nseg$num_seg
  
  summary_df[rownames(in_tab_mean_num_sites), 'mean_num_sites'] <- in_tab_mean_num_sites$num_sites
  
  rownames(summary_df) <- NULL
  
  return(summary_df)
}

datasets <- c('Ellegaard_2019', 'Ellegaard_2020', 'Sun_2022', 'Wu_2021', 'Zhang_2022')

min_sites <- 200

raw_pi <- list()

for (d in datasets) {
 
  print(d)
  
  pi_file <- paste('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/bcftools_out_merged_pi/',
                   d,
                   '.summary.tsv.gz',
                   sep = '')
  
  pi_tab <- read.table(file = pi_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  pi_tab <- pi_tab[which(pi_tab$num_sites >= min_sites), ]
  
  colnames(pi_tab)[which(colnames(pi_tab) == 'contig')] <- 'gene'
  
  pi_tab_within <- pi_tab[which(pi_tab$comparison_type == 'within'), ]
  pi_tab_within_mean <- create_mean_tab(in_tab = pi_tab_within)
  pi_tab_within_mean$type <- 'Within'
  
  pi_tab_between <- pi_tab[which(pi_tab$comparison_type == 'between'), ]
  pi_tab_between_mean <- create_mean_tab(in_tab = pi_tab_between)
  pi_tab_between_mean$type <- 'Between'

  combined_raw <- rbind(pi_tab_within_mean, pi_tab_between_mean)
  
  combined_raw$dataset <- d
  
  raw_pi[[d]] <- combined_raw
  
}

combined_tab <- do.call(rbind, raw_pi)

write.table(x = combined_tab,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/mean_pi_and_nseg_per_gene.tsv',
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
