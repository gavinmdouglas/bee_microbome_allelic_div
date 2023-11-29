rm(list = ls(all.names = TRUE))

library(parallel)

compute_within_sample_mean_allele_divergence <- function(divergence_tab,
                                                         relabun_tab) {

  unique_samples <- unique(relabun_tab$sample)
  
  within_sample_divergence <- numeric()
  for (sampleid in unique_samples) {
    sample_strains <- sort(relabun_tab[which(relabun_tab$sample == sampleid), 'strain'])
    
    if (length(sample_strains) > 1) {
      
      per_within_sample_divergence <- numeric()
      
      for (i in 1:(length(sample_strains) - 1)) {
        
        strain1 <- sample_strains[i]
        
        for (j in (i + 1):length(sample_strains)) {
          
          strain2 <- sample_strains[j]
          
          if (strain1 == strain2) { next }
          
          strain_combo <- paste(strain1, strain2, sep = '_')
          
          per_within_sample_divergence <- c(per_within_sample_divergence,
                                            divergence_tab[strain_combo, 'percent_identity'])
          
        }
        
      }
      
      within_sample_divergence <- c(within_sample_divergence,
                                    mean(per_within_sample_divergence))
    } else {
      
      within_sample_divergence <- c(within_sample_divergence, NA)
      
    }

  }

  return(data.frame(sample = unique_samples,
                    mean_within_divergence = within_sample_divergence))

}


compute_between_sample_mean_allele_divergence <- function(divergence_tab,
                                                          relabun_tab) {

  unique_samples <- unique(relabun_tab$sample)
  
  between_sample_divergence <- numeric()
  sample1_set <- character()
  sample2_set <- character()
  
  if (length(unique_samples) > 1) {
    for (i in 1:(length(unique_samples) - 1)) {
      sample_i <- unique_samples[i]
      sample_i_strains <- sort(relabun_tab[which(relabun_tab$sample == sample_i), 'strain'])
      
      for (j in (i + 1):length(unique_samples)) {
        sample_j <- unique_samples[j]
        sample_j_strains <- sort(relabun_tab[which(relabun_tab$sample == sample_j), 'strain'])
        
        per_between_sample_divergence <- numeric()
        
        for (strain1 in sample_i_strains) {
          
          for (strain2 in sample_j_strains) {
            
            if (strain1 == strain2) { next }
            
            strain_combo <- paste(sort(c(strain1, strain2)), collapse = '_')
            
            per_between_sample_divergence <- c(per_between_sample_divergence,
                                               divergence_tab[strain_combo, 'percent_identity'])
          }
          
        }
        
        between_sample_divergence <- c(between_sample_divergence,
                                       mean(per_between_sample_divergence))
        sample1_set <- c(sample1_set, sample_i)
        sample2_set <- c(sample2_set, sample_j)
      }
    }
  }
  
  return(data.frame(sample1 = sample1_set,
                    sample2 = sample2_set,
                    mean_between_divergence = between_sample_divergence))
}

allele_divergence <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/strainfacts_allele_hamming.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
allele_divergence <- allele_divergence[which(allele_divergence$comparable_length >= 100), ]

allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

# Then compute for alleles.
# Run these in parallel as there are so many.
within_allele_div_raw_INPUT <- list()
between_allele_div_raw_INPUT <- list()

for (dataset in names(allele_relabun)) {
  
  for (sp in names(allele_relabun[[dataset]])) {
    
    for (gene in names(allele_relabun[[dataset]][[sp]])) {
      
      gene_dataset <- paste(gene, dataset, sep = '.')
      
      if (! gene_dataset %in% allele_divergence$annot) { next }
      
      allele_divergence_subset <- allele_divergence[which(allele_divergence$annot == gene_dataset), ]
      rownames(allele_divergence_subset) <- paste(allele_divergence_subset$seq1, allele_divergence_subset$seq2, sep = '_')
      
      allele_relabun_subset <- allele_relabun[[dataset]][[sp]][[gene]]
      allele_relabun_subset$strain <- paste('allele', as.character(allele_relabun_subset$strain), sep = '_')
      rownames(allele_relabun_subset) <- paste(dataset, rownames(allele_relabun_subset), sep = '.')
      
      within_allele_div_raw_INPUT[[gene_dataset]] <- list(divergence_tab = allele_divergence_subset,
                                                          relabun_tab = allele_relabun_subset,
                                                          species = sp,
                                                          dataset = dataset,
                                                          gene = gene,
                                                          gene_dataset = gene_dataset)
      
      between_allele_div_raw_INPUT[[gene_dataset]] <- list(divergence_tab = allele_divergence_subset,
                                                           relabun_tab = allele_relabun_subset,
                                                           species = sp,
                                                           dataset = dataset,
                                                           gene = gene,
                                                           gene_dataset = gene_dataset)
    }
    
  }
  
}

within_allele_div_raw <- parallel::mclapply(X = within_allele_div_raw_INPUT,
                                            function(x) {
                                              df_out <- compute_within_sample_mean_allele_divergence(divergence_tab = x$divergence_tab,
                                                                                              relabun_tab = x$relabun_tab)
                                              
                                              df_out$species <- x$species
                                              df_out$dataset <- x$dataset
                                              df_out$gene <- x$gene
                                              
                                              return(df_out)
                                            },
                                            mc.cores = 55)

between_allele_div_raw <- parallel::mclapply(X = between_allele_div_raw_INPUT,
                                             function(x) {
                                               df_out <- compute_between_sample_mean_allele_divergence(divergence_tab = x$divergence_tab,
                                                                                                relabun_tab = x$relabun_tab)
                                               
                                               df_out$species <- x$species
                                               df_out$dataset <- x$dataset
                                               df_out$gene <- x$gene
                                               
                                               return(df_out)
                                             },
                                             mc.cores = 55)

within_allele_div <- do.call(rbind, within_allele_div_raw)
rownames(within_allele_div) <- NULL
within_allele_div <- within_allele_div[which(! is.na(within_allele_div$mean_within_divergence)), ]

between_allele_div <- do.call(rbind, between_allele_div_raw)
rownames(between_allele_div) <- NULL
between_allele_div <- between_allele_div[which(! is.na(between_allele_div$mean_between_divergence)), ]

write.table(x = within_allele_div,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/within_sample_allele_divergence.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')

write.table(x = between_allele_div,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/between_sample_allele_divergence.tsv',
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')
