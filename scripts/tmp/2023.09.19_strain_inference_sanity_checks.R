rm(list = ls(all.names = TRUE))

strain_abun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

min_depth <- 20

species <- names(strain_abun)[3]

sample_filename <- paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/prepped_core_input/',
                           species,
                           '/samples.tsv.gz', sep = '')
  
  sample_in <- read.table(sample_filename,
                          header = FALSE,
                          sep = '\t',
                          stringsAsFactors = FALSE,
                          row.names = 1)
  
  allele_in_file <- '/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/prepped_accessory_input/metagenotypes/Bartonella_apis_glpF_metagenotype.tsv.gz'
  
  allele_in <- read.table(allele_in_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  allele_in_ref <- allele_in[which(allele_in$allele == 'ref'), ]
  allele_in_alt <- allele_in[which(allele_in$allele == 'alt'), ]
  rownames(allele_in_ref) <- NULL
  rownames(allele_in_alt) <- NULL
  
  if (! identical(allele_in_ref[, c('sample', 'position')],
                  allele_in_alt[, c('sample', 'position')])) {
    stop('Error - ref and alt rows in different order.')
    
  }
  
  allele_combined <- data.frame(position = allele_in_ref$position,
                                ref_depth = allele_in_ref$metagenotype,
                                alt_depth = allele_in_alt$metagenotype)
  allele_combined$total_depth <- allele_combined$ref_depth + allele_combined$alt_depth
  allele_combined$alt_freq <- allele_combined$alt_depth / allele_combined$total_depth
  
  metagenotype <- data.frame(sample = sample_in[as.character(metagenotype_ref$sample), 'V2'],
                             position = metagenotype_ref$position,
                             ref_depth = metagenotype_ref$metagenotype,
                             alt_depth = metagenotype_alt$metagenotype)
  
  metagenotype_filename <- paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/prepped_core_input/',
                                 species,
                                 '/metagenotypes.tsv.gz', sep = '')
  
  metagenotype_in <- read.table(metagenotype_filename,
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
  metagenotype_ref <- metagenotype_in[which(metagenotype_in$allele == 'ref'), ]
  metagenotype_alt <- metagenotype_in[which(metagenotype_in$allele == 'alt'), ]
  rownames(metagenotype_ref) <- NULL
  rownames(metagenotype_alt) <- NULL
  
  if (! identical(metagenotype_ref[, c('sample', 'position')],
                  metagenotype_alt[, c('sample', 'position')])) {
    stop('Error - ref and alt rows in different order.')
    
  }
  
  metagenotype <- data.frame(sample = sample_in[as.character(metagenotype_ref$sample), 'V2'],
                             position = metagenotype_ref$position,
                             ref_depth = metagenotype_ref$metagenotype,
                             alt_depth = metagenotype_alt$metagenotype)
  
  metagenotype$total_depth <- metagenotype$ref_depth + metagenotype$alt_depth
  metagenotype$alt_freq <- metagenotype$alt_depth / metagenotype$total_depth
  
  cat('::: {.panel-tabset}\n\n')
  num_to_sample <- min(10, nrow(sample_in))
  for (sample_name in sample(sample_in$V2, num_to_sample)) {
    
    strainfacts_inferred <- length(which(strain_abun[[species]][, sample_name] > 0))
    metagenotype_subset <- metagenotype[which(metagenotype$sample == sample_name), ]
    metagenotype_subset <- metagenotype_subset[which(metagenotype_subset$alt_freq > 0), ]
    metagenotype_subset <- metagenotype_subset[which(metagenotype_subset$total_depth >= min_depth), ]
    
    if (nrow(metagenotype_subset) < 10) { next }
    
    cat('\n## ', sample_name, '\n')
    
    panel_title <- paste(sample_name, ' - ',
                         as.character(strainfacts_inferred),
                         ' inferred by StrainFacts',
                         sep = '')
    hist(metagenotype_subset$alt_freq,
         main = panel_title,
         xlab = 'Alt. freq.',
         breaks = 50,
         xlim = c(0, 1))
    cat('\n\n')
    
  }
  cat(':::\n\n')
  
}

