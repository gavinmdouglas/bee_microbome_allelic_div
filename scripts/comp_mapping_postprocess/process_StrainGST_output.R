# Clean-up StrainGST output to get per-species relative abundance matrices of strains across samples.
# Split this up by dataset.

rm(list = ls(all.names = TRUE))

library(reshape2)

datasets = c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

samples2dataset_RAW <- list()

for (d in datasets) {
  dataset_SRR_file <- paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/', d, '_SRRs.txt.gz', sep = '')
  samples2dataset_RAW[[d]] <- data.frame(SRR = read.table(dataset_SRR_file, stringsAsFactors = FALSE)$V1,
                                         dataset = d)
}
samples2dataset <- do.call(rbind, samples2dataset_RAW)
rownames(samples2dataset) <- samples2dataset$SRR

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz', stringsAsFactors = FALSE)$V1

straingst_outfiles <- list.files(path = '/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/strainge/straingst_out',
                                 pattern = '.strains.tsv', full.names = TRUE)

straingst <- list()

for (sp in all_species) {

  sp_outfiles <- grep(sp, straingst_outfiles, value = TRUE)
  
  if (length(sp_outfiles) == 0) { next }
  
  sp_samples <- sapply(sp_outfiles, basename)
  sp_samples <- gsub('_.*$', '', sp_samples)
  names(sp_samples) <- NULL
  
  sp_datasets <- samples2dataset[sp_samples, 'dataset']
  
  for (d in unique(sp_datasets)) {
    
    sp_d_i <- which(sp_datasets == d)
    
    tmp <- list()
    for (i in sp_d_i) {
        sp_samp <- sp_samples[i]
        sp_samp_file <- sp_outfiles[i]
        sp_samp_straingst <- read.table(sp_samp_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
        if (nrow(sp_samp_straingst) == 0) { next }
        tmp[[sp_samp]] <- data.frame(sample = sp_samp,
                                     strain = sp_samp_straingst$strain,
                                     relabun = sp_samp_straingst$rapct)
        
        # Normalize relabun to sum to 100%.
        tmp[[sp_samp]]$relabun <- (tmp[[sp_samp]]$relabun / sum(tmp[[sp_samp]]$relabun)) * 100

    }
    
    if (length(tmp) == 0) { next }
    
    sp_d_long <- do.call(rbind, tmp)
    
    sp_d_wide <- data.frame(dcast(data = sp_d_long, formula = strain ~ sample, value.var = "relabun", fill = 0))
    
    rownames(sp_d_wide) <- sp_d_wide$strain
    sp_d_wide <- sp_d_wide[, which(colnames(sp_d_wide) != 'strain'), drop = FALSE]
    
    if (! sp %in% names(straingst)) {
      straingst[[sp]] <- list()
    }
    
    straingst[[sp]][[d]] <- sp_d_wide
    
  }
  
}

saveRDS(object = straingst,
        file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')
