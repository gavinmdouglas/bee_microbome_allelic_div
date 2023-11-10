rm(list = ls(all.names = TRUE))

jaccard_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv.gz',
                              header = TRUE, sep = '\t', stringsAsFactors = FALSE)

jaccard_summary_clean <- jaccard_summary
jaccard_summary$datatype <- tolower(jaccard_summary$datatype)

jaccard_summary <- jaccard_summary[which(jaccard_summary$species != 'All'), ]

# Remove rows with no variation.
jaccard_summary <- jaccard_summary[which(jaccard_summary$sd_jaccard > 0), ]

datatypes <- unique(jaccard_summary$datatype)

raw_cor <- list()

for (datatype1 in datatypes[1:(length(datatypes) - 1)]) {

  for (datatype2 in datatypes[(which(datatypes == datatype1) + 1):length(datatypes)]) {

    for (dataset in unique(jaccard_summary$dataset)) {

      subset1 <- jaccard_summary[which(jaccard_summary$datatype == datatype1 & jaccard_summary$dataset == dataset), , drop = FALSE]
      if (nrow(subset1) == 0) { next }
      rownames(subset1) <- subset1$species

      subset2 <- jaccard_summary[which(jaccard_summary$datatype == datatype2 & jaccard_summary$dataset == dataset), , drop = FALSE]
      if (nrow(subset2) == 0) { next }
      rownames(subset2) <- subset2$species
      
      intersecting_species <- intersect(rownames(subset1), rownames(subset2))
      
      if (length(intersecting_species) < 10) { next }
      
      cor_out <- cor.test(subset1[intersecting_species, 'mean_jaccard'], subset2[intersecting_species, 'mean_jaccard'], method = 'spearman')
      
      raw_cor[[paste(datatype1, datatype2, dataset, sep = ',')]] <- data.frame(datatype1=datatype1,
                                                                               datatype2=datatype2,
                                                                               dataset=dataset,
                                                                               spearman_rho=cor_out$estimate,
                                                                               spearman_p=cor_out$p.value)
    }

  }

}

jaccard_cor <- do.call(rbind, raw_cor)

write.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_cross.datatype_cor.tsv',
            x = jaccard_cor,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')


