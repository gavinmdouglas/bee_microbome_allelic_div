rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)

jaccard_cor <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_cross.datatype_cor.tsv.gz',
                          header = TRUE, sep = '\t', stringsAsFactors = FALSE)

rownames(jaccard_cor) <- paste(jaccard_cor$datatype1, jaccard_cor$datatype2, jaccard_cor$dataset, sep = '|')

jaccard_cor <- jaccard_cor[, -c(1:3)]
jaccard_cor_rounded <- format(round(jaccard_cor, 2), nsmall = 2)

jaccard_cor$kendall_rho[which(jaccard_cor$kendall_p >= 0.05)] <- NA


Heatmap(matrix = as.matrix(jaccard_cor[, 1, drop = FALSE]),
        
        
        col = circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'firebrick')),
        
        show_heatmap_legend = TRUE,

        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        row_names_side = 'left',
  
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(jaccard_cor_rounded[, 1, drop = FALSE][i, j] > 0))
            grid.text(jaccard_cor_rounded[, 1, drop = FALSE][i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
        })


# Run sanity checks on correlations
# (Note that rows with 0 SD not filtered out in these tests, to make sure that doesn't change things.)
jaccard_summary <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv.gz',
                                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

jaccard_summary_clean <- jaccard_summary
jaccard_summary$datatype <- tolower(jaccard_summary$datatype)

jaccard_summary <- jaccard_summary[which(jaccard_summary$species != 'All'), ]

datasets <- unique(jaccard_summary$dataset)
datatypes <- unique(jaccard_summary$datatype)

# Zhang 2022 strainfacts vs allele
for (d1 in datatypes) {

  for(d2 in datatypes) {
    
    if (d1 == d2) { next }
    
    for (dataset in datasets) {
      
      tmp1 <- jaccard_summary[which(jaccard_summary$dataset == dataset & jaccard_summary$datatype == d1), ]
      tmp2 <- jaccard_summary[which(jaccard_summary$dataset == dataset & jaccard_summary$datatype == d2), ]
      
      rownames(tmp1) <- tmp1$species
      rownames(tmp2) <- tmp2$species
      
      intersecting_species <- intersect(tmp1$species, tmp2$species)
      
      tmp1 <- tmp1[intersecting_species, ]
      tmp2 <- tmp2[intersecting_species, ]
      
      print(paste(d1, d2, dataset))
      print(cor.test(tmp1$mean_jaccard, tmp2$mean_jaccard, method = 'kendall'))
      
    }
  }
}
    
