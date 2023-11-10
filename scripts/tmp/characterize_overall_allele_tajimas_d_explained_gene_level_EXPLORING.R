rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ComplexHeatmap)

glmm_gene_summary_out <- readRDS(file = '/data1/gdouglas/tmp/glmm_gene_summary_by_species.rds')

COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")
COG_category_info[which(COG_category_info$COG_category == ''), 'COG_category'] <- '-'

COG_map <- readRDS('/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds')
COG_descrip <- read.table('/data1/gdouglas/db/COG_definitions/cog-20.def.tab',
                          header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
COG_category_descrip <- read.table('/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv',
                                   header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

COG_category_descrip$clean <- paste(rownames(COG_category_descrip), COG_category_descrip$V2, sep = ' - ')
COG_category_descrip['-', 'clean'] <- 'No COG annotation'

sig_gene_tallies_down <- data.frame(matrix(0,
                                          nrow = nrow(COG_category_descrip),
                                          ncol = length(glmm_gene_summary_out)))
rownames(sig_gene_tallies_down) <- rownames(COG_category_descrip)
colnames(sig_gene_tallies_down) <- names(glmm_gene_summary_out)

sig_gene_tallies_up <- data.frame(matrix(0,
                                        nrow = nrow(COG_category_descrip),
                                        ncol = length(glmm_gene_summary_out)))
rownames(sig_gene_tallies_up) <- rownames(COG_category_descrip)
colnames(sig_gene_tallies_up) <- names(glmm_gene_summary_out)


all_tested_genes <- data.frame(matrix(0,
                                      nrow = nrow(COG_category_descrip),
                                      ncol = length(glmm_gene_summary_out)))
rownames(all_tested_genes) <- rownames(COG_category_descrip)
colnames(all_tested_genes) <- names(glmm_gene_summary_out)

for (sp in colnames(sig_gene_tallies_up)) {
  sp_coef <- glmm_gene_summary_out[[sp]]$coefficients$cond
  sig_genes_lower <- rownames(sp_coef[which(sp_coef[, 4] < 0.05 & sp_coef[, 1] < 0), ])
  sig_genes_lower <- gsub('^gene', '', sig_genes_lower[which(sig_genes_lower != '(Intercept)')])
  
  sig_genes_higher <- rownames(sp_coef[which(sp_coef[, 4] < 0.05 & sp_coef[, 1] > 0), ])
  sig_genes_higher <- gsub('^gene', '', sig_genes_higher[which(sig_genes_higher != '(Intercept)')])

  sig_genes_lower_COG_categories <- COG_category_info[sig_genes_lower, 'COG_category']
  if (length(which(is.na(sig_genes_lower_COG_categories))) > 0) {
    sig_genes_lower_COG_categories[which(is.na(sig_genes_lower_COG_categories))] <- '-'
  }
  
  sig_genes_higher_COG_categories <- COG_category_info[sig_genes_higher, 'COG_category']
  if (length(which(is.na(sig_genes_higher_COG_categories))) > 0) {
    sig_genes_higher_COG_categories[which(is.na(sig_genes_higher_COG_categories))] <- '-'
  }
  
  for (down_COG_category in sig_genes_lower_COG_categories) {
    categories <- strsplit(down_COG_category, ',')[[1]]
    count <- 1 / length(categories)
    for (tmp in categories) {
      sig_gene_tallies_down[tmp, sp] <- sig_gene_tallies_down[tmp, sp] + count
    }
  }
  
  for (up_COG_category in sig_genes_higher_COG_categories) {
    categories <- strsplit(up_COG_category, ',')[[1]]
    count <- 1 / length(categories)
    for (tmp in categories) {
      sig_gene_tallies_up[tmp, sp] <- sig_gene_tallies_up[tmp, sp] + count
    }
  }
  
  
  tested_genes <-  gsub('^gene', '', rownames(sp_coef)[which(rownames(sp_coef) != '(Intercept)')])
  tested_genes_COG_categories <- COG_category_info[tested_genes, 'COG_category']
  if (length(which(is.na(tested_genes_COG_categories))) > 0) {
    tested_genes_COG_categories[which(is.na(tested_genes_COG_categories))] <- '-'
  }
  
  for (tested_gene_category in tested_genes_COG_categories) {
    categories <- strsplit(tested_gene_category, ',')[[1]]
    count <- 1 / length(categories)
    for (tmp in categories) {
      all_tested_genes[tmp, sp] <- all_tested_genes[tmp, sp] + count
    }
  }
}

colnames(sig_gene_tallies_up) <- paste('up', colnames(sig_gene_tallies_up), sep = '.')
colnames(sig_gene_tallies_down) <- paste('down', colnames(sig_gene_tallies_down), sep = '.')

all_sig_gene_tallies <- cbind(sig_gene_tallies_down, sig_gene_tallies_up)
all_sig_gene_tallies <- all_sig_gene_tallies[-which(rowSums(all_sig_gene_tallies) == 0), ]
all_sig_gene_tallies <- round(all_sig_gene_tallies, 1)

all_sig_gene_percents <- sweep(x = all_sig_gene_tallies, MARGIN = 2, STATS = colSums(all_sig_gene_tallies), FUN = '/') * 100
all_sig_gene_percents[is.na(all_sig_gene_percents)] <- 0

column_names_prepped <- c(names(glmm_gene_summary_out), names(glmm_gene_summary_out))
Heatmap(matrix = as.matrix(all_sig_gene_percents),
        col = circlize::colorRamp2(c(0, 60), c('white', 'firebrick')),
        
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(title = 'Percent'),
        column_labels = gt_render(paste('*', column_names_prepped, '*', sep = '')),
        show_column_names = TRUE,
        column_names_rot = 45,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        column_split = c(rep(x = 'Significantly lower', 16),
                         rep(x = 'Significantly higher', 16)),
        column_gap = unit(10, "mm"),
        
        row_names_side = 'left',
       
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(all_sig_gene_tallies[i, j] > 0))
            grid.text(all_sig_gene_tallies[i, j], x, y, gp = gpar(fontsize = 10))
        },
        width = 1)






Heatmap(matrix = as.matrix(all_sig_gene_tallies),
        col = circlize::colorRamp2(c(0, 80), c('white', 'firebrick')),
        
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(title = 'Count'),
        column_labels = gt_render(paste('*', column_names_prepped, '*', sep = '')),
        show_column_names = TRUE,
        column_names_rot = 45,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        column_split = c(rep(x = 'Significantly lower', 16),
                         rep(x = 'Significantly higher', 16)),
        column_gap = unit(10, "mm"),
        
        row_names_side = 'left',
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(all_sig_gene_tallies[i, j] > 0))
            grid.text(all_sig_gene_tallies[i, j], x, y, gp = gpar(fontsize = 10))
        },
        width = 1)



sig_gene_tallies_up_prep <- sig_gene_tallies_up[-which(rowSums(sig_gene_tallies_up) == 0), ]
sig_gene_tallies_up_prep <- round(sig_gene_tallies_up_prep, 1)

sig_gene_tallies_up_prep <- sig_gene_tallies_up_prep[, names(sort(colSums(sig_gene_tallies_up_prep), decreasing = TRUE))]
column_names_prepped <- gsub('^up.', '', colnames(sig_gene_tallies_up_prep))
column_names_prepped <- gsub('_', ' ', column_names_prepped)
column_names_prepped <- gsub(' sp$', ' sp.', column_names_prepped)
column_names_prepped <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', column_names_prepped)

sig_gene_tallies_up_prep <- sig_gene_tallies_up_prep[-which(rownames(sig_gene_tallies_up_prep) == 'Z'), ]

COG_cat_order <- sort(rownames(sig_gene_tallies_up_prep))
COG_cat_order <- c(COG_cat_order[-which(COG_cat_order %in% c('-', 'S', 'R'))], 'R', 'S', '-')
sig_gene_tallies_up_prep <- sig_gene_tallies_up_prep[COG_cat_order, ]

rownames(sig_gene_tallies_up_prep) <- COG_category_descrip[rownames(sig_gene_tallies_up_prep), 'clean']

percent_sig_up_tested_genes <- round(100 * (colSums(sig_gene_tallies_up)/colSums(all_tested_genes)), 2)
percent_sig_up_tested_genes <- percent_sig_up_tested_genes[colnames(sig_gene_tallies_up_prep)]
percent_sig_up_tested_genes_annot <- HeatmapAnnotation(`Percent of tested genes` = anno_barplot(percent_sig_up_tested_genes,
                                                                                                add_numbers = TRUE, 
                                                                                                height = unit(2, "cm")))
upper_heatmap <- Heatmap(matrix = as.matrix(sig_gene_tallies_up_prep),
        col = circlize::colorRamp2(c(0, 30), c('white', 'firebrick')), 
        bottom_annotation = percent_sig_up_tested_genes_annot,
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(title = 'Count'),
        column_labels = gt_render(paste('*', column_names_prepped, '*', sep = '')),
        show_column_names = TRUE,
        column_names_rot = 45,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_gap = unit(10, "mm"),

        row_names_side = 'left',

        cell_fun = function(j, i, x, y, width, height, fill) {
          if(! is.na(sig_gene_tallies_up_prep[i, j] > 0))
            grid.text(sig_gene_tallies_up_prep[i, j], x, y, gp = gpar(fontsize = 10))
        },
        width = 1)

draw(upper_heatmap, padding = unit(c(2, 70, 2, 2), "mm"))
