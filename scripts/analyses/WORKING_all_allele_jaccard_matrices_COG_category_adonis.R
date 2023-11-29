rm(list = ls(all.names = TRUE))

# Per-species, get mean and SD of all pairwise Jaccard distances between sample pairs.
# Do this separately for all three datatypes: strains, accessory genes, and accessory gene alleles.
# Note that this should be based on presence/absence in all cases, to make them easier to compare.

library(vegan)

Ellegaard2019_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Ellegaard2019_metadata.tsv.gz',
                                     header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

Ellegaard2020_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Ellegaard2020_metadata.tsv.gz',
                                     header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

Zhang_metadata <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/metadata/Zhang2022_metadata.tsv.gz',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
Zhang_metadata$Day <- gsub('Day ', '', Zhang_metadata$Day)
Zhang_metadata$Day <- factor(Zhang_metadata$Day, levels = c('7', '11', '19'))
Zhang_metadata$Treatment <- factor(Zhang_metadata$Treatment,
                                   levels = c('Control', 'Antibiotic'))

dataset_meta <- list(Ellegaard2019 = Ellegaard2019_metadata[, -c(1, 2)],
                     Ellegaard2020 = Ellegaard2020_metadata[, -c(1, 2)],
                     Zhang2022 = Zhang_metadata[, -1])

run_adonis2_w_metadata <- function(dist_in, meta_in) {
  
  dist_samples <- rownames(as.matrix(dist_in))
  
  if (length(setdiff(dist_samples, rownames(meta_in))) > 0) { stop('Sample mismatch!') }
  
  meta_in <- meta_in[dist_samples, , drop = FALSE]
  
  vars <- character()
  for (meta_colname in colnames(meta_in)) {
    meta_breakdown <- table(meta_in[, meta_colname])
    meta_breakdown <- meta_breakdown[meta_breakdown > 5]
    if (length(meta_breakdown) > 1) {
      vars <- c(vars, meta_colname)
    }
  }
  
  if (length(vars) > 0) {
    adonis_formula <- as.formula(paste('dist_in ~ ', paste(vars, collapse = ' + '), sep = ''))
    return(adonis2(formula = adonis_formula, data = meta_in))
  } else {
    return(NULL)
  }
}

# Calculate all these measures for each dataset separately.
datasets_w_metadata <- c('Ellegaard2019', 'Ellegaard2020', 'Zhang2022')

sample_ids <- list()
for (d in datasets_w_metadata) {
  sample_ids[[d]] <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/', d, '_SRRs.txt.gz', sep = ''),
                                header = FALSE, stringsAsFactors = FALSE)$V1
  
}

all_jaccard_raw <- list()

# For strains.
# Only consider samples intersecting between StrainFacts and StrainGST.
# Keep track of this for subsequent allele-based analyses.
strainfacts_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
straingst_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')

strain_species <- intersect(names(strainfacts_relabun), names(straingst_relabun))
strain_samples <- list()

for (sp in strain_species) {
  tmp_datasets <- intersect(names(strainfacts_relabun[[sp]]), names(straingst_relabun[[sp]]))
  tmp_datasets <- intersect(datasets_w_metadata, tmp_datasets)
  if (length(tmp_datasets) == 0) { next }
  
  for (tmp_d in tmp_datasets) {
    tmp_intersecting_samp <- intersect(colnames(strainfacts_relabun[[sp]][[tmp_d]]), colnames(straingst_relabun[[sp]][[tmp_d]]))
    if (length(tmp_intersecting_samp) == 0) { next }
    
    if (! sp %in% names(strain_samples))   {
      strain_samples[[sp]] <- list()
    }
    
    strain_samples[[sp]][[tmp_d]] <- tmp_intersecting_samp
  }
}

# Prep Adonis2 tables.
R2_tables <- list()
for (d in datasets_w_metadata) {
  R2_tables[[d]] <- data.frame(matrix(NA,
                                      nrow = length(straingst_relabun),
                                      ncol = ncol(dataset_meta[[d]])))
  rownames(R2_tables[[d]]) <- names(straingst_relabun)
  colnames(R2_tables[[d]]) <- colnames(dataset_meta[[d]])
}
P_tables <- R2_tables

COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")
COG_category_info[which(COG_category_info$COG_category == ''), 'COG_category'] <- '-'

COG_descrip <- read.table('/data1/gdouglas/db/COG_definitions/cog-20.def.tab',
                          header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
COG_category_descrip <- read.table('/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv',
                                   header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
categories_to_ignore <- c('A', 'B', 'Y', 'Z')

all_allele_jaccard <- readRDS(file = '/home/gdouglas/tmp/per_allele_jaccard_for_metadata_datasets.rds')

COG_categories <- setdiff(c(rownames(COG_category_descrip), '-'), categories_to_ignore)

R2_raw <- list()
P_raw <- list()

for (d in datasets_w_metadata) {
  print(d)
  
  R2_raw[[d]] <- list()
  P_raw[[d]] <- list()
  
  for (sp in names(all_allele_jaccard[[d]])) {
    print(sp)
    sp_COG_annot <- COG_category_info[intersect(names(all_allele_jaccard[[d]][[sp]]), rownames(COG_category_info)), ]
    
    missing_genes <- setdiff(names(all_allele_jaccard[[d]][[sp]]), rownames(COG_category_info))
    if (length(missing_genes) > 0) {
      sp_COG_annot[missing_genes, 'COG_category'] <- '-'
    }

    for (COG_category in COG_categories) {
      print(COG_category)
      COG_gene_hits <- rownames(sp_COG_annot)[grep(COG_category, sp_COG_annot$COG_category)]
      
      if (length(COG_gene_hits) < 50) { next }
      
      sample_dist_sum <- data.frame(matrix(0, nrow = nrow(dataset_meta[[d]]), ncol = nrow(dataset_meta[[d]])))
      rownames(sample_dist_sum) <- rownames(dataset_meta[[d]])
      colnames(sample_dist_sum) <- rownames(dataset_meta[[d]])
      
      sample_tally <- sample_dist_sum
      
      COG_category_samples <- character()
      
      for (g in COG_gene_hits) {
       
        gene_allele_jaccard <- all_allele_jaccard[[d]][[sp]][[g]]
        
        if (class(gene_allele_jaccard) == 'logical') {
          next
        } else {
          gene_allele_jaccard <- as.matrix(gene_allele_jaccard)
          sample_dist_sum[rownames(gene_allele_jaccard), colnames(gene_allele_jaccard)] <- sample_dist_sum[rownames(gene_allele_jaccard), colnames(gene_allele_jaccard)] + gene_allele_jaccard
          sample_tally[rownames(gene_allele_jaccard), colnames(gene_allele_jaccard)] <- sample_tally[rownames(gene_allele_jaccard), colnames(gene_allele_jaccard)] + 1
          COG_category_samples <- unique(c(COG_category_samples, rownames(gene_allele_jaccard)))
        }

      }
      
      sample_dist_sum <- sample_dist_sum[COG_category_samples, COG_category_samples]
      sample_tally <- sample_tally[COG_category_samples, COG_category_samples]
      
      # Ignore sample pairs with < 50 genes.
      sample_dist_sum[sample_tally < 50] <- NA
      
      if (sum(colSums(! is.na(sample_dist_sum))) < 20) { next }
      
      # Remove samples with most to fewest NAs, until no NAs are left.
      NAs_per_samples <- colSums(is.na(sample_dist_sum))

      while(max(NAs_per_samples) != 0) {

        sample_i_to_rm <- which.max(NAs_per_samples)
        sample_dist_sum <- sample_dist_sum[-sample_i_to_rm, -sample_i_to_rm, drop = FALSE]
        sample_tally <- sample_tally[-sample_i_to_rm, -sample_i_to_rm, drop = FALSE]

        NAs_per_samples <- colSums(is.na(sample_dist_sum))

      }

      mean_dist <- as.dist(sample_dist_sum / sample_tally)

      adonis2_out <- run_adonis2_w_metadata(dist_in = mean_dist, meta_in = dataset_meta[[d]])
      vars_present <- intersect(rownames(adonis2_out), colnames(R2_tables[[d]]))

      if (is.null(vars_present)) { next }
      
      test_id <- paste(d, sp, COG_category)
      
      dummy_df <- data.frame(dataset = d, species = sp, COG_category = COG_category)

      for (possible_col in colnames(R2_tables[[d]])) {
        dummy_df[, possible_col] <- NA
      }
      
      R2_raw[[d]][[test_id]] <- dummy_df
      P_raw[[d]][[test_id]] <- dummy_df

      for (var_present in vars_present) {
        R2_raw[[d]][[test_id]][, var_present] <- adonis2_out[var_present, 'R2']
        P_raw[[d]][[test_id]][, var_present] <- adonis2_out[var_present, 'Pr(>F)']
      }
    }
  }
}

combined_R2 <- list()
combined_P <- list()

for (d in datasets_w_metadata) {
  combined_R2[[d]] <- do.call(rbind, R2_raw[[d]])
  combined_P[[d]] <- do.call(rbind, P_raw[[d]])
}

prep_heatmap_tabs <- function(P_tab,
                              R2_tab,
                              var_cols = c('Day', 'Treatment'),
                              prefix_to_rm = '') {
  P_tab <- P_tab[, var_cols]
  P_tab_char <- P_tab[, var_cols]
  R2_tab <- format(round(R2_tab[, var_cols], 2), nsmall = 2)
  
  P_tab_char[P_tab < 0.05] <- 'Sig.'
  P_tab_char[P_tab >= 0.05] <- 'Non-sig.'

  rownames(P_tab_char) <- gsub(prefix_to_rm, '', rownames(P_tab_char))
  rownames(R2_tab) <- gsub(prefix_to_rm, '', rownames(R2_tab))

  return(list(P = P_tab_char, R2 = R2_tab))
  
}

# Ellegaard2019
Ellegaard2019_heatmap_input <- prep_heatmap_tabs(P_tab = combined_P$Ellegaard2019,
                                                 R2_tab = combined_R2$Ellegaard2019,
                                                 var_cols = colnames(dataset_meta$Ellegaard2019),
                                                 prefix_to_rm <- 'Ellegaard2019 ')

Sig_cols <- c('red', 'grey95')
names(Sig_cols) <- c('Sig.', 'Non-sig.')

Ellegaard2019_cluster <- hclust(dist(Ellegaard2019_heatmap_input$R2))

Ellegaard2019_allele_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(Ellegaard2019_heatmap_input$P),

                        heatmap_legend_param = list(title = ''),
                        
                        col = Sig_cols,
                        
                        na_col = 'grey70',
                        
                        show_heatmap_legend = TRUE,
                        
                        cluster_rows = Ellegaard2019_cluster,
                        cluster_columns = FALSE,
                        cluster_column_slices = FALSE,
                        cluster_row_slices = FALSE,
                        
                        row_names_side = 'left',
                        row_dend_side = 'right',
                        column_names_rot = 45,
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(! is.na(Ellegaard2019_heatmap_input$R2[i, j] > 0))
                            grid.text(Ellegaard2019_heatmap_input$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
                        })

draw(Ellegaard2019_allele_adonis_heatmap, padding = unit(c(2, 20, 2, 2), "mm"), column_title = "Ellegaard 2019 (Allele-based Adonis)")

# Ellegaard2020
Ellegaard2020_heatmap_input <- prep_heatmap_tabs(P_tab = combined_P$Ellegaard2020,
                                                 R2_tab = combined_R2$Ellegaard2020,
                                                 var_cols = colnames(dataset_meta$Ellegaard2020),
                                                 prefix_to_rm <- 'Ellegaard2020 ')

Sig_cols <- c('red', 'grey95')
names(Sig_cols) <- c('Sig.', 'Non-sig.')

Ellegaard2020_cluster <- hclust(dist(Ellegaard2020_heatmap_input$R2))

Ellegaard2020_allele_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(Ellegaard2020_heatmap_input$P),
  
  heatmap_legend_param = list(title = ''),
  
  col = Sig_cols,
  
  na_col = 'grey70',
  
  show_heatmap_legend = TRUE,
  
  cluster_rows = Ellegaard2020_cluster,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  
  row_names_side = 'left',
  row_dend_side = 'right',
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(! is.na(Ellegaard2020_heatmap_input$R2[i, j] > 0))
      grid.text(Ellegaard2020_heatmap_input$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
  })

draw(Ellegaard2020_allele_adonis_heatmap, padding = unit(c(2, 20, 2, 2), "mm"), column_title = "Ellegaard 2020 (Allele-based Adonis)")


# Zhang2022
Zhang2022_heatmap_input <- prep_heatmap_tabs(P_tab = combined_P$Zhang2022,
                                                 R2_tab = combined_R2$Zhang2022,
                                                 var_cols = colnames(dataset_meta$Zhang2022),
                                                 prefix_to_rm <- 'Zhang2022 ')

Sig_cols <- c('red', 'grey95')
names(Sig_cols) <- c('Sig.', 'Non-sig.')

Zhang2022_allele_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(Zhang2022_heatmap_input$P),
  
  heatmap_legend_param = list(title = ''),
  
  col = Sig_cols,
  
  na_col = 'grey70',
  
  show_heatmap_legend = TRUE,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  
  row_names_side = 'left',
  row_dend_side = 'right',
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(! is.na(Zhang2022_heatmap_input$R2[i, j] > 0))
      grid.text(Zhang2022_heatmap_input$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
  })

draw(Zhang2022_allele_adonis_heatmap, padding = unit(c(2, 20, 2, 2), "mm"), column_title = "Zhang 2022 (Allele-based Adonis)")
