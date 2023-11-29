rm(list = ls(all.names = TRUE))

# Per-species, get mean and SD of all pairwise Jaccard distances between sample pairs.
# Do this separately for all three datatypes: strains, accessory genes, and accessory gene alleles.
# Note that this should be based on presence/absence in all cases, to make them easier to compare.

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')
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

strainfacts_R2_tables <- R2_tables
strainfacts_P_tables <- P_tables

# StrainFacts
for (sp in names(strain_samples)) {

  for (d in names(strain_samples[[sp]])) {
    
    sample_subset <- strain_samples[[sp]][[d]]
    
    if (length(sample_subset) < 20) { next }
    
    tmp_binary <- strainfacts_relabun[[sp]][[d]][, sample_subset, drop = FALSE]
    tmp_binary[tmp_binary > 0] <- 1
    tmp_jaccard <- stats::dist(t(tmp_binary), method = "binary", diag = FALSE, upper = FALSE)
    
    adonis2_out <- run_adonis2_w_metadata(dist_in = tmp_jaccard, meta_in = dataset_meta[[d]])
    vars_present <- intersect(rownames(adonis2_out), colnames(R2_tables[[d]]))
    
    for (var_present in vars_present) {
      strainfacts_R2_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'R2']
      strainfacts_P_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'Pr(>F)']
    }
  }
}


# StrainGST
for (sp in names(strain_samples)) {
  
  for (d in names(strain_samples[[sp]])) {
    
    tmp_binary <- straingst_relabun[[sp]][[d]][, strain_samples[[sp]][[d]], drop = FALSE]
    tmp_binary[tmp_binary > 0] <- 1
    tmp_jaccard <- stats::dist(t(tmp_binary), method = "binary", diag = FALSE, upper = FALSE)
    all_jaccard_raw[[paste('straingst', sp, d, sep = '_')]] <- data.frame(species = sp,
                                                                          datatype = 'StrainGST',
                                                                          dataset = d,
                                                                          mean_jaccard = mean(tmp_jaccard),
                                                                          sd_jaccard = sd(tmp_jaccard))
  }
}

straingst_R2_tables <- R2_tables
straingst_P_tables <- P_tables

# straingst
for (sp in names(strain_samples)) {
  
  for (d in names(strain_samples[[sp]])) {
    
    sample_subset <- strain_samples[[sp]][[d]]
    
    if (length(sample_subset) < 20) { next }
    
    tmp_binary <- straingst_relabun[[sp]][[d]][, sample_subset, drop = FALSE]
    tmp_binary[tmp_binary > 0] <- 1
    tmp_jaccard <- stats::dist(t(tmp_binary), method = "binary", diag = FALSE, upper = FALSE)
    
    adonis2_out <- run_adonis2_w_metadata(dist_in = tmp_jaccard, meta_in = dataset_meta[[d]])
    vars_present <- intersect(rownames(adonis2_out), colnames(R2_tables[[d]]))
    
    for (var_present in vars_present) {
      straingst_R2_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'R2']
      straingst_P_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'Pr(>F)']
    }
  }
}


# For accessory genes.
gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1

gene_R2_tables <- R2_tables
gene_P_tables <- P_tables

for (sp in names(strain_samples)) {
  for (d in names(strain_samples[[sp]])) {
    sp_strain_dataset_samples <- strain_samples[[sp]][[d]]
    if (length(sp_strain_dataset_samples) == 0) { next }
    sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
    sp_dataset_gene_cooccur <- gene_cooccur[sp_accessory, sp_strain_dataset_samples, drop = FALSE]
    sp_dataset_gene_cooccur <- sp_dataset_gene_cooccur[which(rowSums(sp_dataset_gene_cooccur) > 0), ]
    tmp_gene_binary <- stats::dist(t(sp_dataset_gene_cooccur), method = "binary", diag = FALSE, upper = FALSE)

    adonis2_out <- run_adonis2_w_metadata(dist_in = tmp_gene_binary, meta_in = dataset_meta[[d]])
    vars_present <- intersect(rownames(adonis2_out), colnames(R2_tables[[d]]))
    
    for (var_present in vars_present) {
      gene_R2_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'R2']
      gene_P_tables[[d]][sp, var_present] <- adonis2_out[var_present, 'Pr(>F)']
    }

  }
}


# Plot as heatmap.
prep_datatype_heatmap_tabs <- function(strainfacts_P_tab,
                                       strainfacts_R2_tab,
                                       
                                       straingst_P_tab,
                                       straingst_R2_tab,
                                       
                                       gene_P_tab,
                                       gene_R2_tab) {
  
  # Make sure that all columns match.
  var_cols <- colnames(strainfacts_P_tab)
  
  if (! identical(colnames(strainfacts_R2_tab), var_cols)) {
    stop('Column mismatch')
  } else if (! identical(colnames(straingst_P_tab), var_cols)) {
    stop('Column mismatch')
  } else if (! identical(colnames(straingst_R2_tab), var_cols)) {
    stop('Column mismatch')
  } else if (! identical(colnames(gene_P_tab), var_cols)) {
    stop('Column mismatch')
  } else if (! identical(colnames(gene_R2_tab), var_cols)) {
    stop('Column mismatch')
  }

  colnames(strainfacts_R2_tab) <- paste('strainfacts', colnames(strainfacts_R2_tab))
  colnames(strainfacts_P_tab) <- paste('strainfacts', colnames(strainfacts_P_tab))

  colnames(straingst_R2_tab) <- paste('straingst', colnames(straingst_R2_tab))
  colnames(straingst_P_tab) <- paste('straingst', colnames(straingst_P_tab))

  colnames(gene_R2_tab) <- paste('gene', colnames(gene_R2_tab))
  colnames(gene_P_tab) <- paste('gene', colnames(gene_P_tab))

  combined_R2 <- cbind(strainfacts_R2_tab, straingst_R2_tab)
  combined_R2 <- cbind(combined_R2, gene_R2_tab)

  combined_P <- cbind(strainfacts_P_tab, straingst_P_tab)
  combined_P <- cbind(combined_P, gene_P_tab)
  
  rownames(combined_R2) <- cleanup_species_names(rownames(combined_R2))
  rownames(combined_P) <- cleanup_species_names(rownames(combined_P))

  combined_P_char <- combined_P
  combined_P_char[combined_P < 0.05] <- 'Sig.'
  combined_P_char[combined_P >= 0.05] <- 'Non-sig.'
  
  combined_R2 <- format(round(combined_R2, 2), nsmall = 2)
  combined_R2[combined_R2 == '  NA'] <- ''
  combined_R2[combined_R2 == 'NA'] <- ''

  column_splits_vec <- c(rep('StrainFacts', length(var_cols)),
                         rep('StrainGST', length(var_cols)),
                         rep('Gene', length(var_cols)))
   column_splits_vec <- factor(column_splits_vec, levels = c('StrainFacts', 'StrainGST', 'Gene'))
  
  return(list(P = combined_P_char,
              R2 = combined_R2,
              column_labels = rep(var_cols, 3),
              column_splits = column_splits_vec))
}

heatmap_input_obj <- list()
for (d in datasets_w_metadata) {
  heatmap_input_obj[[d]] <- prep_datatype_heatmap_tabs(strainfacts_P_tab = strainfacts_P_tables[[d]],
                                                       strainfacts_R2_tab = strainfacts_R2_tables[[d]],
                                                       
                                                       straingst_P_tab = straingst_P_tables[[d]],
                                                       straingst_R2_tab = straingst_R2_tables[[d]],
                                                       
                                                       gene_P_tab = gene_P_tables[[d]],
                                                       gene_R2_tab = gene_R2_tables[[d]])
}

library(ComplexHeatmap)

Sig_cols <- c('red', 'grey95')
names(Sig_cols) <- c('Sig.', 'Non-sig.')

# Ellegaard2019
Ellegaard2019_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(heatmap_input_obj$Ellegaard2019$P),
  
  heatmap_legend_param = list(title = ''),
  
  col = Sig_cols,
  
  na_col = 'grey70',
  
  show_heatmap_legend = TRUE,
  
  column_labels = heatmap_input_obj$Ellegaard2019$column_labels,
  column_split = heatmap_input_obj$Ellegaard2019$column_splits,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  
  row_labels = gt_render(paste('*', rownames(heatmap_input_obj$Ellegaard2019$P), '*', sep = '')),
  
  row_names_side = 'left',
  row_dend_side = 'right',
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(! is.na(heatmap_input_obj$Ellegaard2019$R2[i, j] > 0))
      grid.text(heatmap_input_obj$Ellegaard2019$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
  })

draw(Ellegaard2019_adonis_heatmap, column_title = "Ellegaard 2019")


# Ellegaard2020
Ellegaard2020_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(heatmap_input_obj$Ellegaard2020$P),
  
  heatmap_legend_param = list(title = ''),
  
  col = Sig_cols,
  
  na_col = 'grey70',
  
  show_heatmap_legend = TRUE,
  
  column_labels = heatmap_input_obj$Ellegaard2020$column_labels,
  column_split = heatmap_input_obj$Ellegaard2020$column_splits,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  
  row_labels = gt_render(paste('*', rownames(heatmap_input_obj$Ellegaard2020$P), '*', sep = '')),
  
  row_names_side = 'left',
  row_dend_side = 'right',
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(! is.na(heatmap_input_obj$Ellegaard2020$R2[i, j] > 0))
      grid.text(heatmap_input_obj$Ellegaard2020$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
  })

draw(Ellegaard2020_adonis_heatmap, column_title = "Ellegaard 2020")



# Zhang2022
Zhang2022_adonis_heatmap <- ComplexHeatmap::Heatmap(
  as.matrix(heatmap_input_obj$Zhang2022$P),
  
  heatmap_legend_param = list(title = ''),
  
  col = Sig_cols,
  
  na_col = 'grey70',
  
  show_heatmap_legend = TRUE,
  
  column_labels = heatmap_input_obj$Zhang2022$column_labels,
  column_split = heatmap_input_obj$Zhang2022$column_splits,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  
  row_labels = gt_render(paste('*', rownames(heatmap_input_obj$Zhang2022$P), '*', sep = '')),
  
  row_names_side = 'left',
  row_dend_side = 'right',
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(! is.na(heatmap_input_obj$Zhang2022$R2[i, j] > 0))
      grid.text(heatmap_input_obj$Zhang2022$R2[i, j], x, y, gp = gpar(fontsize = 10), just = 'centre')
  })

draw(Zhang2022_adonis_heatmap, column_title = "Zhang2022")
