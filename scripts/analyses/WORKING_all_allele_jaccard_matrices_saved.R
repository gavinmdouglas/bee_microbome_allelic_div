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

allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

all_allele_jaccard <- list()

for (d in datasets_w_metadata) {
  
  print(d)
  
  all_allele_jaccard[[d]] <- list()
  
  for (sp in names(allele_relabun[[d]])) {
    
    print(sp)
    
    sp_strain_dataset_samples <- strain_samples[[sp]][[d]]
    
    if (length(sp_strain_dataset_samples) == 0) { next }
    
    raw_out <- parallel::mclapply(X = names(allele_relabun[[d]][[sp]]), mc.cores = 40,
                                  FUN = function(gene) {
                                    
                                    sample_subset <- intersect(allele_relabun[[d]][[sp]][[gene]]$sample, sp_strain_dataset_samples)
                                    if (length(sample_subset) <= 2) { return(NA) }
                                    
                                    allele_binary <- allele_relabun[[d]][[sp]][[gene]][which(allele_relabun[[d]][[sp]][[gene]]$sample %in% sample_subset), ]
                                    if (length(unique(allele_binary$strain)) == 1) { return(NA) } # Skip if there is only one allele.
                                    
                                    allele_binary[which(allele_binary$community > 0), 'community'] <- 1
                                    allele_binary_wide <- reshape2::dcast(data = allele_binary, formula = sample ~ strain, value.var = 'community', fill = 0, drop = FALSE)
                                    rownames(allele_binary_wide) <- allele_binary_wide$sample
                                    allele_binary_wide <- allele_binary_wide[, -1, drop = FALSE]
                                    return(stats::dist(allele_binary_wide, method = "binary", diag = FALSE, upper = FALSE))
                                  })
    
    names(raw_out) <- names(allele_relabun[[d]][[sp]])
    
    all_allele_jaccard[[d]][[sp]] <- raw_out

  }
  
}

saveRDS(object = all_allele_jaccard, file = '/home/gdouglas/tmp/per_allele_jaccard_for_metadata_datasets.rds')
