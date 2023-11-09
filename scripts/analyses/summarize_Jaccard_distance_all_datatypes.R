rm(list = ls(all.names = TRUE))

# Per-species, get mean and SD of all pairwise Jaccard distances across samples.
# Do this separately for all three datatypes: species, strains, accessory genes, and accessory gene alleles.
# (Although note that for species that this is based on all species of course).
# Note that this should be based on presence/absence in all cases, to make them easier to compare.

# For accessory gene alleles, also parse out the mean and sd of the numbers of alleles per sample
# (summarized as the mean of all genes per species).

# Calculate all these measures for each dataset separately.
datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')

sample_ids <- list()
for (d in datasets) {
  sample_ids[[d]] <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/SRRs/', d, '_SRRs.txt.gz', sep = ''),
                                header = FALSE, stringsAsFactors = FALSE)$V1
  
}

all_jaccard_raw <- list()

# For species.
species_90per <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz',
                            row.names = 1, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

for (d in datasets) {
  species_90per_dataset <- species_90per[intersect(rownames(species_90per), sample_ids[[d]]), ]
  species_90per_dataset_jaccard <- c(stats::dist(t(species_90per_dataset), method = "binary", diag = FALSE, upper = FALSE))
  all_jaccard_raw[[paste('species', d, sep = '_')]] <- data.frame(species = 'All',
                                                                  datatype = 'Species',
                                                                  dataset = d,
                                                                  mean_jaccard = mean(species_90per_dataset_jaccard),
                                                                  sd_jaccard = sd(species_90per_dataset_jaccard))
}


# For strains.
# Only consider samples intersecting between StrainFacts and StrainGST.
# Keep track of this for subsequent allele-based analyses.
strainfacts_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
straingst_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')

strain_species <- intersect(names(strainfacts_relabun), names(straingst_relabun))
strain_samples <- list()

for (sp in strain_species) {
  tmp_datasets <- intersect(names(strainfacts_relabun[[sp]]), names(straingst_relabun[[sp]]))
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

# StrainFacts
for (sp in names(strain_samples)) {

  for (d in names(strain_samples[[sp]])) {
    
    tmp_binary <- strainfacts_relabun[[sp]][[d]][, strain_samples[[sp]][[d]], drop = FALSE]
    tmp_binary[tmp_binary > 0] <- 1
    tmp_jaccard <- stats::dist(t(tmp_binary), method = "binary", diag = FALSE, upper = FALSE)
    all_jaccard_raw[[paste('strainfacts', sp, d, sep = '_')]] <- data.frame(species = sp,
                                                                      datatype = 'StrainFacts',
                                                                      dataset = d,
                                                                      mean_jaccard = mean(tmp_jaccard),
                                                                      sd_jaccard = sd(tmp_jaccard))
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
# For accessory genes.
gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1

for (sp in names(strain_samples)) {
  for (d in names(strain_samples[[sp]])) {
    sp_strain_dataset_samples <- strain_samples[[sp]][[d]]
    if (length(sp_strain_dataset_samples) == 0) { next }
    sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
    sp_dataset_gene_cooccur <- gene_cooccur[sp_accessory, sp_strain_dataset_samples, drop = FALSE]
    sp_dataset_gene_cooccur <- sp_dataset_gene_cooccur[which(rowSums(sp_dataset_gene_cooccur) > 0), ]
    tmp_gene_binary <- stats::dist(t(sp_dataset_gene_cooccur), method = "binary", diag = FALSE, upper = FALSE)
    all_jaccard_raw[[paste('gene', sp, d, sep = '_')]] <- data.frame(species = sp,
                                                                     datatype = 'Gene',
                                                                     dataset = d,
                                                                     mean_jaccard = mean(tmp_gene_binary),
                                                                     sd_jaccard = sd(tmp_gene_binary))
  }
}


# For accessory gene alleles also get the mean number of alleles per sample per gene/dataset.
allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

allele_count_raw <- list()

for (d in datasets) {
 
  for (sp in names(allele_relabun[[d]])) {
   
    sp_strain_dataset_samples <- strain_samples[[sp]][[d]]
    
    if (length(sp_strain_dataset_samples) == 0) { next }
    
    raw_out <- parallel::mclapply(X = names(allele_relabun[[d]][[sp]]), mc.cores = 40,
                                  FUN = function(gene) {
                                    
                                    sample_subset <- intersect(allele_relabun[[d]][[sp]][[gene]]$sample, sp_strain_dataset_samples)
                                    if (length(sample_subset) <= 2) { return(NA) }
                                    
                                    allele_binary <- allele_relabun[[d]][[sp]][[gene]][which(allele_relabun[[d]][[sp]][[gene]]$sample %in% sample_subset), ]
                                    if (length(unique(allele_binary$strain)) == 1) { return(NA) } # Skip if there is only one allele.
                                    
                                    allele_binary[which(allele_binary$community > 0), 'community'] <- 1
                                    allele_binary_wide <- reshape2::dcast(data = allele_binary, formula = sample ~ strain, value.var = 'community', fill = 0, drop = FALSE)[, -1, drop = FALSE]
                                    allele_binary_dist <- c(stats::dist(allele_binary_wide, method = "binary", diag = FALSE, upper = FALSE))
                                    alleles_per_sample <- rowSums(allele_binary_wide)
                                    
                                    return(list(mean_dist = mean(allele_binary_dist),
                                                mean_per_sample = mean(alleles_per_sample)))
                                    
                                  })
    
    all_mean_jaccard_dist <- sapply(1:length(raw_out), function(i) { if(identical(NA, raw_out[[i]])) { return(NA) } else { return(raw_out[[i]]$mean_dist) }})
    
    # Skip if there are fewer than 50 genes for which Jaccard distances based on alleles could be computed.
    nonNA_count <- length(which(! is.na(all_mean_jaccard_dist)))
    if (nonNA_count < 50) { next }
    
    all_jaccard_raw[[paste('allele', sp, d, sep = '_')]] <- data.frame(species = sp,
                                                                       datatype = 'Allele',
                                                                       dataset = d,
                                                                       mean_jaccard = mean(all_mean_jaccard_dist, na.rm = TRUE),
                                                                       sd_jaccard = sd(all_mean_jaccard_dist, na.rm = TRUE))
    
    all_mean_per_sample <- sapply(1:length(raw_out), function(i) { if(identical(NA, raw_out[[i]])) { return(NA) } else { return(raw_out[[i]]$mean_per_sample) }})

    allele_count_raw[[paste('allele', sp, d, sep = '_')]] <- data.frame(species = sp,
                                                                        dataset = d,
                                                                        num_genes_considered = length(which(! is.na(all_mean_per_sample))),
                                                                        mean_per_sample = mean(all_mean_per_sample, na.rm = TRUE),
                                                                        sd_per_sample = sd(all_mean_per_sample, na.rm = TRUE))
     
  }
  
}

all_jaccard <- do.call(rbind, all_jaccard_raw)
allele_count <- do.call(rbind, allele_count_raw)

rownames(all_jaccard) <- NULL
rownames(allele_count) <- NULL

write.table(x = all_jaccard,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(x = allele_count,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/accessory_allele_per_sample_summary.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
