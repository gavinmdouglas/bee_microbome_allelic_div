rm(list = ls(all.names = TRUE))

# Per-species, get mean and SD of all pairwise Jaccard distances across samples.
# Do this separately for all three datatypes: species, strains, accessory genes, and accessory gene alleles.
# (Although note that for species that this is based on all species of course).
# Note that this should be based on presence/absence in all cases, to make them easier to compare.

# For accessory gene alleles, also parse out the mean and sd of the numbers of alleles per sample
# (summarized as the mean of all genes per species).


# For species.
species_90per <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/species_presence_core_90percent.tsv.gz',
                            row.names = 1, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
species_90per_jaccard <- c(stats::dist(t(species_90per), method = "binary", diag = FALSE, upper = FALSE))
species_jaccard_summary <- data.frame(species = 'All',
                                      datatype = 'Species',
                                      mean_jaccard = mean(species_90per_jaccard),
                                      sd_jaccard = sd(species_90per_jaccard))


# For strains.
strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
strain_binary <- lapply(strain_relabun, function(x) { x[x > 0] <- 1; return(x) })
strain_binary_dist <- lapply(strain_binary, function(x) { c(stats::dist(t(x), method = "binary", diag = FALSE, upper = FALSE)) })
strain_jaccard_summary <- data.frame(species = names(strain_binary_dist),
                                     datatype = 'Strain',
                                     mean_jaccard = NA,
                                     sd_jaccard = NA)
rownames(strain_jaccard_summary) <- names(strain_binary_dist)
for (sp in names(strain_binary_dist)) {
  strain_jaccard_summary[sp, 'mean_jaccard'] <- mean(strain_binary_dist[[sp]])
  strain_jaccard_summary[sp, 'sd_jaccard'] <- sd(strain_binary_dist[[sp]])
}

rownames(strain_jaccard_summary) <- NULL

# For accessory genes.
gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1
gene_jaccard_summary <- data.frame(species = names(strain_binary_dist),
                                   datatype = 'Gene',
                                   mean_jaccard = NA,
                                   sd_jaccard = NA)
rownames(gene_jaccard_summary) <- names(strain_binary_dist)
for (sp in names(strain_binary_dist)) {
  sp_strain_samples <- colnames(strain_binary[[sp]])
  sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
  sp_accessory_called <- intersect(sp_accessory, rownames(gene_cooccur))
  sp_gene_presence <- gene_cooccur[sp_accessory_called, sp_strain_samples]
  sp_gene_presence_dist <-  c(stats::dist(t(sp_gene_presence), method = "binary", diag = FALSE, upper = FALSE))
  gene_jaccard_summary[sp, 'mean_jaccard'] <- mean(sp_gene_presence_dist)
  gene_jaccard_summary[sp, 'sd_jaccard'] <- sd(sp_gene_presence_dist)
}

rownames(gene_jaccard_summary) <- NULL

# For accessory gene alleles.
allele_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds')

allele_jaccard_summary <- data.frame(species = names(allele_relabun),
                                     datatype = 'Allele',
                                     mean_jaccard = NA,
                                     sd_jaccard = NA)
rownames(allele_jaccard_summary) <- names(allele_relabun)


allele_count_summary <- data.frame(species = names(allele_relabun),
                                   mean_per_sample = NA,
                                   sd_per_sample = NA)
rownames(allele_count_summary) <- names(allele_relabun)


for (sp in names(allele_relabun)) {
  
  print(sp)
  
  sp_strain_samples <- colnames(strain_binary[[sp]])
  
  raw_out <- parallel::mclapply(X = names(allele_relabun[[sp]]), mc.cores = 40,
                                           FUN = function(gene) {
                                             
                                             sample_subset <- intersect(allele_relabun[[sp]][[gene]]$sample, sp_strain_samples)
                                             if (length(sample_subset) <= 2) { return(NA) }
                                             
                                             allele_binary <- allele_relabun[[sp]][[gene]][which(allele_relabun[[sp]][[gene]]$sample %in% sample_subset), ]
                                             if (length(unique(allele_binary$strain)) == 1) { return(NA) } # Skip if there is only one allele.
                                             
                                             allele_binary[which(allele_binary$community > 0), 'community'] <- 1
                                             allele_binary_wide <- reshape2::dcast(data = allele_binary, formula = sample ~ strain, value.var = 'community', fill = 0, drop = FALSE)[, -1, drop = FALSE]
                                             allele_binary_dist <- c(stats::dist(allele_binary_wide, method = "binary", diag = FALSE, upper = FALSE))
                                             alleles_per_sample <- rowSums(allele_binary_wide)
                                             
                                             return(list(mean_dist = mean(allele_binary_dist),
                                                         mean_per_sample = mean(alleles_per_sample)))
                                             
                                            })

  all_mean_jaccard_dist <- sapply(1:length(raw_out), function(i) { if(identical(NA, raw_out[[i]])) { return(NA) } else { return(raw_out[[i]]$mean_dist) }})
  allele_jaccard_summary[sp, 'mean_jaccard'] <- mean(all_mean_jaccard_dist, na.rm = TRUE)
  allele_jaccard_summary[sp, 'sd_jaccard'] <- sd(all_mean_jaccard_dist, na.rm = TRUE)
  
  all_mean_per_sample <- sapply(1:length(raw_out), function(i) { if(identical(NA, raw_out[[i]])) { return(NA) } else { return(raw_out[[i]]$mean_per_sample) }})
  allele_count_summary[sp, 'mean_per_sample'] <- mean(all_mean_per_sample, na.rm = TRUE)
  allele_count_summary[sp, 'sd_per_sample'] <- sd(all_mean_per_sample, na.rm = TRUE)

}

rownames(allele_jaccard_summary) <- NULL
rownames(allele_count_summary) <- NULL

combined_dist_summary <- do.call(rbind,
                                 list(species_jaccard_summary,
                                      strain_jaccard_summary,
                                      gene_jaccard_summary,
                                      allele_jaccard_summary))

write.table(x = combined_dist_summary,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/jaccard_dist_summary_by_datatype.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(x = allele_count_summary,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/accessory_allele_per_sample_summary.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
