rm(list = ls(all.names = TRUE))

# Per-species, summarize to what extent the variation in the gene presence/absence matrix (i.e., restricted to genes of that species)
# are is associated with variation in the presence/absence profile of strains (*not* relative abundances).
# This is primarily meant to be a sanity check, as of course these should agree to *some* extent.

library(vegan)

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1

strain_binary <- lapply(strain_relabun, function(x) { x[x > 0] <- 1; return(x) })

# Then compute Jaccard distances.
strain_binary_dist <- lapply(strain_binary, function(x) { stats::dist(t(x), method = "binary", diag = FALSE, upper = FALSE) })

# Table to retain key information.
mantel_summary_tab <- data.frame(matrix(NA,
                                        nrow = length(strain_relabun),
                                        ncol = 6))
colnames(mantel_summary_tab) <- c('species', 'mantel_kendall', 'mantel_p',
                                  'num_strains', 'num_samples', 'num_genes')
rownames(mantel_summary_tab) <- names(strain_relabun)

mantel_out <- list()
for (sp in names(strain_binary_dist)) {

  sp_strain_samples <- colnames(strain_binary[[sp]])
  sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
  sp_accessory_called <- intersect(sp_accessory, rownames(gene_cooccur))
  sp_gene_presence <- gene_cooccur[sp_accessory_called, sp_strain_samples]
  sp_gene_presence_dist <-  stats::dist(t(sp_gene_presence), method = "binary", diag = FALSE, upper = FALSE)
  
  mantel_out[[sp]] <- vegan::mantel(xdis = strain_binary_dist[[sp]],
                                    ydis = sp_gene_presence_dist,
                                    method = 'kendall',
                                    parallel = 40)
  
  mantel_summary_tab[sp, 'species'] <- sp
  mantel_summary_tab[sp, 'mantel_kendall'] <- mantel_out[[sp]]$statistic
  mantel_summary_tab[sp, 'mantel_p'] <- mantel_out[[sp]]$signif
  mantel_summary_tab[sp, 'num_strains'] <- nrow(strain_binary[[sp]])
  mantel_summary_tab[sp, 'num_samples'] <- ncol(strain_binary[[sp]])
  mantel_summary_tab[sp, 'num_genes'] <- length(sp_accessory_called)

}

write.table(x = mantel_summary_tab, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_presence_mantel.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
