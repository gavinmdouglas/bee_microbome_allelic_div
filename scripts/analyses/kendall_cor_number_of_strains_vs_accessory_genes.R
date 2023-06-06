rm(list = ls(all.names = TRUE))

# Compute Kendall correlation between the number of accessory genes and number of strains per sample.
# Also, just provide the means and sd's of these values per species, which is easier to interpret.
# Note that these values are restricted to samples for which there is at least one strain called.

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1

strain_binary <- lapply(strain_relabun, function(x) { x[x > 0] <- 1; return(x) })

# Table to retain key information.
summary_tab <- data.frame(matrix(NA,
                                        nrow = length(strain_relabun),
                                        ncol = 7))
colnames(summary_tab) <- c('species', 'cor_kendall', 'cor_p',
                                  'mean_strains_per_sample', 'sd_strains_per_sample',
                                  'mean_acc.genes_per_sample', 'sd_acc.genes_per_sample')
rownames(summary_tab) <- names(strain_relabun)

for (sp in names(strain_binary)) {

  sp_strain_samples <- colnames(strain_binary[[sp]])
  sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
  sp_accessory_called <- intersect(sp_accessory, rownames(gene_cooccur))
  sp_gene_presence <- gene_cooccur[sp_accessory_called, sp_strain_samples]

  summary_tab[sp, 'species'] <- sp
  
  num_genes <- colSums(sp_gene_presence > 0)
  num_strains <- colSums(strain_binary[[sp]] > 0)
  cor_out <- cor.test(num_strains, num_genes, method = 'kendall', exact = FALSE)
  summary_tab[sp, c('cor_kendall', 'cor_p')] <- c(cor_out$estimate, cor_out$p.value)
  
  summary_tab[sp, c('mean_strains_per_sample', 'sd_strains_per_sample')] <- c(mean(num_strains), sd(num_strains))
  summary_tab[sp, c('mean_acc.genes_per_sample', 'sd_acc.genes_per_sample')] <- c(mean(num_genes), sd(num_genes))

}

write.table(x = summary_tab, file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_count_cor_summary.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
