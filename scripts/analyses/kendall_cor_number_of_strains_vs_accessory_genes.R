rm(list = ls(all.names = TRUE))

# Compute Kendall correlation between the number of accessory genes and number of strains per sample.
# Also, just provide the means and sd's of these values per species, which is easier to interpret.
# Note that these values are restricted to samples for which there is at least one strain called.

strain_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1
raw_out <- list()

for (sp in names(strain_relabun)) {
  
  for (d in names(strain_relabun[[sp]])) {

  sample_subset <- colnames(strain_relabun[[sp]][[d]])
  gene_subset <- grep(sp, passed_accessory_genes, value = TRUE)
  sp_accessory_called <- intersect(gene_subset, rownames(gene_cooccur))
  sp_gene_presence <- gene_cooccur[sp_accessory_called, sample_subset]
  sp_gene_presence <- sp_gene_presence[which(rowSums(sp_gene_presence) > 0), ]
  
  num_genes <- colSums(sp_gene_presence > 0)
  num_strains <- colSums(strain_relabun[[sp]][[d]] > 0)
  
  cor_out <- cor.test(num_strains, num_genes, method = 'kendall', exact = FALSE)
  
  raw_out[[paste(sp, d, sep = '_')]] <- data.frame(species = sp,
                                                   dataset = d,
                                                   cor_kendall = cor_out$estimate,
                                                   cor_p = cor_out$p.value,
                                                   mean_strains_per_sample = mean(num_strains),
                                                   sd_strains_per_sample = sd(num_strains),
                                                   mean_acc.genes_per_sample = mean(num_genes),
                                                   sd_acc.genes_per_sample = sd(num_genes))
  }
}

summary_tab <- do.call(rbind, raw_out)

write.table(x = summary_tab,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_count_cor_summary.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
