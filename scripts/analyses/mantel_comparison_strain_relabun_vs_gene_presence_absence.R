rm(list = ls(all.names = TRUE))

# Per-species, summarize to what extent the variation in the gene presence/absence matrix (i.e., restricted to genes of that species)
# are is associated with variation in the presence/absence profile of strains (*not* relative abundances).
# This is primarily meant to be a sanity check, as of course these should agree to *some* extent.

library(vegan)

strainfacts_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')

straingst_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')

gene_cooccur <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv.gz',
                           header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)

passed_accessory_genes <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/all_species_pangenome_reference.trimmed.nocore.singletons.presentspecies.bed',
                                     header = FALSE, sep = '\t', stringsAsFactors = FALSE)$V1

intersecting_species <- intersect(names(strainfacts_relabun), names(straingst_relabun))
all_datasets <- c("Ellegaard2019", "Ellegaard2020", "Sun2022", "Wu2021", "Zhang2022")

raw_out <- list()
for (sp in intersecting_species) {
  
  print(sp)
  
  for (d in all_datasets) {

    print(d)
    
    if (! d %in% names(strainfacts_relabun[[sp]]) | ! d %in% names(straingst_relabun[[sp]])) { next }
    
    sp_strain_samples <- intersect(colnames(strainfacts_relabun[[sp]][[d]]), colnames(straingst_relabun[[sp]][[d]]))

    if (length(sp_strain_samples) <= 2) { next }
    
    sp_accessory <- grep(sp, passed_accessory_genes, value = TRUE)
    sp_accessory_called <- intersect(sp_accessory, rownames(gene_cooccur))
    if (length(sp_accessory_called) <= 50) { next }
    
    sp_gene_presence <- gene_cooccur[sp_accessory_called, sp_strain_samples]
    sp_gene_presence <- sp_gene_presence[which(rowSums(sp_gene_presence) > 0), , drop = FALSE]
    sp_gene_presence_dist <-  stats::dist(t(sp_gene_presence), method = "binary", diag = FALSE, upper = FALSE)
    
    strainfacts_binary <- strainfacts_relabun[[sp]][[d]][, sp_strain_samples, drop = FALSE]
    strainfacts_binary <- strainfacts_binary[which(rowSums(strainfacts_binary) > 0), , drop = FALSE]
    strainfacts_binary[strainfacts_binary > 0] <- 1
  
    straingst_binary <- straingst_relabun[[sp]][[d]][, sp_strain_samples, drop = FALSE]
    straingst_binary <- straingst_binary[which(rowSums(straingst_binary) > 0), , drop = FALSE]
    straingst_binary[straingst_binary > 0] <- 1
    
    if (nrow(strainfacts_binary) <= 2) { next }
    if (nrow(straingst_binary) <= 2) { next }
    
    strainfacts_dist <- stats::dist(t(strainfacts_binary), method = 'binary')
    straingst_dist <- stats::dist(t(straingst_binary), method = 'binary')
    
    strainfacts_vs_gene_mantel <- vegan::mantel(xdis = strainfacts_dist,
                                                ydis = sp_gene_presence_dist,
                                                method = 'kendall',
                                                parallel = 40)
    
    straingst_vs_gene_mantel <- vegan::mantel(xdis = straingst_dist,
                                              ydis = sp_gene_presence_dist,
                                              method = 'kendall',
                                              parallel = 40)

    raw_out[[paste(sp, d)]] <- data.frame(species=sp,
                                          dataset=d,
                                          num_samples=length(sp_strain_samples),
                                          num_genes=length(sp_accessory_called),
                                          strainfacts_mantel_kendall=strainfacts_vs_gene_mantel$statistic,
                                          strainfacts_mantel_kendall_p=strainfacts_vs_gene_mantel$signif,
                                          strainfacts_num_strains=nrow(strainfacts_binary),
                                          straingst_mantel_kendall=straingst_vs_gene_mantel$statistic,
                                          straingst_mantel_kendall_p=straingst_vs_gene_mantel$signif,
                                          straingst_num_strains=nrow(straingst_binary))
  }

}

mantel_summary_tab <- do.call(rbind, raw_out)
  
write.table(x = mantel_summary_tab,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/strain_vs_gene_presence_mantel.tsv',
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
