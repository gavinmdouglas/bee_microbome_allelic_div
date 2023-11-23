rm(list = ls(all.names = TRUE))

library(ape)
library(parallel)
library(TreeDist)

compare_gene_vs_sp_tree <- function(gene_treefile, sp_tree, species_name) {

  gene_tree <- ape::read.tree(gene_treefile)
  
  gene_name <- gsub('.tree$', '', basename(gene_treefile))
  
  df_to_return <- data.frame(species=species_name,
                            gene=gene_name,
                            num_tips = length(gene_tree$tip.label),
                            num_internal_nodes = length(gene_tree$node.label),
                            tree_dist = NA)
  
  if (length(gene_tree$tip.label) < 3) { return(df_to_return) }
  if (length(setdiff(gene_tree$tip.label, rownames(sp_gene_to_genome))) > 0) { stop('Error: Tip IDs not found in Panaroo table.') }
  if (length(which(duplicated(sp_gene_to_genome[gene_tree$tip.label, 'gff_file']))) > 0) {  return(df_to_return) }
  
  gene_tree$tip.label <- sp_gene_to_genome[gene_tree$tip.label, 'gff_file']
  
  missing_ref_genomes <- setdiff(sp_tree$tip.label, gene_tree$tip.label)
  prepped_sp_tree <- ape::drop.tip(phy = sp_tree, tip = missing_ref_genomes)
  
  tree_dist <- TreeDist::DifferentPhylogeneticInfo(tree1 = gene_tree,
                                                   tree2 = prepped_sp_tree,
                                                   normalize = TRUE)
  df_to_return$tree_dist <- tree_dist

  return(df_to_return)

}

ref_treefiles <- list.files('/data1/gdouglas/projects/honey_bee/ref_genomes/core_genome_trees/iqtree_out/',
                            pattern = '.treefile$', full.names = TRUE)

raw_dist_by_sp <- list()

for (sp_treefile in ref_treefiles) {
  
  sp <- gsub('.treefile$', '', basename(sp_treefile))

  print(sp)
  
  sp_gene_info_file <- paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_panaroo/', sp, '/gene_data.csv.gz', sep = '')
  sp_gene_to_genome <- read.table(file = sp_gene_info_file, header = TRUE, sep = ',', stringsAsFactors = FALSE, quote = '', comment.char = '')
  sp_gene_to_genome <- sp_gene_to_genome[, c('gff_file', 'annotation_id')]
  rownames(sp_gene_to_genome) <- sp_gene_to_genome$annotation_id

  unique_ref_genomes <- unique(sp_gene_to_genome$gff_file)

  sp_tree <- ape::read.tree(sp_treefile)
  
  genomes_to_remove <- setdiff(sp_tree$tip.label, unique_ref_genomes)
  if (length(genomes_to_remove) > 0) {
    sp_tree <- ape::drop.tip(phy = sp_tree, tip = genomes_to_remove)
  }
  
  path_to_treefiles <- paste('/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/ortholog_pandora_prgs/',
                             sp,
                             '.ortholog_fastas_muscle_aligned_trees/',
                             sep = '')
  
  sp_gene_treefiles <- list.files(path_to_treefiles,
                                  pattern = '.tree$', full.names = TRUE)
  
  raw_out <- parallel::mclapply(X = sp_gene_treefiles,
                                FUN = compare_gene_vs_sp_tree,
                                sp_tree = sp_tree,
                                species_name = sp,
                                mc.cores = 50)
 
  raw_dist_by_sp[[sp]] <- do.call(rbind, raw_out)

}

combined_dist <- do.call(rbind, raw_dist_by_sp)

rownames(combined_dist) <- NULL

write.table(x = combined_dist,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_vs_species_dist.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
