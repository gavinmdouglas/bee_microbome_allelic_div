rm(list = ls(all.names = TRUE))

# Compare gene vs species tree distances (subset to the same tips).
# Also, compute Fritz and Purvis D statistic for each accessory gene
# within the range of >=25% and <= 75% of genomes.

library(ape)
library(caper)
library(parallel)
library(TreeDist)

compare_gene_vs_sp_tree <- function(gene_treefile, sp_tree, species_name) {
  
  gene_tree <- ape::read.tree(gene_treefile)

  gene_name <- gsub('.tree$', '', basename(gene_treefile))

  df_to_return <- data.frame(species=species_name,
                            gene=gene_name,
                            num_tips = length(gene_tree$tip.label),
                            num_internal_nodes = length(gene_tree$node.label),
                            tree_dist = NA,
                            tree_dist_subset5 = NA,
                            prevalence = NA,
                            D = NA,
                            D_p1 = NA,
                            D_p0 = NA)

  if (length(gene_tree$tip.label) < 3) { return(df_to_return) }
  if (length(setdiff(gene_tree$tip.label, rownames(sp_gene_to_genome))) > 0) { stop('Error: Tip IDs not found in Panaroo table.') }
  if (length(which(duplicated(sp_gene_to_genome[gene_tree$tip.label, 'gff_file']))) > 0) {  return(df_to_return) }

  gene_tree$tip.label <- sp_gene_to_genome[gene_tree$tip.label, 'gff_file']

  missing_ref_genomes <- setdiff(sp_tree$tip.label, gene_tree$tip.label)
  prepped_sp_tree <- ape::drop.tip(phy = sp_tree, tip = missing_ref_genomes)
  
  tree_dist <- TreeDist::TreeDistance(tree1 = gene_tree,
                                      tree2 = prepped_sp_tree)
  
  df_to_return$tree_dist <- tree_dist
  
  # If more than 5 tips, also try subsetting the tree to 5 tips to see how this affects things
  # (which could be useful for enabling more standardized comparisons across species!).
  # Do this for 10 replicates (or however are possible) and take the mean.
  if (length(gene_tree$tip.label) > 5) {
    
    five_tip_combos <- combn(gene_tree$tip.label, 5)
    
    if (ncol(five_tip_combos) > 10) {
      five_tip_combos <- five_tip_combos[, sample(x = 1:ncol(five_tip_combos), size = 10)] 
    }
    
    rep_dist <- numeric()
    for (i in 1:ncol(five_tip_combos)) {
      rep_tips <- five_tip_combos[, i]
      tips_to_drop <- setdiff(gene_tree$tip.label, rep_tips)

      rep_gene_tree <- ape::drop.tip(phy = gene_tree, tip = tips_to_drop)
      rep_sp_tree <- ape::drop.tip(phy = prepped_sp_tree, tip = tips_to_drop)
      
      rep_dist <- c(rep_dist,
                    TreeDist::TreeDistance(tree1 = rep_gene_tree,
                                          tree2 = rep_sp_tree))
      
    }
    
    df_to_return$tree_dist_subset5 <- mean(rep_dist)
    
  } else if (length(gene_tree$tip.label) == 5) {
     df_to_return$tree_dist_subset5 <- tree_dist
  }
  
  df_to_return$prevalence <- length(gene_tree$tip.label) / length(sp_tree$tip.label)

  if (df_to_return$prevalence >= 0.15 & df_to_return$prevalence <= 0.85) {
    in_data <- data.frame(tips = sp_tree$tip.label, presence = 0)
    in_data[which(in_data$tips %in% gene_tree$tip.label), 'presence'] <- 1

    sp_tree_D_prep <- midpoint_root_make_tree_binary_and_no_zero_branch(sp_tree)

    raw_phyloD <- caper::phylo.d(data = in_data,
                                 phy = sp_tree_D_prep,
                                 names.col = tips,
                                 binvar = presence)

    df_to_return$D <- raw_phyloD$DEstimate
    df_to_return$D_p1 <- raw_phyloD$Pval1
    df_to_return$D_p0 <- raw_phyloD$Pval0
  }

  return(df_to_return)

}

midpoint_root_make_tree_binary_and_no_zero_branch <- function(in_tree) {

  if (! ape::is.rooted(in_tree)) {
    in_tree <- phytools::midpoint.root(in_tree)
  }

  if (! is.binary(in_tree)) {
    in_tree <- TreeTools::MakeTreeBinary(in_tree)
  }

  if (min(in_tree$edge.length) == 0) {
    min_non0 <- min(in_tree$edge.length[which(in_tree$edge.length > 0)])
    in_tree$edge.length[which(in_tree$edge.length == 0)] <- min_non0
  }
  
  return(in_tree)
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
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_vs_species_dist_and_D.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
