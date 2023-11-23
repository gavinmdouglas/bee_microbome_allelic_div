rm(list = ls(all.names = TRUE))

library(reshape2)
library(ape)
library(castor)

# Many genomes were dropped because they are redundant. Identify these before processing trees.
all_species <- read.table('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt', stringsAsFactors = FALSE)$V1

genomes_to_exclude <- list()
for (sp in all_species) {
  sp_genomes <- gsub('\\.fna$', '', list.files(path = paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes/', sp, sep = ''), pattern = '.fna'))
  
  sp_genomes_retained_files <- read.table(paste('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/strainge/kmer_comparison/references_to_keep/', sp, '.txt', sep = ''),
                                          stringsAsFactors = FALSE)$V1
  sp_genomes_retained <- gsub('\\.hdf5$', '', basename(sp_genomes_retained_files))
  
  if (length(setdiff(sp_genomes_retained, sp_genomes)) > 0) {
    stop('Error - Retained genomes not found in orig set.')
  }
  
  genomes_to_exclude[[sp]] <- setdiff(sp_genomes, sp_genomes_retained)
}

trees <- list()
tree_paths <- list.files('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strain_trees/',
                         pattern = '.tree',
                         full.names = TRUE)

for (tree_path in tree_paths) {
  tree_id <- gsub('.tree', '', basename(tree_path))
  sp_id <- gsub('\\..*$', '', tree_id)
  trees[[tree_id]] <- ape::read.tree(tree_path)
  
  if (length(genomes_to_exclude[[sp_id]]) > 0) {
    trees[[tree_id]] <- ape::drop.tip(trees[[tree_id]], genomes_to_exclude[[sp_id]])
  }
}

species <- sort(unique(sapply(tree_paths, function(x) { strsplit(basename(x), '\\.')[[1]][1] })))

# Compute observed percentile vs. real data

strain_dists_RAW <- list()

for (tree_id in names(trees)) {
  
  print(tree_id)
  
  tree_id_split <- strsplit(tree_id, '\\.')[[1]]
  sp <- tree_id_split[1]
  dataset <- tree_id_split[2]
  
  tip_dists <- data.frame(castor::get_all_pairwise_distances(trees[[tree_id]],
                                                             only_clades = 1:length(trees[[tree_id]]$tip.label)))
  tip_labels <- trees[[tree_id]]$tip.label
  tip_labels <- gsub(dataset, '', tip_labels)
  tip_labels <- gsub('^\\.', '', tip_labels)
  colnames(tip_dists) <- tip_labels
  rownames(tip_dists) <- tip_labels
  
  tip_dists$tip <- tip_labels
  tip_dists_long <- reshape2::melt(tip_dists, id.vars = 'tip')
  tip_dists_long <- tip_dists_long[which(tip_dists_long$tip != tip_dists_long$variable), , drop = FALSE]
  
  tip_dists_long$pair <- apply(tip_dists_long[, c(1, 2)], 1, function(row) paste(sort(row), collapse = ','))
  tip_dists_long <- tip_dists_long[, c('pair' , 'value'), drop = FALSE]
  tip_dists_long <- tip_dists_long[which(! duplicated(tip_dists_long)), , drop = FALSE]
  rownames(tip_dists_long) <- tip_dists_long$pair

  strainfacts_strain_abun <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/community_relabun/', dataset, '.', sp, '.tsv.gz', sep = ''),
                                        header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
  
  all_strainfacts_strains <- sort(rownames(strainfacts_strain_abun))
  all_strainge_strains <- sort(setdiff(tip_labels, all_strainfacts_strains))
  
  possible_strainfacts_combos <- list()
  possible_strainge_combos <- list()
  
  out_df <- data.frame(species = sp,
                       dataset = dataset,
                       sample = colnames(strainfacts_strain_abun),
                       num_strainfacts_strains = NA,
                       num_strainge_strains = NA,
                       mean_obs_dist = NA,
                       num_background_combos = NA,
                       mean_background_dist = NA,
                       background_dist_percentile = NA)
  rownames(out_df) <- out_df$sample

  for (sample_id in colnames(strainfacts_strain_abun)) {

    strainfacts_strains <- rownames(strainfacts_strain_abun)[which(strainfacts_strain_abun[, sample_id] > 0)]

    strainge_strains <- read.table(paste('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/strainge/straingst_out/', sample_id, '_', sp, '.strains.tsv', sep = ''),
                                   header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)$strain

    strain_combos <- expand.grid(strainge_strains, strainfacts_strains)
    strain_combos <- paste(strain_combos$Var1, strain_combos$Var2, sep = ',')
    if (length(setdiff(strain_combos, tip_dists_long$pair)) > 0) { stop('Observed strains not in expected set -- must be an issue with strain pair ID format.') }

    # Compute observed mean distance.
    obs_mean_dist <- mean(tip_dists_long[strain_combos, 'value'])

    # Then compute background distribution of mean distances possible with all possible combinations of strain IDs, of the same number from each tool.
    num_strainfacts_strains <- length(strainfacts_strains)
    num_strainfacts_strains_char <- as.character(num_strainfacts_strains)
    if (! num_strainfacts_strains_char %in% names(possible_strainfacts_combos)) {
      possible_strainfacts_combos[[num_strainfacts_strains_char]] <- combn(all_strainfacts_strains, num_strainfacts_strains)
    }

    num_strainge_strains <- length(strainge_strains)
    num_strainge_strains_char <- as.character(num_strainge_strains)

    if (! num_strainge_strains_char %in% names(possible_strainge_combos)) {
      possible_strainge_combos[[num_strainge_strains_char]] <- combn(all_strainge_strains, num_strainge_strains)
    }

    num_strainfacts_background_combos <- ncol(possible_strainfacts_combos[[num_strainfacts_strains_char]])
    num_strainge_background_combos <- ncol(possible_strainge_combos[[num_strainge_strains_char]])
    
    num_background_combos <- num_strainfacts_background_combos * num_strainge_background_combos
    
    if (num_background_combos > 1024) {
      
      if (num_strainfacts_background_combos >= 32 & num_strainge_background_combos >= 32) {
        num_background_combos <- 1024
        rep_strainfacts_indices <- sample(1:ncol(possible_strainfacts_combos[[num_strainfacts_strains_char]]), 32, replace = FALSE)
        rep_strainge_indices <- sample(1:ncol(possible_strainge_combos[[num_strainge_strains_char]]), 32, replace = FALSE)
        
      } else if (num_strainfacts_background_combos < num_strainge_background_combos) {
        num_strainge_background_combos <- ceiling(1024 / num_strainfacts_background_combos)
        num_background_combos <- num_strainfacts_background_combos * num_strainge_background_combos
        rep_strainfacts_indices <- 1:ncol(possible_strainfacts_combos[[num_strainfacts_strains_char]])
        rep_strainge_indices <- sample(1:ncol(possible_strainge_combos[[num_strainge_strains_char]]), num_strainge_background_combos, replace = FALSE)
        
      } else {
        num_strainfacts_background_combos <- ceiling(1024 / num_strainge_background_combos)
        num_background_combos <- num_strainfacts_background_combos * num_strainge_background_combos
        rep_strainfacts_indices <- sample(1:ncol(possible_strainfacts_combos[[num_strainfacts_strains_char]]), num_strainfacts_background_combos, replace = FALSE)
        rep_strainge_indices <- 1:ncol(possible_strainge_combos[[num_strainge_strains_char]])
      }
      
    } else {
      rep_strainfacts_indices <- 1:ncol(possible_strainfacts_combos[[num_strainfacts_strains_char]])
      rep_strainge_indices <- 1:ncol(possible_strainge_combos[[num_strainge_strains_char]])
    }
    
    background_means <- numeric()
    
    for (rep_i in rep_strainfacts_indices) {
      
      for (rep_j in rep_strainge_indices) {
        
        rep_strain_combos <- expand.grid(possible_strainge_combos[[num_strainge_strains_char]][, rep_j],
                                         possible_strainfacts_combos[[num_strainfacts_strains_char]][, rep_i])
        rep_strain_combos <- paste(rep_strain_combos$Var1, rep_strain_combos$Var2, sep = ',')
        
        background_means <- c(background_means, mean(tip_dists_long[rep_strain_combos, 'value']))
      }
      
    }
    
    background_dist_P <- length(which(background_means <= obs_mean_dist)) / num_background_combos
    
    if (background_dist_P == 0) {
      background_dist_P <- (length(which(background_means <= obs_mean_dist)) + 1) / num_background_combos
    }
    
    mean_background_dist <- mean(background_means)
    
    out_df[sample_id, c('num_strainfacts_strains',
                        'num_strainge_strains',
                        'mean_obs_dist',
                        'num_background_combos',
                        'mean_background_dist',
                        'background_dist_percentile')] <- c(num_strainfacts_strains,
                                                            num_strainge_strains,
                                                            obs_mean_dist,
                                                            num_background_combos,
                                                            mean_background_dist,
                                                            background_dist_P * 100)
    
  }
  
  strain_dists_RAW[[tree_id]] <- out_df
  
}

strain_dists_summary <- do.call(rbind, strain_dists_RAW) 

write.table(x = strain_dists_summary,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/strainfacts_vs_strainge_dist_summary.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
