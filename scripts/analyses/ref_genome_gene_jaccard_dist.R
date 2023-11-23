rm(list = ls(all.names = TRUE))

# Compute distance in accessory gene content between reference genomes of the same species.

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

raw_dist_by_sp <- list()

for (sp in all_species) {

  if (sp == 'Apilactobacillus_apinorum') { next }
  
  print(sp)
  
  sp_panaroo_file <- paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_panaroo/', sp, '/gene_presence_absence.csv.gz', sep = '')
  sp_gene_to_genome <- read.table(file = sp_panaroo_file, header = TRUE, sep = ',', stringsAsFactors = FALSE, quote = '', comment.char = '')
  rownames(sp_gene_to_genome) <- sp_gene_to_genome$Gene
  sp_gene_to_genome <- sp_gene_to_genome[, -c(1, 2, 3)]
  
  # Ignore core genes.
  sp_gene_to_genome <- sp_gene_to_genome[which(rowSums(sp_gene_to_genome == '') > 0), ]

  # Ignore any multi-copy genes:
  multi_copy_rows_i <- integer()
  for (row_i in 1:nrow(sp_gene_to_genome)) {
    if (length(grep(';', as.character(sp_gene_to_genome[row_i, , drop = FALSE]))) > 0) {
      multi_copy_rows_i <- c(multi_copy_rows_i, row_i)
    }
  }
  
  if (length(multi_copy_rows_i) > 0) {
    sp_gene_to_genome <- sp_gene_to_genome[-multi_copy_rows_i, ]
  }
                                                                              
  sp_gene_to_genome_binary <- data.frame(matrix(0,
                                                nrow = nrow(sp_gene_to_genome),
                                                ncol = ncol(sp_gene_to_genome)))
  rownames(sp_gene_to_genome_binary) <- rownames(sp_gene_to_genome)
  colnames(sp_gene_to_genome_binary) <- colnames(sp_gene_to_genome)

  sp_gene_to_genome_binary[sp_gene_to_genome != ''] <- 1
  
  acc_gene_jaccard <- as.matrix(stats::dist(x = t(sp_gene_to_genome_binary), method = 'binary'))
  
  num_genomes <- nrow(acc_gene_jaccard)
  
  sp_combined_dist <- data.frame(matrix(NA,
                                        nrow = choose(n = num_genomes, k = 2),
                                        ncol = 4))
  colnames(sp_combined_dist) <- c('species', 'genome1', 'genome2', 'acc_gene_jaccard')
  sp_combined_dist$species <- sp

  row_i <- 1
  for (i in 1:(num_genomes - 1)) {
  
    genome1 <- colnames(acc_gene_jaccard)[i]

    for (j in (i + 1):num_genomes) {
    
      genome2 <- colnames(acc_gene_jaccard)[j]
      
      sp_combined_dist[row_i, c('genome1', 'genome2')] <- c(genome1, genome2)
      sp_combined_dist[row_i, 'acc_gene_jaccard'] <- acc_gene_jaccard[genome1, genome2]
      
      row_i <- row_i + 1
      
    }
  }
  
  raw_dist_by_sp[[sp]] <- sp_combined_dist
}

combined_dist <- do.call(rbind, raw_dist_by_sp)

rownames(combined_dist) <- NULL

write.table(x = combined_dist,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/pairwise_genome_acc_gene_jaccard.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
