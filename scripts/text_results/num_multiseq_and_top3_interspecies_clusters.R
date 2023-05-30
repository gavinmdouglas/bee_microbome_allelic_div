rm(list = ls(all.names = TRUE))

# Get number of multi-seq (and multi-species) clusters.
# Also, determine the top three species-species pairs.

cluster_map <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/cluster_members.tsv.gz",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cluster_tallies <- table(cluster_map$cluster)

singleton_clusters <- names(cluster_tallies)[which(cluster_tallies == 1)]
multi_seq_clusters <- names(cluster_tallies)[which(cluster_tallies > 1)]

# Total clusters:
length(cluster_tallies)

# Number multi-seq clusters:
length(multi_seq_clusters)

# Percentage of multi-seq clusters
round((length(multi_seq_clusters) / length(cluster_tallies)) * 100, 2)

# Number of genes in multi-seq clusters
sum(cluster_tallies[multi_seq_clusters])

cluster_map_multi_seq <- cluster_map[which(cluster_map$cluster %in% multi_seq_clusters), ]

species <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz",
                      header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
cluster_map_multi_seq$species <- NA
for (sp in species) {
  cluster_map_multi_seq[grep(sp, cluster_map_multi_seq$gene), "species"] <- sp
}

cross_species_clusters <- c()
cross_species_cluster_members <- c()
for (multi_seq_cluster in multi_seq_clusters) {
  cluster_map_subset <- cluster_map_multi_seq[which(cluster_map_multi_seq$cluster == multi_seq_cluster), ]
  if (length(unique(cluster_map_subset$species)) > 1) {
    cross_species_clusters <- c(cross_species_clusters, multi_seq_cluster)
    cross_species_cluster_members <- c(cross_species_cluster_members, paste(sort(unique(cluster_map_subset$species)), collapse = ","))
  }
}

# Number of cross-species clusters:
length(cross_species_clusters)

head(sort(table(cross_species_cluster_members), decreasing = TRUE), n = 10)

