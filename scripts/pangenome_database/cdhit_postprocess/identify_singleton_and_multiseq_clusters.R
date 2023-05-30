rm(list = ls(all.names = TRUE))

# Figure out which genes are singletons and which are in multi-seq clusters, the latter of which should be ignored.

# Note that the raw CD-HIT output was processed with these two scripts:
#   get_cluster_member_breakdown.py
#   get_cluster_rep_id_map.py

cluster_map <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/cluster_members.tsv.gz",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cluster_tallies <- table(cluster_map$cluster)

singleton_clusters <- names(cluster_tallies)[which(cluster_tallies == 1)]
multi_seq_clusters <- names(cluster_tallies)[which(cluster_tallies > 1)]

singleton_cluster_members <- cluster_map[which(cluster_map$cluster %in% singleton_clusters), "gene"]
multi.seq_cluster_members <- cluster_map[which(cluster_map$cluster %in% multi_seq_clusters), "gene"]

write.table(x = singleton_cluster_members,
            file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/singleton_members.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(x = multi.seq_cluster_members,
            file = "/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/cdhit_out/multi.seq_members.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
