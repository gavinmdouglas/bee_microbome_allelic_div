rm(list = ls(all.names = TRUE))

pairwise_matrix_mantel <- function(dist_mats) {
  
  if (length(dist_matrices) <= 1) {
    return(NULL)
  }
  
  sample_subset <- colnames(as.matrix(dist_mats[[1]]))
  for (i in 2:length(dist_mats)) {
    sample_subset <- intersect(sample_subset, colnames(as.matrix(dist_mats[[i]])))
  }
  
  all_mats <- list()
  for (i in 1:length(dist_mats)) {
    all_mats[[i]] <- as.matrix(dist_mats[[i]])[sample_subset, sample_subset]
  }
  
  mantel_r <- numeric()
  rep_i <- character()
  rep_j <- character()
  for (i in 1:(length(all_mats) - 1)) {
    for (j in (i + 1):length(dist_mats)) {
      mantel_r <- c(mantel_r, mantel(all_mats[[i]], all_mats[[j]])$statistic)
      rep_i <- c(rep_i, names(dist_mats)[i])
      rep_j <- c(rep_j,names(dist_mats)[j])
    }
  }
  
  return(data.frame(rep_i=rep_i, rep_j=rep_j, mantel_r=mantel_r))
}

d <- 'Ellegaard2019'
sp <- 'Bartonella_apis'

dist_matrices <- list()
for (subsample in c("subsample20", "subsample50")) {
  for (rep in replicates) {
    if (! is.null(strain_abun_0.1[[d]][[sp]][[subsample]][[rep]])) {
      dist_matrices[[paste(subsample, rep)]] <- vegdist(x=t(strain_abun_0.1[[d]][[sp]][[subsample]][[rep]]), method="bray")
    }
  }
}
pairwise_mantel_summary <- pairwise_matrix_mantel(dist_matrices)
g <- graph_from_data_frame(pairwise_mantel_summary, directed=FALSE)
# Plot the graph
plot(g, edge.width=E(g)$mantel_r * 2)
pairwise_mantel_summary
g <- graph_from_data_frame(pairwise_mantel_summary, directed=FALSE)
# Plot the graph

E(g)$color <- ifelse(E(g)$mantel_r > 0, "blue", "red")
E(g)$width <- abs(E(g)$mantel_r) * 5
plot(g)

closeness_centrality <- closeness(g)
closeness_centrality
eigenvector_centrality <- eigen_centrality(g)$vector
eigenvector_centrality

page_rank(g)$vector

pairwise_mantel_summary_long <- data.frame(rep_id = c(pairwise_mantel_summary$rep_i, pairwise_mantel_summary$rep_j),
                                           mantel_r = c(pairwise_mantel_summary$mantel_r, pairwise_mantel_summary$mantel_r))

mantel_r_mean <- aggregate(x = mantel_r ~ rep_id, FUN = mean, data = pairwise_mantel_summary_long)

mantel_r_mean[order(mantel_r_mean$mantel_r, decreasing = TRUE), ][1, 'rep_id']