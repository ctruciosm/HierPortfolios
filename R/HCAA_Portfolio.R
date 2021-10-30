HCAA_Portfolio = function(covar, clustering.method = "ward"){
  if (clustering.method %in% c("single", "complete", "average", "ward")) {
    if (clustering.method == "ward") {
      clustering.method = "ward.D2"
    }
  } else {
    return("ERROR: clustering.method argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  corre <- cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  clustering <- fastcluster::hclust(euclidean_distance , method = clustering.method, members = NULL)
  n_cols <- ncol(corre)
  fun_clus_num <- function(x,k) list(cluster = cutree(fastcluster::hclust(as.dist(x) , method = "ward.D2", members = NULL), k)) 
  gap <- clusGap(as.matrix(euclidean_distance), FUN = fun_clus_num, K.max = floor(n_cols/2), B = 100)
  n_clusters <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = "Tibs2001SEmax")
  n_clusters <- max(2, n_clusters)
  elements_in_cluster <- matrix(cutree(clustering, 2:n_clusters), ncol = n_clusters - 1)
  # Using cluster's hierarchy 
  indexa <- elements_in_cluster[,1] == 1
  indexb <- elements_in_cluster[,1] == 2
  weights <- rep(1, n_cols)
  weights[indexa] <- weights[indexa]/2
  weights[indexb] <- weights[indexb]/2
  if (n_clusters > 2) {
    for (i in 2:(n_clusters - 1)) {
      groups <- unique(elements_in_cluster[,i - 1])
      for (j in 1:i) {
        index <- elements_in_cluster[, i - 1] == groups[j]
        if (length(unique(elements_in_cluster[index, i])) > 1) {
          sub_groups <- unique(elements_in_cluster[index, i])
          indexa <- elements_in_cluster[ ,i] == sub_groups[1]
          indexb <- elements_in_cluster[ ,i] == sub_groups[2]
          weights[indexa] <- weights[indexa]/2
          weights[indexb] <- weights[indexb]/2
          break
        } 
      }
    }
  }
  for (j in 1:n_clusters) {
    index <- elements_in_cluster[, n_clusters - 1] == j
    weights[index] <- weights[index]/sum(index)
  }
  return(weights)
}