HERC_Portfolio = function(covar, clustering.method = "ward"){
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
  weights_in_cluster <- elements_in_cluster
  weights_erc <- optimalPortfolio(Sigma = covar, control = list(type = 'erc', constraint = 'lo'))
  risk_contribution <- weights_erc * covar %*% weights_erc  / as.numeric(weights_erc %*% covar %*% weights_erc)
  # Bisection 1
  indexa <- elements_in_cluster[,1] == 1
  indexb <- elements_in_cluster[,1] == 2
  risk_contribution_clustera <- sum(risk_contribution[indexa])
  risk_contribution_clusterb <- sum(risk_contribution[indexb])
  # We use cluster's variance as a measure of risk
  w <- risk_contribution_clustera/(risk_contribution_clustera + risk_contribution_clusterb)
  weights_in_cluster[, 1] <- 1
  weights_in_cluster[indexa, 1] <- w
  weights_in_cluster[indexb, 1] <- 1 - w
  # Other Bisections
  if (n_clusters > 2) {
    for (i in 2:(n_clusters - 1)) {
      groups <- unique(elements_in_cluster[,i - 1])
      for (j in 1:i) {
        index <- elements_in_cluster[, i - 1] == groups[j]
        if (length(unique(elements_in_cluster[index, i])) > 1) {
          sub_groups <- unique(elements_in_cluster[index, i])
          indexa <- elements_in_cluster[ ,i] == sub_groups[1]
          indexb <- elements_in_cluster[ ,i] == sub_groups[2]
          risk_contribution_clustera <- sum(risk_contribution[indexa])
          risk_contribution_clusterb <- sum(risk_contribution[indexb])
          w <- risk_contribution_clustera/(risk_contribution_clustera + risk_contribution_clusterb)
          weights_in_cluster[,i] <- 1
          weights_in_cluster[indexa, i] <- w
          weights_in_cluster[indexb, i] <- 1 - w
          break
        } 
      }
    }
  }
  # Naive Risk Parity for elements in each cluster: inverse of variance
  naive_rp_weights <- rep(1,n_cols)
  for (j in 1:n_clusters) {
    index <- elements_in_cluster[, n_clusters - 1] == j
    covar_elements <- as.matrix(covar[index, index])
    naive_rp_weights[index] <- 1/diag(covar_elements) / sum(1/diag(covar_elements))
  }
  weights <- apply(weights_in_cluster, 1, prod) * naive_rp_weights
  return(weights)
}