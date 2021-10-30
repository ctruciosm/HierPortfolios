#' Hierarchical Risk Parity portfolio
#'
#' Performs the Hierarchical Risk Parity portfolio proposed strategy by De Prado (2016). Several linkage methods for the hierarchical clustering can be used, by default the `single` linkage is used.
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into distance matrix.
#' @param clustering.method Linkage method used in the hierarchical clustering. Allowed options are `single`, `complete`, `average` or `ward`. Default option is `single`.
#' @return portfolio weights
#' @seealso `HCAA_porfolio` and `HERC_Portfolio`
#' @export
HRP_Portfolio = function(covar, clustering.method= "single"){
  if (clustering.method %in% c("single", "complete", "average", "ward")) {
    if (clustering.method == "ward") {
      clustering.method = "ward.D2"
    }
  } else {
    return("ERROR: cclustering.method argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  # Stage 1: Tree clustering
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  # Stage 2: Quasi-Diagonalisation
  clusters_order <- fastcluster::hclust(euclidean_distance , method = clustering.method, members = NULL)$order
  # Stage 3: Recursive Bisection
  weights <- rep(1,ncol(covar))
  index <- list(clusters_order)
  while (length(index) > 0) {
    new.index <- list()
    for (i in index) {
      middle <- floor(length(i)/2)
      indexa <- i[1:middle]
      indexb <- i[-c(1:middle)]
      covar_clustera <- as.matrix(covar[indexa, indexa])
      covar_clusterb <- as.matrix(covar[indexb, indexb])
      weightsa <- 1/diag(covar_clustera) / sum(1/diag(covar_clustera))
      weightsb <- 1/diag(covar_clusterb) / sum(1/diag(covar_clusterb))
      variance_clustera <- weightsa %*% covar_clustera %*% weightsa
      variance_clusterb <- weightsb %*% covar_clusterb %*% weightsb
      alpha <- as.numeric(1 - variance_clustera/(variance_clustera + variance_clusterb))
      weights[indexa] <- weights[indexa] * alpha
      weights[indexb] <- weights[indexb] * (1 - alpha)
      if (length(indexa) > 1) new.index <- c(new.index, list(indexa))
      if (length(indexb) > 1) new.index <- c(new.index, list(indexb))
      index = new.index
    }
  }
  return(weights)
}
