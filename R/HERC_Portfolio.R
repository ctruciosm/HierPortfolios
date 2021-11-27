#' Hierarchical Equal Risk Contribution
#'
#' Performs the Hierarchical Equal Risk Contribution portfolio strategy proposed by Raffinot (2018). Several linkage methods for the hierarchical clustering can be used, by default the "ward" linkage is used. This function uses the variance as risk measure. The number of clusters is selected using the Gap index of Tibshirani et al. (2001). The implemenation follows Sjostrand and Nina (2020).
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
#' @param linkage Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "ward".
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @param clusters Numbers of clusters. If NULL (default), the gap index is applied.
#' @return portfolio weights.
#' @seealso \code{HRP_Portfolio} and \code{HCAA_Portfolio}
#' @references Raffinot, Thomas. "The hierarchical equal risk contribution portfolio." Available at SSRN 3237540 (2018).
#' @references Sjostrand, Daniel, and Nima Behnejad. "Exploration of Hierarchical Clustering in Long-Only Risk-Based Portfolio Optimization." Master Thesis, Copenhagen Business School (2020).
#' @references Tibshirani, Robert, Guenther Walther, and Trevor Hastie. "Estimating the number of clusters in a data set via the gap statistic." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 63.2 (2001): 411-423.
#' @aliases HERC_Portfolio
#' @keywords HERC
#' @examples
#' covar <- cov(daily_returns)
#' HERC_Portfolio(covar)
#' @export
HERC_Portfolio = function(covar, linkage = "ward", graph = FALSE, clusters = NULL) {
  if (linkage %in% c("single", "complete", "average", "ward")) {
    if (linkage == "ward") {
      linkage = "ward.D2"
    }
  } else {
    return("ERROR: linkage argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  clustering <- fastcluster::hclust(euclidean_distance , method = linkage, members = NULL)
  n_cols <- ncol(corre)
  if (is.null(clusters)) {
    fun_clus_num <- function(x,k) list(cluster = stats::cutree(fastcluster::hclust(stats::as.dist(x) , method = linkage, members = NULL), k))
    gap <- cluster::clusGap(as.matrix(euclidean_distance), FUN = fun_clus_num, K.max = floor(n_cols/2), B = 100)
    n_clusters <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = "Tibs2001SEmax")
    n_clusters <- max(2, n_clusters)
  } else{
    n_clusters <- clusters
  }
  print(sprintf("Number of clusters: %i", n_clusters))
  elements_in_cluster <- matrix(stats::cutree(clustering, 2:n_clusters), ncol = n_clusters - 1)
  weights_in_cluster <-  rep(1,n_cols)
# Bisection 1
  indexa <- elements_in_cluster[,1] == 1
  indexb <- elements_in_cluster[,1] == 2
  covar_clustera <- as.matrix(covar[indexa, indexa])
  covar_clusterb <- as.matrix(covar[indexb, indexb])
  weightsa <- 1/diag(covar_clustera) / sum(1/diag(covar_clustera))
  weightsb <- 1/diag(covar_clusterb) / sum(1/diag(covar_clusterb))
  # HERC does not mention how to compute the cluster variances.
  # Since HERC is based on HRP and HCAA, we follow the HRP idea and compute the cluster variances as:
  risk_contribution_clustera <- as.numeric(weightsa %*% covar_clustera %*% weightsa)
  risk_contribution_clusterb <- as.numeric(weightsb %*% covar_clusterb %*% weightsb)
  alphaa <- 1 - risk_contribution_clustera/(risk_contribution_clustera + risk_contribution_clusterb)
  alphab <- 1 - alphaa
  weights_in_cluster[indexa] <- alphaa*weights_in_cluster[indexa]
  weights_in_cluster[indexb] <- alphab*weights_in_cluster[indexb]
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
          covar_clustera <- as.matrix(covar[indexa, indexa])
          covar_clusterb <- as.matrix(covar[indexb, indexb])
          weightsa <- 1/diag(covar_clustera) / sum(1/diag(covar_clustera))
          weightsb <- 1/diag(covar_clusterb) / sum(1/diag(covar_clusterb))
          risk_contribution_clustera <- as.numeric(weightsa %*% covar_clustera %*% weightsa)
          risk_contribution_clusterb <- as.numeric(weightsb %*% covar_clusterb %*% weightsb)
          alphaa <- 1 - risk_contribution_clustera/(risk_contribution_clustera + risk_contribution_clusterb)
          alphab <- 1 - alphaa
          weights_in_cluster[indexa] <- alphaa*weights_in_cluster[indexa]
          weights_in_cluster[indexb] <- alphab*weights_in_cluster[indexb]
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
  weights <- weights_in_cluster * naive_rp_weights
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - HERC")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}
