#' Hierarchical Clustering-Based Asset Allocation
#'
#' Performs the Hierarchical Clustering-Based Asset Allocation strategy proposed by Raffinot (2017). Several linkage methods for the hierarchical clustering can be used, by default the "ward" linkage is used. The numbers of clusters is selected using the Gap index of Tibshirani et al. (2001).
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
#' @param linkage Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "ward".
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @param clusters Numbers of clusters. If NULL (default), the gap index is applied.
#' @return portfolio weights.
#' @seealso \code{HRP_Portfolio}
#' @references Raffinot, Thomas. "Hierarchical clustering-based asset allocation." The Journal of Portfolio Management 44.2 (2017): 89-99.
#' @references Tibshirani, Robert, Guenther Walther, and Trevor Hastie. "Estimating the number of clusters in a data set via the gap statistic." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 63.2 (2001): 411-423.
#' @aliases HCAA_Portfolio
#' @keywords HCAA
#' @examples
#' covar <- cov(daily_returns)
#' HCAA_Portfolio(covar)
#' @export
HCAA_Portfolio = function(covar, linkage = "ward", graph = FALSE, clusters = NULL) {
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
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - HCAA")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}
