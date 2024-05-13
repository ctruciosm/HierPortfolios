#' Hierarchical Equal Risk Contribution
#'
#' Performs the Hierarchical Equal Risk Contribution portfolio strategy proposed by Raffinot (2018). Several linkage methods for the hierarchical clustering can be used, by default the "ward" linkage is used. This function uses the variance as risk measure. The number of clusters is selected using the Gap index of Tibshirani et al. (2001). The implemenation follows Sjostrand and Nina (2020).
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
#' @param linkage Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "ward".
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @param clusters Numbers of clusters. If NULL (default), the gap index is applied.
#' @return portfolio weights.
#' @seealso \code{HRP_Portfolio}, \code{HCAA_Portfolio} and \code{DHRP_Portfolio}
#' @references Raffinot, Thomas. "The hierarchical equal risk contribution portfolio." Available at SSRN 3237540 (2018).
#' @references Tibshirani, Robert, Guenther Walther, and Trevor Hastie. "Estimating the number of clusters in a data set via the gap statistic." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 63.2 (2001): 411-423.
#' @aliases HERC_Portfolio
#' @keywords HERC
#' @examples
#' covar <- cov(daily_returns)
#' HERC_Portfolio(covar)
#' @author Carlos Trucios and Moon Jun Kwon
#' @export
HERC_Portfolio = function(covar, linkage = "ward", graph = FALSE, clusters = NULL) {
  if (linkage %in% c("single", "complete", "average", "ward")) {
    if (linkage == "ward") {
      linkage = "ward.D2"
    }
  } else {
    return("ERROR: linkage argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  # Stage 1: Tree clustering
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  n_cols <- ncol(corre)
  # Stage 2: Quasi-Diagonalisation
  clustering <- fastcluster::hclust(euclidean_distance , method = linkage, members = NULL)
  if (is.null(clusters)) {
    fun_clus_num <- function(x,k) list(cluster = stats::cutree(fastcluster::hclust(stats::as.dist(x) , method = linkage, members = NULL), k))
    gap <- cluster::clusGap(as.matrix(euclidean_distance), FUN = fun_clus_num, K.max = floor(n_cols/2), B = 100)
    n_clusters <- cluster::maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = "Tibs2001SEmax")
    n_clusters <- max(2, n_clusters)
  } else{
    n_clusters <- clusters
  }
  elements_in_cluster <- matrix(stats::cutree(clustering, 2:n_clusters), ncol = n_clusters - 1);
  risk_in_cluster <- matrix(0, ncol = n_clusters, nrow = n_cols)
  w <- rep(NA, n_cols)

  assets_in_cluster <- elements_in_cluster[, n_clusters - 1]
  for (i in 1:n_clusters) {
    index_assets_in_cluster <- which(assets_in_cluster == i)
    if (length(index_assets_in_cluster) == 1) {
      w[index_assets_in_cluster] <- 1
      risk_in_cluster[index_assets_in_cluster[1], n_clusters] <- w[index_assets_in_cluster] * covar[index_assets_in_cluster, index_assets_in_cluster] * w[index_assets_in_cluster]
    } else {
      w[index_assets_in_cluster] <- (1/diag(covar[index_assets_in_cluster, index_assets_in_cluster])) / sum(1/diag(covar[index_assets_in_cluster, index_assets_in_cluster]))
      risk_in_cluster[index_assets_in_cluster[1], n_clusters] <- w[index_assets_in_cluster] %*% covar[index_assets_in_cluster, index_assets_in_cluster] %*% w[index_assets_in_cluster]
    }
  }

  for (i in (n_clusters - 1):min(2, (n_clusters - 1))) {
    risk_in_cluster[, i] = risk_in_cluster[, i + 1]
    for (j in 1:i) {
      aux2_assets_in_cluster <- which(elements_in_cluster[, i - 1] == j)
      if (length(unique(elements_in_cluster[aux2_assets_in_cluster, i])) > 1) {
        risk_in_cluster[aux2_assets_in_cluster, i] <- 0
        risk_in_cluster[aux2_assets_in_cluster[1], i] <- sum(risk_in_cluster[aux2_assets_in_cluster, i + 1])
      }
    }
  }
  risk_in_cluster <- matrix(risk_in_cluster[, -1], ncol = n_clusters - 1)
  indexa <- elements_in_cluster[,1] == 1
  indexb <- elements_in_cluster[,1] == 2
  alpha <- matrix(1, ncol = n_clusters, nrow = nrow(covar))
  alpha[, n_clusters] <- 0
  alpha[indexa, 1] = 1 - sum(risk_in_cluster[indexa, 1])/sum(risk_in_cluster[, 1])
  alpha[indexb, 1] = 1 - sum(risk_in_cluster[indexb, 1])/sum(risk_in_cluster[, 1])

  if (n_clusters > 2) {
    for (i in 2:(n_clusters - 1)) {
      groups <- unique(elements_in_cluster[, i - 1])
      for (j in 1:i) {
        if (alpha[groups[j], n_clusters] == 0) {
          index <- elements_in_cluster[, i - 1] == groups[j]
          if (length(unique(elements_in_cluster[index, i])) > 1) {
            sub_groups <- unique(elements_in_cluster[index, i])
            indexa <- elements_in_cluster[ , i] == sub_groups[1]
            indexb <- elements_in_cluster[ , i] == sub_groups[2]
            alpha[indexa, i] = 1 - sum(risk_in_cluster[indexa, i])/sum(risk_in_cluster[index, i])
            alpha[indexb, i] = 1 - sum(risk_in_cluster[indexb, i])/sum(risk_in_cluster[index, i])
          } else {
            alpha[index, n_clusters ] <- 1
          }
        }
      }
    }
  }
  alpha <- matrix(alpha[, -n_clusters], ncol = n_clusters - 1)
  weights <- apply(alpha, 1, prod)*w
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - HERC")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}



