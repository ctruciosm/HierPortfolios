#' Hierarchical Risk Parity
#'
#' Performs the Hierarchical Risk Parity portfolio proposed strategy by De Prado (2016). Several linkage methods for the hierarchical clustering can be used, by default the "single" linkage is used.
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
#' @param linkage Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "single".
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @return portfolio weights
#' @seealso \code{HCAA_Portfolio}, \code{HERC_Portfolio} and \code{DHRP_Portfolio}
#' @references De Prado, Marcos Lopez. "Building diversified portfolios that outperform out of sample." The Journal of Portfolio Management 42.4 (2016): 59-69.
#' @aliases HRP_Portfolio
#' @keywords HRP
#' @examples
#' covar <- cov(mldp_returns)
#' HRP_Portfolio(covar)
#' @author Carlos Trucios
#' @export
HRP_Portfolio = function(covar, linkage = "single", graph = FALSE) {
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
  # Stage 2: Quasi-Diagonalisation
  clustering <- fastcluster::hclust(euclidean_distance , method = linkage, members = NULL)
  clusters_order <- clustering$order
  # Stage 3: Recursive Bisection
  weights <- rep(1,ncol(covar))
  index <- list(clusters_order)
  while (length(index) > 0) {
    new_index <- list()
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
      if (length(indexa) > 1) new_index <- c(new_index, list(indexa))
      if (length(indexb) > 1) new_index <- c(new_index, list(indexb))
      index = new_index
    }
  }
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - HRP")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}
