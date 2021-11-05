#' Hierarchical Risk Parity
#'
#' Performs the Hierarchical Risk Parity portfolio proposed strategy by De Prado (2016). Several linkage methods for the hierarchical clustering can be used, by default the "single" linkage is used.
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into distance matrix.
#' @param clustering.method Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "single".
#' @param graph To plot de dendogram set this value to TRUE. By default this value is equal to FALSE.
#' @return portfolio weights
#' @seealso \code{HCAA_porfolio} and \code{HERC_Portfolio}
#' @references De Prado, Marcos Lopez. "Building diversified portfolios that outperform out of sample." The Journal of Portfolio Management 42.4 (2016): 59-69.
#' @aliases HRP
#' @keywords HRP
#' @examples
#' covar <- cov(daily_returns)
#' HRP_Portfolio(covar)
#' @export
HRP_Portfolio = function(covar, clustering.method= "single", graph = FALSE) {
  if (clustering.method %in% c("single", "complete", "average", "ward")) {
    if (clustering.method == "ward") {
      clustering.method = "ward.D2"
    }
  } else {
    return("ERROR: clustering.method argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  # Stage 1: Tree clustering
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  # Stage 2: Quasi-Diagonalisation
  clustering <- fastcluster::hclust(euclidean_distance , method = clustering.method, members = NULL)
  clusters_order <- clustering$order
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
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogam - HRP")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}
