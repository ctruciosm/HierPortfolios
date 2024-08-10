#' Constrained Hierarchical Risk Parity
#'
#' Performs the Constrained Hierarchical Risk Parity portfolio strategy proposed by Pfitzinger and Katzke (2019).
#'
#' @param covar Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into a distance matrix.
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @param tau Parameter to evaluate asset similarity at the cluster edges. Default value is 1.
#' @param UB Upper bound for weights. By default this value is equal to NULL
#' @param LB Lower bound for weights. By default this value is equal to NULL
#' @return portfolio weights
#' @seealso \code{HCAA_Portfolio}, \code{HRP_Portfolio} and \code{HERC_Portfolio}
#' @references Pfitzinger, J., and Katzke, N. A constrained hierarchical risk parity algorithm with cluster-based capital allocation (2019). Working Paper.
#' @aliases DHRP_Portfolio
#' @keywords DHRP
#' @examples
#' covar <- cov(mldp_returns)
#' DHRP_Portfolio(covar)
#' @author Carlos Trucios and Moon Jun Kwon
#' @export
DHRP_Portfolio = function(covar, graph = FALSE, tau = 1, UB = NULL, LB = NULL) {

  if (is.null(UB)) UB <- rep(1, nrow(covar))
  if (is.null(LB)) LB <- rep(0, nrow(covar))
  # Stage 1: Tree clustering
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  # Stage 2: Divisive clustering
  distmat <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
  clustering <- cluster::diana(distmat)
  clusters_order <- clustering$order
  clusters_height <- clustering$height
  # Stage 3: Recursive Bisection
  weights <- rep(1,ncol(covar))
  index <- list(clusters_order)
  while (length(index) > 0) {
    new_index <- list()
    for (i in index) {
      aux <- floor(length(i)/2 - (length(i)/2 - 1) * tau):floor(length(i)/2 + (length(i)/2 - 1) * tau)
      index_aux <- clusters_order %in% i
      split_tau <- aux[which.max(clusters_height[index_aux][aux])]
      indexa <- i[1:split_tau]
      indexb <- i[-c(1:split_tau)]
      covar_clustera <- as.matrix(covar[indexa, indexa])
      covar_clusterb <- as.matrix(covar[indexb, indexb])
      weightsa <- 1/diag(covar_clustera) / sum(1/diag(covar_clustera))
      weightsb <- 1/diag(covar_clusterb) / sum(1/diag(covar_clusterb))
      variance_clustera <- weightsa %*% covar_clustera %*% weightsa
      variance_clusterb <- weightsb %*% covar_clusterb %*% weightsb
      alpha <- as.numeric(1 - variance_clustera/(variance_clustera + variance_clusterb))
      v_alpha <- c(alpha, 1 - alpha)

      LBsub <- c(sum(LB[indexa]), sum(LB[indexb])) / c(prod(weights[indexa]), prod(weights[indexb]))
      UBsub <- c(sum(UB[indexa]), sum(UB[indexb])) / c(prod(weights[indexa]), prod(weights[indexb]))
      maxit <- 100
      niter <- 0
      while (any(v_alpha > UBsub | v_alpha < LBsub) && niter < maxit) {
        alpha_tilde <-  sapply(sapply(v_alpha, min, UBsub), max, LBsub)
        aux <- which(alpha_tilde != UBsub & alpha_tilde != LBsub)
        alpha_tilde[aux] <- alpha_tilde[aux] + (1 - sum(alpha_tilde)) * alpha_tilde[aux] / sum(alpha_tilde[aux])
        alpha_tilde <- alpha_tilde / sum(alpha_tilde)
        v_alpha <- alpha_tilde
        niter <- niter + 1
      }

      weights[indexa] <- weights[indexa] * v_alpha[1]
      weights[indexb] <- weights[indexb] * v_alpha[2]
      if (length(indexa) > 1) new_index <- c(new_index, list(indexa))
      if (length(indexb) > 1) new_index <- c(new_index, list(indexb))
      index = new_index
    }

  }
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - DHRP")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}





