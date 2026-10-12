#' Hierarchical Risk Budgeting
#'
#' Computes portfolio weights using the Hierarchical Risk Budgeting (HRB) strategy proposed by Koumou (2024). The method combines hierarchical clustering with a risk-budgeting procedure. The clustering is based on asset correlations, and several linkage methods are available. By default, single linkage is used.
#'
#' @param covar Covariance matrix of returns. The covariance matrix is converted into a correlation matrix, which is then used to construct the distance matrix for hierarchical clustering.
#' @param linkage Linkage method used in the hierarchical clustering. Allowed options are "single", "complete", "average" or "ward". Default option is "single".
#' @param graph To plot de dendrogram set this value to TRUE. By default this value is equal to FALSE.
#' @param upsilon Pparameter controlling the sensitivity of the allocation procedure to rare assets. The default is \code{2}, as in the empirical application of Koumou (2024)
#' @param Sigma If "full", the portfolio volatility  is used as R(w), if "diag" the measure used in the empirical application of the original paper, which ignores dependence between assets is used. By default "diag" is used.
#' @return portfolio weights
#' @seealso \code{HRP_Portfolio}, \code{HCAA_Portfolio}, \code{HERC_Portfolio} and \code{DHRP_Portfolio}
#' @references Koumou, Gilles Boevi. "Hierarchical Risk Budgeting". Journal of Financial Data Science 6.2 (2024): 35-53
#' @aliases HRB_Portfolio
#' @keywords HRB
#' @examples
#' covar <- cov(mldp_returns)
#' HRB_Portfolio(covar)
#' @author Carlos Trucios
#' @export
HRB_Portfolio = function(covar, linkage = "single", graph = FALSE, upsilon = 2, Sigma = c("diag", "full")) {
  if (linkage %in% c("single", "complete", "average", "ward")) {
    if (linkage == "ward") {
      linkage = "ward.D2"
    }
  } else {
    return("ERROR: linkage argument only supports 'single', 'complete', 'average' or 'ward' options")
  }
  Sigma <- match.arg(Sigma)

  n <- ncol(covar)
  corre <- stats::cov2cor(covar)
  distance <- sqrt(0.5 * (1 - corre))
  euclidean_distance <- stats::dist(distance, method = "euclidean", diag = TRUE, upper = TRUE, p = 2) / n
  clustering <- fastcluster::hclust(euclidean_distance, method = linkage, members = NULL)
  distance_bar <- as.matrix(stats::cophenetic(clustering))
  S_bar <- 1 -  distance_bar

  if (Sigma == "diag") {
    if (isTRUE(all.equal(upsilon, 2))) {
      b <- CVXR::Variable(n)
      objective <- CVXR::Minimize(CVXR::quad_form(b, S_bar))
      constraints <- list(b >= 0, sum(b) == 1)
      problem <- CVXR::Problem(objective, constraints)
      fit <- CVXR::psolve(problem)
      b_star <- CVXR::value(b)
      sigma_hat <- sqrt(diag(covar))
      w <- b_star / sigma_hat
      weights <- w / sum(w)
    } else {
      obj_function <- function(z, covar, S_bar, upsilon) {
        w <- exp(z) / sum(exp(z))
        q <- w * sqrt(diag(covar)) / sum(w * sqrt(diag(covar)))
        xi <- drop(S_bar %*% q)
        H <- sum(q * xi^(upsilon - 1))^(1 / (1 - upsilon))
        if (!is.finite(H)) return(1e10)
        -H
      }
      fit <- stats::optim(par = rep(0, n), fn = obj_function, covar = covar, S_bar = S_bar, upsilon = upsilon, method = "BFGS", control = list(maxit = 5000))
      z <- fit$par
      weights  <- exp(z) / sum(exp(z))
    }
  } else {
    obj_function <- function(z, covar, S_bar, upsilon) {
      w <- exp(z) / sum(exp(z))
      Rw <- sqrt(drop(t(w) %*% covar %*% w))
      RC <- w * drop(covar %*% w) / Rw
      q <- RC / Rw
      xi <- drop(S_bar %*% q)
      H <- sum(q * xi^(upsilon - 1))^(1 / (1 - upsilon))
      if (!is.finite(H)) return(1e10)
      -H
    }
    fit <- stats::optim(par = rep(0, n), fn = obj_function, covar = covar, S_bar = S_bar, upsilon = upsilon, method = "BFGS", control = list(maxit = 5000))
    z <- fit$par
    weights  <- exp(z) / sum(exp(z))
  }
  if (graph) plot(clustering, xlab = "", ylab = "", main = "Cluster Dendrogram - HRP")
  weights <- data.frame(weights)
  row.names(weights) <- colnames(covar)
  return(weights)
}
