% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCAA_Portfolio.R
\name{HCAA_Portfolio}
\alias{HCAA_Portfolio}
\title{Hierarchical Clutering-Bases Asset Allocation}
\usage{
HCAA_Portfolio(covar, clustering.method = "ward")
}
\arguments{
\item{covar}{Covariance matrix of returns. The covariance matrix will be transformed into correlation matrix and then into distance matrix.}

\item{clustering.method}{Linkage method used in the hierarchical clustering. Allowed options are `single`, `complete`, `average` or `ward`. Default option is `single`.}
}
\value{
portfolio weights
}
\description{
Performs the Hierarchical Clustering-Bases Asset Allocation portfolio strategy proposed by Raffinot (2017). Several linkage methods for the hierarchical clustering can be used, by default the `single` linkage is used.
}
\seealso{
`HRP_porfolio` and `HERC_Portfolio`
}