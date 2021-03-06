% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.r
\name{dgp_shear}
\alias{dgp_shear}
\title{Generate sheared bivariate normal data compared to a bivariate normal distribution
with a certain mean and covariance matrix.}
\usage{
dgp_shear(n, mu, Sigma, v = c(1, 1))
}
\arguments{
\item{n}{sample size.}

\item{mu}{expectation, a vector of size two.}

\item{Sigma}{covariance matrix.}

\item{v}{shearing vector.}
}
\value{
a \eqn{n\times 2} matrix of the data.
}
\description{
Generate sheared bivariate normal data compared to a bivariate normal distribution
with a certain mean and covariance matrix.
}
\author{
Nadja Klein.
}
