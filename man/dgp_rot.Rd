% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.r
\name{dgp_rot}
\alias{dgp_rot}
\title{Generate rotated bivariate normal data compared to a bivariate normal distribution
with a certain mean and covariance matrix.}
\usage{
dgp_rot(n, mu, Sigma, gamma = pi/4)
}
\arguments{
\item{n}{sample size.}

\item{mu}{expectation, a vector of size two.}

\item{Sigma}{covariance matrix.}

\item{gamma}{angle of rotation in radients.}
}
\value{
a \eqn{n\times 2} matrix of the data.
}
\description{
Generate rotated bivariate normal data compared to a bivariate normal distribution
with a certain mean and covariance matrix.
}
\author{
Nadja Klein.
}
