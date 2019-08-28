{}

#' Generate translated bivariate normal data (with a certain mean and covariance matrix) compared to a zero mean distribution.
#' 
#' 
#' @param n sample size.
#' @param mu expectation, a vector of size two.
#' @param Sigma covariance matrix.
#' @param v translation vector.
#' @return a \eqn{n\times 2} matrix of the data. 
#' @author Nadja Klein.
#' @import mvtnorm
#' @export

dgp_trans <- function(n, mu, Sigma, v=c(0,0))
  {
  X <- rmvnorm(n, mu, Sigma)

	v <- matrix(rep(v, n), ncol=2, byrow=TRUE)
	X <- X + v
	return(X)
	}

	
#' Generate scaled bivariate normal data compared to a bivariate normal distribution
#' with unit marginal variances.
#' 
#' 
#' @param n sample size.
#' @param mu expectation, a vector of size two.
#' @param Sigma covariance matrix.
#' @param v scaling vector.
#' @return a \eqn{n\times 2} matrix of the data. 
#' @author Nadja Klein.
#' @import mvtnorm
#' @export

dgp_scale <- function(n, mu, Sigma, v=c(1,1))
  {
  X <- rmvnorm(n, mu, Sigma)

	v <- matrix(rep(v, n), ncol=2, byrow=TRUE)
	X <- X * v
	return(X)
	}

#' Generate sheared bivariate normal data compared to a bivariate normal distribution
#' with a certain mean and covariance matrix.
#' 
#' 
#' @param n sample size.
#' @param mu expectation, a vector of size two.
#' @param Sigma covariance matrix.
#' @param v shearing vector.
#' @return a \eqn{n\times 2} matrix of the data. 
#' @author Nadja Klein.
#' @import mvtnorm
#' @export

dgp_shear <- function(n, mu, Sigma, v=c(1,1))
  {
  X <- rmvnorm(n, mu, Sigma)

	v <- matrix(v, ncol=2, byrow=TRUE)
  X <- t(apply(X, 1, function(x) t(v) %*% x))
	return(X)
	}
	

#' Generate rotated bivariate normal data compared to a bivariate normal distribution
#' with a certain mean and covariance matrix.
#' 
#' 
#' @param n sample size.
#' @param mu expectation, a vector of size two.
#' @param Sigma covariance matrix.
#' @param gamma angle of rotation in radients.
#' @return a \eqn{n\times 2} matrix of the data. 
#' @author Nadja Klein.
#' @import mvtnorm
#' @export

dgp_rot <- function(n, mu, Sigma, gamma=pi/4)
  {
  X <- rmvnorm(n, mu, Sigma)

	v <- matrix(c(cos(gamma), -sin(gamma), sin(gamma), cos(gamma)), ncol=2, byrow=TRUE)
  X <- t(apply(X, 1, function(x) t(v) %*% x))
	return(X)
	}
	

#' Generate bivariate data with different margins and dependence structures from Archimedian copulas.
#' 
#' 
#' @param n sample size.
#' @param family the copula of type Archimedian.
#' @param margins the margins to be specified.
#' @param paramMargins parameters of marginal distributions.
#' @param rho the copula parameter
#' @param ... further arguments for function \code{\link{mvdc}}.
#' @return a \eqn{n x 2} matrix of the data which is an object of class \code{\link{mvdc}}.
#' @seealso \code{\link{rMvdc}} and \code{\link{mvdc}} in \code{copula-package} for details. 
#' @author Nadja Klein.
#' @import copula
#' @export

dgp_cop <- function(n, family="clayton", margins=c("norm", "norm"),
										paramMargins=list(
																	list(mean = 0, sd = 1), 
																	list(mean = 0, sd = 1)
																	), rho=2, ...)
  {
	dc <- mvdc(copula = archmCopula(family = family, param = rho),
              margins = margins,
              paramMargins = paramMargins, ...)
  X <- rMvdc(n,dc)

	return(X)
	}
	
