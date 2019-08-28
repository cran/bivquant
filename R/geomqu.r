#' @include utils.r 
{} 

#'  On Affine Equivariant Multivariate Quantiles
#' 
#' This function fits empirical bivariate quantiles as proposed by  Chakraborty, B. (2001) <https://doi.org/10.1023/A:1012478908041>.
#'
#'
#' 
#' @param data must be a (n,d) matrix of observations. 
#' @param u a (m,d) matrix with m directions.
#' @param p defines the p-norm, must be in \eqn{[1,\infty)}.
#' @param alpha must be missing or a vector of length d+1 with distict values from 1 to n.
#' @param keep_optim boolean keep results from optimization in result object, default is FALSE.
#' @param keep_data boolean keep data parameter in result object, default is TRUE.
#' @return an object of class \code{geomqu} with methods \code{\link{plot.geomqu}} and \code{\link{plot.geomqu}}.
#' @author Nadja Klein.
#' @references Chakraborty, B. (2001). On affine equivariant multivariate quantiles. 
#' \emph{Annals of the Institute of Statistical Mathematics}, \bold{53}, 380--. <https://doi.org/10.1023/A:1012478908041>.
#' @export
geomqu <- function(data, u, p, alpha, keep_optim=FALSE, keep_data=TRUE) {

  if(missing(alpha)) {
    alpha <- find_good_alpha(data)$alpha
  }

  if (!is.matrix(u)) {
    u <- t(as.matrix(u))
  }

  result <- apply(u, 1, geomqu.intern,
                  data=data, p=p, alpha=alpha)
  names(result) <- NULL

  ret_list <- structure(list(), class = "geomqu")
  ret_list$p <- p
  ret_list$theo <- F
  ret_list$alpha <- alpha
	if (keep_data)
    ret_list$data <- data
  if (keep_optim)
    ret_list$optim_results <- result
  qs <- matrix(unlist(lapply(result, function(x) c(x$par, x$convergence))),
               ncol=3, byrow=T)
  ret_list$res <- data.frame(u, qs)
  names(ret_list$res) <- c(paste0("u", 1:ncol(u)),
                           paste0("quantile", 1:ncol(u)),
                           "convergence")
  return(ret_list)
}

# geomqu <- function(
    # data, u, p, a,
    # optim.method=c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG"),
    # optim.keep_res=FALSE,
    # keep_data=TRUE) {

  # optim.method <- match.arg(optim.method)
  # if(missing(a)) {
    # a <- find_good_alpha(data)$a
  # }

  # if (!is.matrix(u)) {
    # u <- t(as.matrix(u))
  # }

  # result <- apply(u, 1, geomqu.intern,
                  # data=data, p=p, alpha=a, optim.method=optim.method)

  # return(result)

  # ret_list <- list()
  # ret_list$u <- u
  # ret_list$record_size <- nrow(data)
  # ret_list$a <- a
  # ret_list$Q <- matrix(unlist(lapply(result, function(x) x$par)),
                       # ncol=length(result[[1]]$par), byrow=T)
  # if (keep_data)
    # ret_list$data <- data
  # if (optim.keep_res)
    # ret_list$optim_results <- result

  # class(ret_list) <- "geomqu"
  # return(ret_list)

# }



#'  On Affine Equivariant Multivariate Quantiles
#' 
#' This function fits empirical bivariate quantiles as proposed by  Chakraborty, B. (2001).
#'  2 dimensional data and with p=2
#'
#' @param data must be a (n,2) matrix of observations. 
#' @param probs vector of probs which is used to calculate the us
#' @param k number of us per prob
#' @param alpha missing or a vector of length 3 with distict values form 1 to n
#' @seealso \code{\link{geomqu}} for details. 
#' @author Nadja Klein.
#' @references Chakraborty, B. (2001). On affine equivariant multivariate quantiles. 
#' \emph{Annals of the Institute of Statistical Mathematics}, \bold{53}, 380--403.
#' @export
#' @examples
#' \dontrun{
# generate data
#' require("MASS")
#' require("mvtnorm")
#' set.seed(42)
#' n <- 50
#' mu <- c(6, 10)
#' rho <- 0.5
#' Sigma <- matrix(c(
#'     1.0, rho,
#'     rho, 1.0
#'   ),
#'   ncol=2, byrow=TRUE)
#' 
#' X <- rmvnorm(n, mu, Sigma)
#' result <- geomqu2d_norm2(X, probs=c(0.8,0.9), k=8)
#' plot(result)
 
#' } 

geomqu2d_norm2 <- function(data, probs, alpha, k=8) {
  # if alpha is missing, find good one
  if(missing(alpha)) {
    alpha <- find_good_alpha(data)$alpha
  }

  # calculate us
  us <- lapply(probs, function(a, n) {
    radius <- 2*a - 1
    if (radius == 0)
      return (matrix(c(0,0), nrow=1))
    else
      ball_2d(radius, n)
  }, n=k)

  if (length(us) == 1) {
    us <- us[[1]]
  } else {
    us <- Reduce(rbind, us[2:length(us)], us[[1]])
  }
  ps <- unlist(lapply(probs, function(p, n) {if (p==.5) p else rep_len(p, n)}, n=k))

  # calculate quantiles
  quantiles <- geomqu(data, us, 2, alpha)

  quantiles$res$prob <- ps
  return(quantiles)
}


#'  On Affine Equivariant Multivariate Quantiles
#' 
#' This function generates the directions u for which the geometric quantiles  are computed. 
#'
#' The Eulicdean norm of u has to be \eqn{\le 1}, see Chakraborty, B. (2001).
#'
#' 
#' @param r length of u computed as \eqn{\sqrt{u_1^2+u_2^2}}.  
#' @param n number of u's, default is 10.
#' @return a matrix of dimension (n,2) containing rowwise the distrinct u's. 
#' @seealso \code{\link{geomqu}}.
#' @author Nadja Klein.
#' @references Chakraborty, B. (2001). On affine equivariant multivariate quantiles. 
#' \emph{Annals of the Institute of Statistical Mathematics}, \bold{53}, 380--403. <https://doi.org/10.1023/A:1012478908041>.

#' @export

ball_2d <- function(r, n=10) {
  alpha <- seq(0, 2*pi-2*pi/n, length.out=n)
  x <- r * cos(alpha)
  y <- r * sin(alpha)
  matrix(c(x, y), ncol=2)
}

#'  On Affine Equivariant Multivariate Quantiles
#' 
#' This function returns the transformation matrix X(a), see Chakraborty, B. (2001). 
#'
#' 
#' @param X must be a (n,2) matrix.  
#' @param eps stopping criterion.
#' @return a list with a the selected rows of the data, X is the X(a) matrix, X_reduced is the data matrix X without the selected columns. 
#' @seealso \code{\link{geomqu}}.
#' @author Nadja Klein, Paul Wiemann.
#' @references Chakraborty, B. (2001). On affine equivariant multivariate quantiles. 
#' \emph{Annals of the Institute of Statistical Mathematics}, \bold{53}, 380--403. <https://doi.org/10.1023/A:1012478908041>.

#' @export

find_good_alpha <- function(X, eps=0.025) {
  Sigma_hat <- cov(X)
  Sigma_hat_inv <- solve(Sigma_hat)
  candidates <- combn(nrow(X), ncol(X)+1)
  ratios <- rep_len(Inf, ncol(candidates))
  for (i in 1:ncol(candidates)) {
    candidate <- candidates[,i]
    X_alpha.c <- X_alpha(X, candidate)
    eigen_values <- eigen(t(X_alpha.c) %*% Sigma_hat_inv %*% X_alpha.c)$values
    geom_mean.c <- exp(mean(log(eigen_values)))
    mean.c <- mean(eigen_values)
    ratios[i] <- mean.c / geom_mean.c
    if (ratios[i] < 1 + eps)
      break
  }
  if (i == ncol(candidates))
    j <- which.min(ratios)
  else
    j <- i
  return(list(
    alpha = candidates[,j],
    X = X_alpha(X, candidates[,j]),
    X_reduced = X[-candidates[,j],]
  ))
}
