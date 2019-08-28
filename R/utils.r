empcdf <- function(q, y) {
	# if(is.vector(q) == 1) {
		# q <- matrix(q, nrow = 1)
	# }
	# if(NCOL(q) != NCOL(y)) {
		# stop("please specify at least one quantile for each marginal")
	# }
	# func <- apply(y, FUN = ecdf, MARGIN = 2)
	# ret <- sapply(1:NCOL(q), FUN = function(x) {func[[x]](q[,x])})
	# return(ret)
	arg <- q
	dat <- y
	ret <- myecdf_mult(arg, dat)
	return(ret)
}

empqf <- function(p, y) {
	# if(is.vector(p) == 1) {
		# p <- matrix(p, nrow = 1)
	# }
	# if(NCOL(p) != NCOL(y)) {
		# stop("please specify at least one quantile for each marginal")
	# }
	# ret <- p
	# for(i in 1:NCOL(p)) {
		# ret[, i] <- sapply(1:NROW(p), FUN = function(x) {emp.qu(p[x,i], y[,i])})
	# }
	arg <- p
	dat <- y
	ret <- myempqu_mult(arg, dat)
	return(ret)
}

# univariate case first

# new implementation of CDF
myecdf <- function(arg, dat=NULL)
  {
  if(is.null(dat))
    dat <- arg
#  tab <- cumsum(table(dat))/length(dat)
#  ff <- function(x) tab[max(which(as.numeric(names(tab))<=x))]
#  datsort <- sort(dat)
  ff <- function(x) length(which(dat<=x))/length(dat)
  as.numeric(sapply(arg, ff))
  }

# new implementation of empirical quantile function
myempqu <- function(arg, dat)
  {
#  tab <- cumsum(table(dat))/length(dat)
#  ff <- function(x) as.numeric(names(tab)[min(which(tab>=x))])
  tab <- myecdf(dat)
  ff <- function(x) min(dat[which(tab>=x)])
  sapply(arg, ff)
  }
	
# and now the multivariate case
myecdf_mult <- function(arg, dat)
  {
  if(is.vector(dat))
    stop("dat should be a matrix")
  if(is.vector(arg) && length(arg==ncol(dat)))
    arg <- matrix(arg, nrow=1)
  res <- arg
  for(i in 1:ncol(arg))
    res[,i] <- myecdf(arg[,i], dat[,i])
  res
  }

myempqu_mult <- function(arg, dat)
  {
  if(is.vector(dat))
    stop("dat should be a matrix")
  if(is.vector(arg) && length(arg==ncol(dat)))
    arg <- matrix(arg, nrow=1)
  res <- arg
  for(i in 1:ncol(arg))
    res[,i] <- myempqu(arg[,i], dat[,i])
  res
  }


  
emp.qu <- function (arg, dat) 
  {
	d <- length(arg)
	if (d > 1) 
	  {
		val <- matrix(0, d, 1)
		n <- dim(dat)[1]
		for (kk in 1:d) 
		  {
			or <- rank(dat[, kk])
			uni <- or/(n + 1)
			cateca <- c(uni, arg[kk])
			or2 <- rank(cateca)
			inter <- or2[length(cateca)]
			if (inter == (n + 1)) 
					inter <- n
			so <- sort(dat[, kk])
			val[kk] <- so[inter]
		  }
	  } else {
		n <- length(dat)
		or <- rank(dat)
		uni <- or/(n + 1)
		cateca <- c(uni, arg)
		or2 <- rank(cateca)
		inter <- or2[length(cateca)]
		if (inter == (n + 1)) 
				inter <- n
		val <- sort(dat)[inter]
	  }
	return(val)
}



#===================================================================================================================
#Chakraborty, 2001


lp_norm <- function(x, p) {
  x_abs <- abs(x)
  if (p == Inf)
    return(max(x_abs))
  else
    return(sum(x_abs^p)^(1/p))
}

# phi as definded at (2.1)
# u, s vector
# p \in N
phi <- function(u, s, p) {
  lp_norm(s, p) + t(u)%*%s
}

# as above but u must be an 2x1 matrix
phi_transu <- function(u, s, p) {
  lp_norm(s, p) + u%*%s
}

# d phi/ds
dphi_transu <- function(u, s, p) {
  if (all(s==0))
    return(u)
  sign(s)*abs(s)^(p-1)/lp_norm(s, p)^(p-1) + u
}

# computes X_alpha for a given alpha as describe on page 382 last paragraph
X_alpha <- function(data, alpha) {
  d <- length(alpha)-1
  apply(data[alpha[2:(d+1)],], 1, function(x) x - data[alpha[1],])
}

# described on page 383 formula after first paragraph
v <- function(u, X_a_inv, q) {
  if (all(u == 0))
    return(u)
  u_ <- X_a_inv %*% u
  return(u_ * lp_norm(u, q) / lp_norm(u_, q))
}

# function to be minimised (def at (2.3))
f <- function(Q, X_a_inv, X_reduced, u, pp) {
  X_ <- apply(X_reduced, 1, function(x) x - Q)
  X__ <- apply(X_, 2, function(x) X_a_inv %*% x)
  sum(apply(X__, 2, function(x) {phi_transu(u, x, pp)}))
}

# df/dQ
df <- function(Q, X_a_inv, X_reduced, u, pp) {
  X_ <- apply(X_reduced, 1, function(x) x - Q)
  X__ <- apply(X_, 2, function(x) X_a_inv %*% x)
  apply(apply(X__, 2, function(x) {dphi_transu(u, x, pp)}), 1, sum) %*% (-1 * X_a_inv)
}


# ableitung von f nach Q
dQf <- function(Q, X_a_inv, X_reduced, u, pp) {
  X_ <- apply(X_reduced, 1, function(x) x - Q)
  X__ <- apply(X_, 2, function(x) X_a_inv %*% x)
  apply(apply(X__, 2, function(x) {dphi_transu(u, x, pp)}), 1, sum) %*% (-1 * X_a_inv)
}



# calculates chakraborty quantile for one u
geomqu.intern <- function(u, data, p, alpha) {
  if (p == 1)
    q = Inf
  else
    q = p/(p-1)
  u <- as.matrix(u)

  X_a <- X_alpha(data, alpha)
  X_a_inv <- solve(X_a)

  X_reduced <- data[-alpha,]
  u_ <- v(u, X_a_inv, q)

  res <- optim(apply(data, 2, median), f, dQf,
               X_a_inv=X_a_inv, X_reduced=X_reduced, u=t(u_), pp=p)
  if (res$convergence != 0)
    warning("optimization did not converge")
  res
}

polygon.geomqu <- function(x) {
  if (!inherits(x, "geomqu"))
    stop(gettextf("'x' must inherit from class %s", dQuote("geomqu")),
         domain = NA)
  stopifnot(ncol(x$Q)==2)
  polygon(x$Q)
}


#====================================================================================
#helper functions for theoretical quantiles of a multivariate normal distribution

dQintegrand <- function(x, Q, U, mu, sigma) {
  a <- as.numeric(exp(-1/2 * t(x-mu) %*% solve(sigma) %*% (x-mu)))
  b <- (2*pi)^(length(x)/2) * sqrt(det(sigma))
  c <- (Q - x)/euclidean_norm(x - Q) - U
  return(a/b * c)
}

integrand <- function(x, Q, U, mu, sigma) {
  a <- as.numeric(exp(-1/2 * t(x-mu) %*% solve(sigma) %*% (x-mu)))
  b <- (2*pi)^(length(x)/2) * sqrt(det(sigma))
  c <- as.numeric(euclidean_norm(x - Q) + t(U) %*% (x - Q))
  return(a/b * c)
}

integral <- function(Q, u, mu, sigma, ...) {
  res <- cubature::adaptIntegrate(integrand,
                                  Q=Q, U=u, mu=mu, sigma=sigma, ...)
  res$integral
}

dQintegral <- function(Q, u, mu, sigma, ...) {
  res <- cubature::adaptIntegrate(dQintegrand,
                                  Q=Q, U=u, mu=mu, sigma=sigma, ...)
  res$integral
}


# calculates the euclidean norm
euclidean_norm <- function(x) {
  return(as.numeric(sqrt(t(x) %*% x)))
}

# density of the multivariate normal distribution
mnorm_dens <- function(x, mu=c(0,0), sigma=diag(c(1,1))) {
  sigma <- as.matrix(sigma)
  a <- (2*pi)^(length(mu)/2) * sqrt(det(sigma))
  b <- as.numeric(exp(
    -1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu)
  ))
  return(b/a)
}

# sucht einen parameter bei dem f(par) > 0 ist.
#
# es muss gelten:
# - f >= 0
# - f hat einen kompakten Traeger
# - es existert ein x in [lower, upper] mit f(x) > 0
find_non_zero <- function(f, lowerLimit, upperLimit, tol, ...) {
  limits <- list(
    list(l=lowerLimit, u=upperLimit)
  )
  c <- 0
  while (length(limits) > 0) {
    c <- c + 1
    l <- limits[[1]]$l
    u <- limits[[1]]$u
    limits <- tail(limits, n=length(limits)-1)
    m <- (l + u)/2
    fm <- f(m, ...)
    if (fm > 0) {
      return(fm)
    } else {
      if (any(u-m > tol))
        limits[[length(limits) + 1]] <- list(l=m, u=u)
      if (any(m-l > tol))
        limits[[length(limits) + 1]] <- list(l=l, u=m)
      stopifnot(length(limits)>0)
    }
  }
}

# versucht den Traeger einer Funktion zu bestimmten (+- tol)
#
# es gelten die gleichen Einschraenkungen wie bei find_non_zero
# anstelle von lower und upper limit kann auch start uebergeben werden.
# es muss dann gelten: f(start) > 0
find_limits <- function(f, lowerLimit, upperLimit, start=NA, tol=10, ...) {
  if (any(is.na(start))) {
    start <- find_non_zero(f, lowerLimit, upperLimit, tol, ...)
  }
  stopifnot(any(f(start, ...) > 0))
  # upper
  upper <- start
  repeat{
    if (all((f(upper + tol, ...)) == 0)) {
      break
    }
    upper <- upper + tol
  }

  # upper
  lower <- start
  repeat{
    if (all(f(lower - tol, ...) == 0)) {
      break
    }
    lower <- lower - tol
  }

  return (list(
    lower=lower,
    upper=upper,
    start=start,
    tol=tol
  ))
}


