
#' @include geomqu.r 
{}


#' Plot the estimated bivariate quantiles (for the  CDF (cumulative distrubtion function)-based quantiles) and the data for a fitted bivquant object. See ?bivquant for more details.
#' 
#' 
#' @param x a fitted \code{bivquant} model.
#' @param mains title of figures.
#' @param ... further arguments for function \code{\link{plot}}.
#' @seealso \code{\link{plot}} for details. 
#' @S3method plot bivquant
#' @method plot bivquant
#' @author Nadja Klein.


plot.bivquant <- function(x, 
													mains=c("CDF-based quantiles (original scale)", "CDF-based quantile (unit square)"), 
													...) 
	{
	transformed <- x$transformed
	if(transformed==FALSE)
    {
		plot(x$y, main=mains[1], ...)
	  tauseq <- x$tau
	  for(l in 1:length(tauseq))
	    lines(x$bivqu[[l]][,1],x$bivqu[[l]][,2],lwd=2)
	  } else { 
		opar <- par(mfrow=c(1,2),no.readonly =TRUE)
		on.exit(par(opar))
		plot(x$y, main=mains[1], ...)
	  tauseq <- x$tau
	  for(l in 1:length(tauseq))
	    lines(x$bivqu[[l]][,1],x$bivqu[[l]][,2],lwd=2)
    plot(x$tildey, main=mains[2], ...)
	  for(l in 1:length(tauseq))
	    lines(x$tildebivqu[[l]][,1],x$tildebivqu[[l]][,2],lwd=2)
		}
	}



#' Plot the estimated bivariate quantiles (for the geometric quantiles of Chakraborty, B. (2001).) and the data. See ?geomqu for more details.
#' 
#' 
#' @param x a fitted \code{geomqu} model.
#' @param ... further arguments for function \code{\link{plot}}.
#' @seealso \code{\link{plot}} for details. 
#' @S3method plot geomqu
#' @method plot geomqu
#' @author Nadja Klein.


plot.geomqu <- function(x, ...) 
	{
  stopifnot(length(grep("quantile", names(x$res)))==2)
  plot_args <- list(xlab="", ylab="", xlim=range(x$res$quantile1),
                ylim=range(x$res$quantile2), main="Geometric Quantiles")
  #plot_args_ellipsis <- list(...)
  #plot_args[names(plot_args_ellipsis)] <- plot_args_ellipsis
  #do.call(plot, c(list(x=1, type="n"), plot_args))
  plot(x$data, col="grey",pch=19)
  if ("prob" %in% names(x$res)) 
		{
    for (prob in unique(x$res$prob)) 
			{
			df <- x$res[x$res$prob == prob,]
      if (prob == .5) 
				{
        polygon(df$quantile1, df$quantile2, pch=19)
				}
      polygon(df$quantile1, df$quantile2, pch=19)
			}
		} else {
		df <- x$res
    polygon(df$quantile1, df$quantile2,pch=19)
		}
	}



