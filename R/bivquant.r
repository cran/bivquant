#' @include utils_lp.r 
#' @include utils.r 
{} 


#'  Bivariate Quantiles
#' 
#' This function fits the empirical bivariate quantiles based on the CDF (cumulative distrubtion function). We use 
#' linear programming. Currently, the solver is \code{\link{lp}} from the package
#' \code{lpSolve}
#'
#'
#' 
#' @param y the responses in a matrix or data frame with 2 columns and rows equal to the number of observations.  
#' @param alphaseq The angles along which the quantile should be computed, can be a vector. If not specified, quantiles will
#'				will be computed for a equidistant grid from 0 to \eqn{\pi/2} of length 10 is used.
#' @param tau The quantile level. If not specified, the median, \eqn{\tau=0.5} will be computed.
#' @param transformed Default is FALSE specifying that quantiles on the original scale are returned. If TRUE, quantiles on the
#'				unit square are returned in addition.
#' @return an object of class \code{bivquant}. 
#' @seealso \code{\link{lp}}.
#' @author Nadja Klein.
#' @details The function imitates rotation around (1,1) in the transformed
#'	coordinate system and thus allows to estimate the marginal quantiles.
#' @references Nadja Klein and Thomas Kneib (2019). Directional Bivariate Quantiles - A Robust Approach based on the Cumulative Distribution Function. To appear in Advances in Statistical Analysis (AStA)
#' @importFrom regpro emp.quantile
#' @importFrom lpSolve lp
#' @export
#' @examples
#' \donttest{
# generate data
#' set.seed(123)
#' n <- 200
#' tauseq <- seq(0.1,0.9,by=0.1)
#' alphas <- seq(0*pi/32,16*pi/32,by=0.5*pi/32)
#'
#' X <- dgp_cop(n, family="clayton", margins=c("norm", "norm"),
#'			 paramMargins=list(list(mean = 4, sd = 1), list(mean = 4, sd = 5)),
#'			 rho=1.75)
#'
#'
# calc & plot quantiles
#' bivqu <- bivquant(X,alpha=alphas,tau=tauseq)
#' plot(bivqu, pch=20,col="grey")
 
#' } 

bivquant <- function(y, alphaseq=NULL, tau=NULL, transformed=FALSE)
  {
	
	if((!is.data.frame(y)) & (!is.matrix(y)))
		stop("y has two be a data frame or matrix with two columns")
	
	n <- dim(y)[1]
	
	if(dim(y)[2]!=2)
		stop("response has to be two dimensional")
		
  yold <- y
	if(any(is.na(y)))
	  {
	  y <- na.omit(y)
		ndiff <- n-dim(y)[1]
		n <- n-ndiff
		warning(paste("Data contain missings, ", ndiff, " observation(s) removed", sep =""))
		}
		
  y <- empcdf(y, y)-1/(2*n)
	
	if(is.null(alphaseq))
	  alphaseq <- seq(0,pi/2,10)
	if(is.null(tau))
	  tau <- 0.5
	if(!transformed%in%c(TRUE,FALSE))
	  stop("transformed has to be either TRUE or FALSE (default).")	
  
  get_biv <- function(y, alphai, taui, transformed=transformed) 
    { 		
		#define subsets
		indexe <- GetIndex_new(y, alphai)
		i2 <- which(indexe$I2==TRUE)
		i1 <- which(indexe$I1==TRUE)
		yorig <- y
		yindexe <- matrix(y[c(i1,i2), ],ncol=2)

		li1 <- length(i1)
		li2 <- length(i2)
		
		if((li1>0) & (li2>0))
		  {
			yi1 <- matrix(yindexe[1:li1,],ncol=2)
			yi2 <- matrix(yindexe[(li1+1):(li1+li2),],ncol=2)
			cvec <- c(0, rep(1-taui, li1) ,rep(taui,  li1), rep(1-taui, li2) ,rep(taui,  li2))
		
			y1 <- (1-yi1[,1])/cos(alphai)
			y2 <- (1-yi2[,2])/sin(alphai)
			b <- c(y1,y2)
  
			A <- cbind(rep(1, n), 
								rbind(diag(li1), matrix(0,ncol=li1,nrow=li2)), rbind(-diag(li1), matrix(0,ncol=li1,nrow=li2)),
								rbind(matrix(0,ncol=li2,nrow=li1),diag(li2)), rbind(matrix(0,ncol=li2,nrow=li1),-diag(li2)))
		
			} else if(li1==0){
			yi2 <- matrix(yindexe[1:li2,],ncol=2)
			cvec <- c(0, rep(1-taui, li2) ,rep(taui,  li2))
		
			y2 <- (1-yi2[,1])/sin(alphai)
			b <- c(y2)
  
			A <- cbind(rep(1, li2), 
								diag(li2),-diag(li2))
			}else {
			yi1 <- matrix(yindexe[1:li1,],ncol=2)
			cvec <- c(0, rep(1-taui, li1) ,rep(taui,  li1))
		
			y1 <- (1-yi1[,1])/cos(alphai)
			b <- c(y1)
  
			A <- cbind(rep(1, li1), 
								diag(li1),-diag(li1))
			}
		qlp <- lp(direction = "min", 
							objective.in=cvec, const.mat=A, const.dir=rep("=",n), const.rhs=b)	
		if(transformed==TRUE)
		  return(list(q=empqf(1-qlp$solution[1]*(c(cos(alphai),sin(alphai))), yold), 
								tildeq=1-qlp$solution[1]*(c(cos(alphai),sin(alphai)))))
		else
			return(empqf(1-qlp$solution[1]*(c(cos(alphai),sin(alphai))), yold))
    }
	
  ret <- list()
  ret$bivqu <- list()
	ret$tildebivqu <- list()
	if(transformed==TRUE)
	  { for(taui in tau)
			{
			tempret <- matrix(unlist(lapply(alphaseq, FUN = get_biv, y=y, taui = taui, transformed = transformed)), ncol = 4, byrow = TRUE)
			ret$bivqu[[which(taui == tau)]] <- tempret[,1:2]
			ret$tildebivqu[[which(taui == tau)]] <- tempret[,3:4]
			}
		} else {
			for(taui in tau)
			{
			ret$bivqu[[which(taui == tau)]] <- matrix(sapply(alphaseq, FUN = get_biv, y=y, taui = taui, transformed = transformed), ncol = 2, byrow = TRUE)
			}
		ret$tildebivqu <- NULL
		}	
  names(ret$bivqu) <- sapply(tau, FUN = function(x){paste("tau = ", x, sep = "")})
  ret$y <- yold
  ret$alpha <- alphaseq
  ret$tau <- tau
	ret$tildey <- y
	ret$transformed <- transformed
  class(ret) <- c("bivquant")
  return(ret)
  }

	

	
