###################################################################################
#function to get the indices
GetIndex_new <- function(y, alpha)
    {
	  if((!((all.equal(alpha,0)==TRUE)|(all.equal(alpha,pi/2)==TRUE))) & (alpha>0) & (alpha<pi/2))
			{
			I1 <- (1-y[,2])/sin(alpha)-(1-y[,1])/cos(alpha) >= 0
			I2 <- I1==FALSE
			
			} else if(all.equal(alpha,0)==TRUE) {
			I1 <- rep(TRUE, dim(y)[1])
	    I2 <- I1==FALSE
			} else if(all.equal(alpha,pi/2)==TRUE) {
			I1 <- rep(FALSE, dim(y)[1])
	    I2 <- I1==FALSE
			}else {
			stop("please specify alpha in [0,pi/2]")
			}
		if((sum(I1+I2)!=(dim(y)[1])))
			stop("something wrong")
		return(list(I1=I1, I2=I2))
	}
##########################################################################################
