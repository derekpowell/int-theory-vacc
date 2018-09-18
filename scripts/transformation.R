boxcox_transform <- function(x, lambda2=NULL) {
	# finds optimal boxcox transformation parameters and performs transformation
  	# (requires AID and car packages)
  	# returns box-cox transformed variable
  minval <- min(x)
  
  if(minval <= 0) {lambda2 <- abs(minval) + .01}
  if(!is.null(lambda2)) {x <- x + lambda2}
  
  bc <- AID::boxcoxnc(x, lambda = seq(-10,10,.01), lambda2=lambda2, plot=FALSE, verbose=FALSE)
    
  return(car::bcPower(x, bc$lambda.hat))
}

normalize <- function(x) {
	# rescale variable to (0,1) range
    return ((x - min(x)) / (max(x) - min(x)))
  }
