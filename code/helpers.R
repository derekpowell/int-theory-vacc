## Helpers for int-theory-vacc project

devtools::source_gist(id = "f1994c0f8325abbc5d300600744af39d", filename="cbrm.R")

rescale_beta <- function(x, lower, upper) {
  # rescales onto the open interval (0,1)
  # rescales over theoretical bounds of measurement, specified by "upper" and "lower"
  
  N <- length(x)
  res <- (x - lower) / (upper - lower)
  res <- (res * (N - 1) + .5) / N
  
  return(as.vector(res))
}


rcor <- function(x,y=NULL){
  r <- cor(x,y)
  return( round(r,3) )
}


