## Inferring Discrete Cognitive Models from Continuous Belief Reports 
## Beta regression of belief reports

## initialize ...

source("cog-model-scoring.R")
source("cog-model-jags-tools.R")
source("cog-model-bidag-functions.R")
source("cog-model-discrete-funcs.R")

rescale_beta <- function(x, lower, upper) {
  # rescales onto the open interval (0,1)
  # rescales over theoretical bounds of measurement, specified by "upper" and "lower"
  
  N <- length(x)
  N <- min(20, N) # bound to [.025,.975] 
  res <- (x - lower) / (upper - lower)
  res <- (res * (N - 1) + .5) / N
  
  return(as.vector(res))
}
