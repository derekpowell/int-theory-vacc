## functions for estimating evidence ratio and augmenting models

inv_logit <- plogis
logit <- qlogis

evid_score_func <- function(data, par) {
  
  inv_logit <- plogis
  logit <- qlogis
  
  pred_y <- with(data,
                 {
                   inv_logit(logit(pre) + evid*par[2] + par[1])
                 }
  )
  pred_beta_A <- pred_y * par[3]
  pred_beta_B <- (1-pred_y) * par[3]
  
  ll <- -1 * sum(dbeta(data$post, pred_beta_A, pred_beta_B, log=TRUE)) # beta regression
  return(ll)
}

get_evid_probs <- function(result){
  # this is a heuristic/hack, could be made better but probably doesn't matter
  evid_ratio <- exp(result$par[1] + result$par[2])
  if (evid_ratio < 1) {
    p0 <- .90
  } else {
    p0 <- .50
  }
  
  if (p0*evid_ratio >= 1) {
    p0 <- p0*.9/(p0*evid_ratio)
  }
  
  p1 <- p0*evid_ratio
  
  c(p0,p1)
}


evid_probs_to_cpt <- function(probs, var_name, evid_name = "evid"){
  p0 <- probs[1]
  p1 <- probs[2]
  
  full_list <- paste0("list(",var_name,"=c('Yes','No'),", evid_name, " = c('Yes', 'No'))")
  make_custom_cpt(c(p1, p0, (1-p1), (1-p0)), c(2,2), eval(parse(text=full_list)), NULL)
}


find_evid_ratio <- function(data){
  optim(par = c(0, 0, 5), fn = evid_score_func, data = data, method="BFGS")
}


make_evid_cpt <- function(data, var_name, evid_name = "evid"){
  result <- find_evid_ratio(data)
  probs <- get_evid_probs(result)
  cpt <- evid_probs_to_cpt(probs, var_name, evid_name)
  
  return(cpt)
}

make_evid_cpt_custom <- function(var_name, evid_ratio, evid_name = "evid"){
  # result <- find_evid_ratio(data)
  result <<- list(par = c(0, log(evid_ratio)) )
  probs <- get_evid_probs(result)
  cpt <- evid_probs_to_cpt(probs, var_name, evid_name)
  
  return(cpt)
}