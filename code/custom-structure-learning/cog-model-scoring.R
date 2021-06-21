## Inferring Discrete Cognitive Models from Continuous Belief Reports 
## Beta regression of belief reports

## ----------------------------------------
## functions to assemble model-fitting expressions 
## ----------------------------------------
 
wrap_paren <- function(string){
  # wrap string in parentheses
  paste0("(",string,")")
}


terms_product <- function(termlist){
  # writes formula for product of terms in termlist
  # termlist: list of strings
  # return: string
  
  output <- ""
  for (i in 1:length(termlist)) {
    output <- paste0(termlist[i],"*",output)
  }
  output <- substr(output, 1, nchar(output)-1)
  return(output)
}


assemble_equation <- function(prevterms, genterms) {
  # assembles the noisy-or-and-nor equation
  # prevterms: probability of non-prevention (p(fail prevent))
  # genterms: probability of non-generation (p(fail generate))
  # return: equation string
  
  gen <- paste0("(1-w_leak)*",terms_product(genterms))
  prev <- terms_product(prevterms)
  paste0(wrap_paren(paste0("1-",gen)), " * ", terms_product(prevterms))
}


assign_pars_exp <- function(pars){
  # creates expressions to initialize named weight parameters from par
  # input: par_weights vector
  # return: list of strings (to be evaluated as expressions)
  
  Map(function(x, i) paste0(x," <- abs(par[",i+2,"])"), pars, seq_along(pars))
}


assign_dpars_exp <- function(predictors){
  # creates expressions to initialize named direction parameters from predictors
  # input: predictors vector
  # return: list of strings (to be evaluated as expressions)
  
  Map(function(x, i) paste0("d_",x," <- ifelse(par[",i+2,"] > 0,1,0)"), predictors, seq_along(predictors))
}


create_init_expression <- function(child, predictors){
  # assembles initialization expressions
  # child: name of predicted variable
  # predictors: list of names of predictor variables (parents)
  # return: list of strings (expressions to evaluate)
  
  par_weights <- sapply(predictors, function(x){paste0("w_",x)})
  par_directions <- sapply(predictors, function(x){paste0("d_",x)})
  parameter_inits <- c(c("phi <- par[1]","w_leak <- par[2]"), assign_pars_exp(par_weights))
  
  inits <- c(parameter_inits, assign_dpars_exp(predictors), paste0("observed <- ","data$",child))
  return(inits)
}


create_equation_expression <- function(predictors){
  # creates equation string
  # predictors: list of names of predictor variables (parents)
  # return: equation string (expression to evaluate)
  
  if (length(predictors) == 0){
    equation_expression <- "w_leak"
  } else {
    par_weights <- sapply(predictors, function(x){paste0("w_",x)})
    par_directions <- sapply(predictors, function(x){paste0("d_",x)})
    
    fail_prev_terms <- sapply(predictors, function(x){
      wrap_paren(paste0("1-",x ,"*", par_weights[x], "*(1-",par_directions[x],")"))
    })
    
    fail_gen_terms <- sapply(predictors, function(x){
      wrap_paren(paste0("1-",x ,"*", par_weights[x], "*", par_directions[x]))
    })
    
    equation_expression <- assemble_equation(fail_prev_terms, fail_gen_terms)
  }
  
  return(equation_expression)
}

score_template <- function(data, par, init, equation) {
  # scoring function template, takes arbitrary init and equation strings
  eval(parse(text = init))

  predicted <- with(data, {
    eval(parse(text = equation))
  })

  alpha <- predicted * phi
  beta <- (1 - predicted) * phi

  ll <- sum(
    dbeta(observed, alpha, beta, log = TRUE)
  ) # beta regression

  if (any(par[-1] > 1 | par[-1] < -1 )){
    ll <- NaN
  }

  return(ll)
}


make_score_func <- function(child, predictors){
  # creates a partial function from score_template
  # for child and predictor variables
  init <- create_init_expression(child,predictors)
  equation <- create_equation_expression(predictors)
  func <- purrr::partial(score_template, init=init, equation=equation)
  
  return(func)
}


make_family_formula <- function(child_node, parent_nodes) {
  # function to build formulas
  # will probably need to build this into myDAGcorescore
  # or add it to BiDAG namespace, etc
  
  if (length(parent_nodes) == 0) {
    f <- as.formula(paste0(child_node,"~ 1"))
  } else {
    rhs <- ""
    for (node in parent_nodes){
      rhs <- paste0(node," + ", rhs)
    }
    f <- paste0(child_node, " ~ ", substr(rhs,1,nchar(rhs)-3))
  }
  
  return(formula(f))
}


# init_optim_pars <- function(predictors){ # old version
#   # creates pars vector to pass to optim, w/ random initializations
#   # predictors: list of predictor variables (parents)
#   # returns: named list of parameters
#   pars <- c(runif(1, min = 5, max = 20), runif(length(predictors) + 1, min = .1, max = .9))
#   # pars <- c(5, rep(.5, length(predictors) + 1))
#   names(pars) <- c(c("phi","w_leak"), sapply(predictors, function(x){paste0("w_",x)}))
#   return(pars)
# }


init_optim_pars <- function(child, predictors, data){
  # creates pars vector to pass to optim, w/ random initializations
  # child: child node being predicted
  # predictors: list of predictor variables (parents)
  # data: data for prediction
  # returns: named list of initial parameters
  
  # Code below adapted from betareg::betareg.R
  # guess for phi based on recommendation of Ferrari & Cribari-Neto, 2004 p.8
  n <- nrow(data)
  k <- ncol(data)
  f <- make_family_formula(child, predictors)
  auxreg <- lm(as.formula(f), data)
  betas <- sapply(auxreg$coefficients, function(x){min(c(max(c(x,-.9)), .9))})
  yhat <- auxreg$fitted.values
  res <- auxreg$residuals
  sigma2 <- sum(res^2)/((n - k) * (1)^2)
  phi_y <- yhat * (1 - yhat)/(n * sigma2) - 1/n
  phi <- ifelse(sum(phi_y) > 1, sum(phi_y), 1)

  pars <- c(phi, max(betas[1], .05), sapply(predictors, function(x){ betas[x] }))
  names(pars) <- c(c("phi","w_leak"), sapply(predictors, function(x){paste0("w_",x)}))

  return(pars)
}


init_optim_bounds <- function(predictors){
  # set bounds for L-BFGS-B algorithm (but that algo is not working currently - 8/16/18, 11:52 AM)
  upper <- c(1e5,rep(1,length(predictors)+1))
  lower <- c(1,rep(-1,length(predictors)+1))

  return(list(upper=upper, lower=lower))
}

bic_score_optim_result <- function(res, data){
  # calculate BIC for output of optim + data 
  # (approximate marginal likelihood P(D|Model))
  # Konishi, Sadanori; Kitagawa, Genshiro (2008). Information criteria and statistical modeling.
  # https://en.wikipedia.org/wiki/Bayesian_information_criterion
  
  k <- length(res$par) # num parameters
  n <- nrow(data) # num of obs
  ll <- res$value # log likelihood
  
  k*log(n) - 2*ll
}

## Scoring Functions Code

# score_family <- function(node, predictors, data, method = "Nelder-Mead") {
#   # return marginal likelihood of data for family
#   # P(data | graph)
#   # node: child node (char)
#   # predictors: parent nodes (char list)
#   # data: dataframe
#   pars <- init_optim_pars(node, predictors, data)
#   score <- make_score_func(node, predictors)
#   result <- suppressWarnings(optim(
#     par = pars,
#     fn = score,
#     method = method, # BFGS need to add error catching!
#     control = list(fnscale = -1),
#     data = data
#   ))
# 
#   b <- -1*bic_score_optim_result(result, data)/2
#   # b <- bic_score_optim_result(result, data)/2
# 
#   return(b)
#   # return(as.numeric(result$value))
# }

# score_family <- function(node, predictors, data, method = "Nelder-Mead") {
#   # return marginal likelihood of data for family
#   # P(data | graph)
#   # node: child node (char)
#   # predictors: parent nodes (char list)
#   # data: dataframe
#   pars <- init_optim_pars(node, predictors, data)
#   score <- make_score_func(node, predictors)
#   result <- suppressWarnings(optim(
#     par = pars,
#     fn = score,
#     method = method, # BFGS need to add error catching!
#     control = list(fnscale = -1),
#     data = data
#   ))
#   
#   b <- -1*bic_score_optim_result(result, data)/2
#   # b <- bic_score_optim_result(result, data)/2
#   
#   return(b)
#   # return(as.numeric(result$value))
# }


# # Model Fitting

find_node_coefs <- function(child, predictors, data) {
  # should rewrite to cover whole discrete space
  score <- make_score_func(child, predictors)
  result <- optim(
    par = init_optim_pars(child, predictors, data),
    fn = score,
    # upper = bounds$upper,
    # lower = bounds$lower,
    # method = "BFGS",
    control = list(fnscale = -1),
    data = data
  )
  
  return(result$par)
}







## Another set of scoring functions that scores all possible causal directions for each predictor
## should be better behaved optimization-wise
## In progress 8/23/18, 3:01 PM

## ----------------------------------------

pass_dpars_exp <- function(predictors, dirs){
  # creates expressions to initialize named direction parameters from predictors
  # input: predictors vector
  # return: list of strings (to be evaluated as expressions)
  
  # Map(function(x, i) paste0("d_",x," <- ifelse(par[",i+2,"] > 0,1,0)"), predictors, seq_along(predictors))
  map2(predictors, dirs, function(p,d){paste0("d_",p," <- ",d)})
}

assign_pars_exp2 <- function(pars){
  # creates expressions to initialize named weight parameters from par
  # input: par_weights vector
  # return: list of strings (to be evaluated as expressions)
  
  Map(function(x, i) paste0(x," <- par[",i+2,"]"), pars, seq_along(pars))
}

create_init_expression2 <- function(child, predictors, dirs){
  # assembles initialization expressions
  # child: name of predicted variable
  # predictors: list of names of predictor variables (parents)
  # dirs: list/vector of directions
  # return: list of strings (expressions to evaluate)

  par_weights <- sapply(predictors, function(x){paste0("w_",x)})
  # par_directions <- sapply(predictors, function(x){paste0("d_",x)})
  parameter_inits <- c(c("phi <- par[1]","w_leak <- par[2]"), assign_pars_exp2(par_weights))

  # inits <- c(parameter_inits, paste0("observed <- ","data$",child))
  inits <- c(parameter_inits, pass_dpars_exp(predictors, dirs), paste0("observed <- ","data$",child))
  return(inits)
}


init_optim_pars2 <- function(child, predictors, data){
  # creates pars vector to pass to optim, w/ random initializations
  # child: child node being predicted
  # predictors: list of predictor variables (parents)
  # data: data for prediction
  # returns: named list of initial parameters
  
  # Code below adapted from betareg::betareg.R
  # guess for phi based on recommendation of Ferrari & Cribari-Neto, 2004 p.8
  n <- nrow(data)
  k <- ncol(data)
  f <- make_family_formula(child, predictors)
  auxreg <- lm(as.formula(f), data)
  betas <- sapply(auxreg$coefficients, function(x){min(c(max(c(x,.1)), .90))})
  yhat <- auxreg$fitted.values
  res <- auxreg$residuals
  sigma2 <- sum(res^2)/((n - k) * (1)^2)
  phi_y <- yhat * (1 - yhat)/(n * sigma2) - 1/n
  phi <- ifelse(sum(phi_y) > 1, sum(phi_y), 1)
  
  pars <- c(phi, .5, sapply(predictors, function(x){ betas[x] }))
  names(pars) <- c(c("phi","w_leak"), sapply(predictors, function(x){paste0("w_",x)}))
  
  return(pars)
}


make_score_func2 <- function(child, predictors, dirs){
  # creates a partial function from score_template
  # for child and predictor variables
  init <- create_init_expression2(child,predictors, dirs)
  equation <- create_equation_expression(predictors)
  func <- purrr::partial(score_template2, init=init, equation=equation)

  return(func)
}

score_template2 <- function(data, par, init, equation) {
  # scoring function template, takes arbitrary init and equation strings
  eval(parse(text = init))
  
  predicted <- with(data, {
    eval(parse(text = equation))
  })
  
  alpha <- predicted * phi
  beta <- (1 - predicted) * phi
  
  ll <- sum(
    dbeta(observed, alpha, beta, log = TRUE)
  ) # beta regression
  
  # if (any(par[-1] > 1 | par[-1] < 0 )){
  #   ll <- NaN
  # }
  
  if(!is.finite(ll)){ ll <- -1e300}
  return(ll)
}


score_family <- function(node, predictors, data) { # comment to switch back to original!
  # return marginal likelihood of data for family
  # P(data | graph)
  #
  # deals with discontinuity by computing over space of discrete parameters
  # that means it's 2^n for number of parents n
  # But avoids discontinuities so can use l-bfgs-b for ~2x speed-up 
  #
  # node: child node (char)
  # predictors: parent nodes (char list)
  # data: dataframe
  pars <- init_optim_pars2(node, predictors, data) # need to update this so no negs

  if (length(predictors > 0 )){
    dirlist <- rep(list(c(0,1)), length(predictors))
    names(dirlist) <- map_chr(predictors, function(x){paste0("d_",x)})
    dirs <- expand.grid(dirlist)
    
    result_list <- list()
    

    # registerDoFuture() 
    
    # result_list <- future_lapply( 1:nrow(dirs), function(i) {
    # foreach(i = 1:nrow(dirs)) %dopar% {
    for (i in 1:nrow(dirs)){
      score <- make_score_func2(node, predictors, dirs[i,])
      
      result <- suppressWarnings(optim(
        par = pars,
        fn = score,
        method = "L-BFGS-B", # BFGS need to add error catching!
        upper = c(Inf, 1, rep(1,length(predictors))), # maybe constrain to .99?
        lower = c(0, 0, rep(0,length(predictors))),
        control = list(fnscale = -1),
        data = data
      ))
      
      result_list[[i]] <- result
    }
    
    ll <- max( map_dbl(result_list, function(x){x$value}) ) # take maximum likelihood score
  
  } else {
    score <- make_score_func2(node, predictors,c())
    
    result <- suppressWarnings(optim(
      par = pars,
      fn = score,
      method = "L-BFGS-B", # BFGS need to add error catching!
      upper = c(Inf, 1, rep(1,length(predictors))),
      lower = c(0, 0, rep(0,length(predictors))),
      control = list(fnscale = -1),
      data = data
    ))
    
    ll <- result$value
  }
  
  
  k <- length(predictors)*2 + 2 # num parameters
  n <- nrow(data) # num of obs
  
  bic_val <- k*log(n) - 2*ll
  
  marglikli <- -1 * bic_val/2
  return(marglikli)
}


find_node_coefs <- function(node, predictors, data) {
  # searches full discrete space
  pars <- init_optim_pars2(node, predictors, data) # need to update this so no negs
  
  if (length(predictors > 0 )){
    dirlist <- rep(list(c(0,1)), length(predictors))
    names(dirlist) <- map_chr(predictors, function(x){paste0("d_",x)})
    dirs <- expand.grid(dirlist)
    
    result_list <- list()
    
    
    # registerDoFuture() 
    
    # result_list <- future_lapply( 1:nrow(dirs), function(i) {
    # foreach(i = 1:nrow(dirs)) %dopar% {
    for (i in 1:nrow(dirs)){
      score <- make_score_func2(node, predictors, dirs[i,])
      
      result <- suppressWarnings(optim(
        par = pars,
        fn = score,
        method = "L-BFGS-B", # BFGS need to add error catching!
        upper = c(Inf, 1, rep(1,length(predictors))),
        lower = c(0, 0, rep(0,length(predictors))),
        control = list(fnscale = -1),
        data = data
      ))
      
      result_list[[i]] <- result
    }
    
    ll_list <- map_dbl(result_list, function(x){x$value})
    ll <- max(ll_list)
    
    ll_index <- match(ll, ll_list)
    
    result <- result_list[[ll_index]]
    
    coefs <- result$par
    dirs <- unlist(dirs[ll_index,])
    
    coefs[-1:-2] <- coefs[-1:-2] * 2*(dirs - .5)
    
    # ll <- max( map_dbl(result_list, function(x){x$value}) ) # take maximum likelihood score
    
  } else {
    score <- make_score_func2(node, predictors,c())
    
    result <- suppressWarnings(optim(
      par = pars,
      fn = score,
      method = "L-BFGS-B", # BFGS need to add error catching!
      upper = c(Inf, 1, rep(1,length(predictors))),
      lower = c(0, 0, rep(0,length(predictors))),
      control = list(fnscale = -1),
      data = data
    ))
    
    coefs <- result$par
  }
  
  return(coefs)
}


coefs_to_dirs <- function(coef_list){
  parent_coefs <- coef_list[-(1:2)]
  dirs <- .5+sign(parent_coefs)*.5 # recode -1,1 to 0,1
  names(dirs) <- gsub("w_","d_",names(dirs)) # rename w_ to d_
  
  return(dirs)
}


cog_model_score_func <- function(node, parents, data, args=NULL){
  
  if (length(args)==0){ # gets passed as empty list I guess? 
    return(score_family(node, parents, data))
    
  } else {
    coefs <- find_node_coefs(node, parents, args$train_data)
    dirs <- coefs_to_dirs(coefs)
    
    score_func <- make_score_func2(node, parents, dirs)
    
    ll <- score_func(data, abs(coefs))
    
    # # can see this works correctly and returns same value if data and args$train_data are the same
    # k <- length(parents)*2 + 2 # num parameters
    # n <- nrow(data) # num of obs
    # 
    # bic_val <- k*log(n) - 2*ll
    # 
    # marglikli <- -1 * bic_val/2
    
    return(ll) # but makes more sense to return log-likelihood for new predictions
  }
  
}

## tests -- this works!
# cog_model_score_func("C", c("A","B"), dsim[1:400,])
# cog_model_score_func("C", c("A","B"), dsim[401:500,], args=list(train_data=dsim[1:400,]))
# cog_model_score_func("C", c("A","B"), dsim[1:400,], args=list(train_data=dsim[1:400,]))
# cog_model_score_func("C", c(), dsim[1:400,])
# cog_model_score_func("C", c(), dsim[1:400,], args=list(train_data=dsim[1:400,]))
