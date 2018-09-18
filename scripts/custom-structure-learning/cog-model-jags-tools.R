# # Doing inference on beta DAGs through rjags
# # ----

do_jags_inference <- function(model, variable.names, data = NULL, iter = 1e4, silent=TRUE) {
  if (silent) {
    pbar <- "none"
  } else {
    pbar <- "text"
  }
  tc <- textConnection(model)
  jags <- rjags::jags.model(tc, data = data, quiet=silent)
  close(tc)
  
  model_samples <- rjags::coda.samples(jags, variable.names, iter*2, progress.bar="none")
  return(model_samples)
}


make_predictions <- function(model, orig_data, nodes_to_predict, nodes_held_out=c(), iter=1e3) {
  # takes a model and dataframe and generates model predictions for variables listed in 
  # nodes_to_predict. By default, uses all remaining nodes in orig_data. Hold out nodes from
  # prediction with nodes_held_out. iter sets mcmc samples
  
  pred_df <- orig_data
  pred_df[,nodes_to_predict] <- NA
  
  for (i in 1:nrow(pred_df)) {
    row <- pred_df[i,]
    row[nodes_held_out] <- NA
    pred <- do_jags_inference(model, nodes_to_predict, data = row, iter = iter)
    point_preds <- rowMeans(t(pred[[1]][(iter/2+1):iter,nodes_to_predict]))
    pred_df[i,nodes_to_predict] <- point_preds
  }
  
  return(pred_df)
}

assemble_fitted_equation <- function(w_leak, prevterms, genterms) {
  # assembles the noisy-or-and-nor equation
  # prevterms: probability of non-prevention (p(fail prevent))
  # genterms: probability of non-generation (p(fail generate))
  # return: equation string
  
  gen <- paste0("(1-",w_leak, ")*",terms_product(genterms))
  prev <- terms_product(prevterms)
  paste0(wrap_paren(paste0("1-",gen)), " * ", terms_product(prevterms))
}


create_fitted_equation_expression <- function(predictors, coefs, digits = 5){
  # creates equation string
  # predictors: list of names of predictor variables (parents)
  # return: equation string (expression to evaluate)
  
  coefs <- as.list(coefs)
  
  w_leak <- as.character(abs(round(coefs[["w_leak"]], digits)))
  
  if (length(predictors) == 0){
    equation_expression <- w_leak
  } else {
    par_weights <- sapply(predictors, function(x){as.character(abs(round(coefs[[paste0("w_",x)]],digits)))})
    par_directions <- sapply(predictors, function(x){ifelse(coefs[[paste0("w_",x)]] > 0, 1, 0)})
    
    fail_prev_terms <- sapply(predictors, function(x){
      wrap_paren(paste0("1-",x ,"*", par_weights[x], "*(1-",par_directions[x],")"))
    })
    
    fail_gen_terms <- sapply(predictors, function(x){
      wrap_paren(paste0("1-",x ,"*", par_weights[x], "*", par_directions[x]))
    })
    
    equation_expression <- assemble_fitted_equation(w_leak, fail_prev_terms, fail_gen_terms)
  }
  
  return(equation_expression)
}

make_jags_mu_eq <- function(to, from, coefs) {
  rhs <- paste0("mu.",to, " = ")
  lhs <- create_fitted_equation_expression(from, coefs)
  
  paste0(rhs,lhs)
}

break_lines <- function(strings){
  output <- ""
  for (x in strings) {
    output <- paste0(output, x,"\n")
  }
  
  substr(output,1,nchar(output)-1)
} 

make_jags_eq <- function(to, from, coefs) {
  # to: response node
  # from: list of predictor nodes
  # coefs: named list from optim results
  # returns jags code
  
  # lhs <- paste0(to,"  ~ ")
  # rhs <- "dbeta("
  # from <- as.list(from)
  mu_str <- make_jags_mu_eq(to, from, coefs)
  phi_str <- paste0("phi.", to, " = ", round(coefs[1], 5))
  
  alpha <- paste0("alpha.",to," = mu.",to, "*","phi.",to)
  beta <- paste0("beta.",to," = (1-(mu.",to ,"))*","phi.",to)
  eq <- paste0(to, " ~ ", "dbeta(","alpha.",to, ", beta.",to,")", " T(0.000001,0.999999)")
  
  all_pars <- break_lines(c(mu_str,phi_str,alpha,beta,eq))
  return(all_pars)
}




write_jags_model <- function(bn.dag, obs_df) {
  # bn.dag: bnlearn dag object
  # obs_df: observed values for variables in bn_df
  # coefs: coefficients from fitting bn in BRMS
  
  ordered_nodes <- bnlearn::node.ordering(bn.dag)
  
  # to <- unique(as.character(bn_df$to))
  # from <- unique(as.character(bn_df$from))
  # exogenous_nodes <- from[!(from %in% to)]
  
  model <- "model{\n"
  
  for (node in ordered_nodes) {
    # node <<- as.character(node)
    # node_parents <- as.character(bn_df[bn_df$to==node,"from"])
    node_parents <- parents(bn.dag, node)
    coefs <- find_node_coefs(node, node_parents, obs_df)
    eq <- make_jags_eq(node, node_parents, coefs)
    model <- paste0(model,eq,"\n")
  } 
  
  model <- paste(model,"}")
  return(model)
}