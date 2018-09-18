# # Author: Derek Powell
# # Updated: 5/14/18, 2:57 PM

# # code for applying virtual evidence to cognitive model

# functions for creating virtual evidence formula links

make_ve_formula_strings <- function(nodes, suffix = "_v") {
  lapply(nodes, function(x){paste0(x, suffix, " ~ ", x)})
}

# functions for fitting the weights of those links
fit_formulas <- function(formula_strings, data, family="gaussian") {
  # take a list of formulas as strings and return fit glm objects
  # formula_strings: strings specifying formulas
  # data: data for analyses
  # family: glm family string or function object
  
  if (family %in% c("gaussian","binomial","gamma")) {# etc, finish later
    lapply(formula_strings, function(x) {
      glm(as.formula(x), data = data, family = family)
    })
  } else if (family=="beta") {
    require(betareg)
    lapply(formula_strings, function(x) {
      betareg::betareg(as.formula(x), data = data)
    })
  }
}

# A function to assemble the results with the original model
# 1. get lm coefficients and make jags code (modify make_jags_eq for lm())

glm_to_jags_eq <- function(fit, digits=5) {
  # fit: glm fit object
  # returns jags code
  
  if (class(fit) == "glm") {
    response_node <- as.character(formula(fit))[2]
    coefs <- round(fit$coefficients, digits)
    linear_eq <- paste0(coefs[1])
    
    if (length(coefs) > 1) {
      
      for (i in 2:length(coefs)) {
        
        coef_name <- names(coefs)[i]
        coef <- coefs[i]
        linear_eq <- paste0(linear_eq, " + ", coef, "*", coef_name)
      }
    }
    if (fit$family$family == "gaussian") {
      precision <- 1/stats::sigma(fit)^2
      paste0(response_node," ~ ", "dnorm(", linear_eq,", ", round(precision,5),")")
    } else if (fit$family$family == "binomial") {
      paste0(response_node, " ~ ", "dbern(ilogit(", linear_eq, "))")
    }
  } else if (class(fit) == "betareg") {
    
    response_node <- as.character(formula(fit))[2]
    coefs <- round(fit$coefficients$mean, digits)
    linear_eq <- paste0(coefs[1])
    
    if (length(coefs) > 1) {
      
      for (i in 2:length(coefs)) {
        
        coef_name <- names(coefs)[i]
        coef <- coefs[i]
        linear_eq <- paste0(linear_eq, " + ", coef, "*", coef_name)
      }
    }
    
    phi <- fit$coefficients$precision
    # do something ...
    paste0(response_node, 
           " ~ ", 
           "dbeta(ilogit(", 
           linear_eq, 
           ")*", 
           phi, 
           ",", 
           "(1-ilogit(", 
           linear_eq, 
           "))*", 
           phi,
           ")",
           " T(0.001,0.999)") # https://stackoverflow.com/questions/47135726/error-slicer-stuck-at-value-with-infinite-density-running-binomial-beta-model
  }
  
  # add support for other families later ...
}


# 2. assemble whole jags model

fit_ve_nodes <- function(ve_df, suffix = "_v", family="gaussian") {
  # estimate coefficients for virtual evidence links
  # ve_df: dataframe with test-retest for virtual evidence data
  # suffix: suffix indicating ve variables
  # family: glm family string
  
  main_nodes <- variable.names(ve_df)[!grepl(suffix, variable.names(ve_df))]
  model_strings <- make_ve_formula_strings(main_nodes, suffix = suffix)
  ve_fits <- fit_formulas(model_strings, ve_df, family)
}


write_jags_model <- function(arcs_df, obs_df, family = "gaussian", ve_df = NULL, ve_suffix = NULL) {
  model_strings <- df_to_flist(arcs_df)
  fits <- fit_formulas(model_strings, obs_df, family)

  if (!is.null(ve_df)) {
    ve_fits <- fit_ve_nodes(ve_df, ve_suffix, family)
    fits <- c(fits, ve_fits)
  }
  
  jags_eqs <- lapply(fits, glm_to_jags_eq)
  model <- "model{\n "
  
  for (eq in jags_eqs) {
    model <- paste0(model, eq, "\n ")
  }
  
  model <- paste0(model,"}")
}


add_ve_to_model <- function(jags_model, ve_df, ve_suffix = "_v", family = "gaussian") {
  ve_fits <- fit_ve_nodes(ve_df, ve_suffix, family)
  jags_eqs <- lapply(ve_fits, glm_to_jags_eq)
  
  model <- substr(jags_model,1,nchar(jags_model)-1)
  for (eq in jags_eqs) {
    model <- paste0(model, eq, "\n ")
  }
  
  model <- paste0(model,"}")
}
   
# # Generate Predictions and score model performance ----------------------------------------
# # later: additional functions to score fit and calculate likelihood, AIC, MAE, and RMSE 
# # For now, will stick with MAE amnd RMSE, though really want AIC

do_jags_inference <- function(model, variable.names, data = NULL, iter = 1e4, silent=TRUE) {
  if (silent) {
    pbar <- "none"
  } else {
    pbar <- "text"
  }
  
  jags <- jags.model(textConnection(model), data = data, quiet=silent)
  model_samples <- coda.samples(jags, variable.names, iter*2, progress.bar="none")
  return(model_samples)
}

make_predictions <- function(model, orig_data, nodes_to_predict, nodes_held_out=c(), iter=1e4) {
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

calc_MAE <- function(true_data, predictions) {
  resids <- true_data - predictions
  MAE <- colMeans(abs(resids))
  return(MAE)
}

calc_RMSE <- function(true_data, predictions) {
  resids <- true_data - predictions
  rmse <- sqrt(colMeans(resids^2))
  return(rmse)
}

# # function to intervention node to cognitive model (or any other node) ...

add_node <- function(jags_model, node_eq) {
  # append an equation to the end of a jags model
  # jags_model: jags model string
  # node_eq: equation to add
  model <- substr(jags_model,1,nchar(jags_model)-1)
  model <- paste0(model, node_eq, "\n}")
  return(model)
}

# # tests -----------------------------------------
# n_rows <- 100
# arc_df <- data.frame(from = c("A","A","B"), to = c("B","C","C"))
# sim_df <- data.frame(A = rnorm(n_rows, 0, 1)) %>%
#   mutate(B = .5 * A + rnorm(n_rows, 0, 1)) %>%
#   mutate(C = .5 * A + .25 * B + rnorm(n_rows, 0, 1)) %>%
#   mutate(
#     A_v = A + rnorm(n_rows,0,.5),
#     B_v = B + rnorm(n_rows,0,.5),
#     C_v = C + rnorm(n_rows,0,.5)
#   )
# 
# jags_model <- write_jags_model(arc_df, sim_df, ve_df = sim_df, ve_suffix = "_v")
# 
# preds <- make_predictions(jags_model, sim_df, c("A","B","C"))
# 
# calc_MAE(sim_df, preds)
# calc_RMSE(sim_df, preds)
