# gmod_tools.R
# Author: Derek Powell
# Created: 4/24/18
# Updated: 5/10/18, 4:04 PM
#
# Tools for working with graphical models. Gathering functions created across notebooks/scripts into a single set of tools
# for the project. May split them back out again if this gets too unweildy
#
# 1. Translating graphs/networks among existing R packages
# 2. Doing MCMC inference on DAGs
# 3. Creating pseudo-dbn networks

# todo: 
# 1. document functions
# 2. refactor and possibly rename funcs/arguments
# 3. explicitly map package functions with :: (done?)
# 4. reduce dependencies(?)

library(tidyverse)
library(stringr)

### ----------------------------------------
# Translation functions
# bnlearn, bdgraph, lavaan, and HydeNet


# basic string functions
cat_string_list <- function(listOfStrings) {
  paste(listOfStrings, collapse = "")
}


cat_plus <- function(x) {
  return(paste(x, "+ "))
}


cat_weight <- function(x, wVal) {
  wVal <- round(wVal, 4) # to keep human-readable
  paste0(wVal, "*", x)
}


cat_weight_list <- function(x, wVals) {
  mapply(cat_weight, x, wVals)
}


write_sum <- function(toBeSummed) {
  firstPart <- sapply(toBeSummed[1:length(toBeSummed) - 1], cat_plus)
  cat_string_list(c(firstPart, toBeSummed[length(toBeSummed)]))
}

cat_and_sort <- function(x, y) {
  sort(c(x, y)) %>% as.character(collapse = "") %>% paste(collapse = "")
}

wrap_brackets <- function(string) {
  return(paste0("[", string, "]"))
}

wrap_space <- function(string) {
  return(paste0(" ", string, " "))
}

get_left_of <- function(pattern, string) {
  trimws(strsplit(string, pattern)[[1]][1])
}

# bnlearn --> arcs_df

bnlearn_to_df <- function(bn.net) {
  arcs <- data.frame(from = bnlearn::arcs(bn.net)[,1],
                     to = bnlearn::arcs(bn.net)[,2])
  return(arcs)
}

# arcs_df --> bnlearn

get_root_nodes <- function(arcs_df) {
  # needs dplyr
  require(dplyr)
  res <- arcs_df %>%
    filter(!(from %in% to))

  unique(res$from)
}

make_bnlearn_node_string <- function(node, children) {
  output <- as.character(node)

  if (length(children) > 0) {
    output <- paste0(output, "|")
    for (c in children) {
      output <- paste0(output, c, ":")
    }
    output <- substr(output, 1, nchar(output) - 1)
  }

  return(wrap_brackets(output))
}

df_to_bnlearn_string <- function(arcs_df) {
  output <- ""
  for (node in unique(unlist(arcs_df))) {
    edges <- arcs_df[which(arcs_df$to == node), "from"]
    nodeString <- make_bnlearn_node_string(node, edges)
    output <- paste0(output, nodeString)
  }
  return(output)
}

# arcs_df --> lavaan

tag_arcs <- function(arcs_df) {
  g <- arcs_df %>%
    mutate(from = as.character(from)) %>%
    mutate(to = as.character(to)) %>%
    mutate(id = map2_chr(from, to, cat_and_sort))
  return(as_tibble(g))
}

make_lavaan_paths <- function(arcs_df) {
  # Translates bnlearn graph into lavaan syntax
  # Accepts both fully directed and partially directed graphs
  # graphStruct: $arcs object from bn object (get with get_graphStruct())

  arcs_df <- as_tibble(arcs_df)
  arcs_df <- arcs_df %>% tag_arcs()
  directed <- arcs_df[!duplicated(arcs_df$id) & !duplicated(arcs_df$id, fromLast = TRUE), ]
  modelStr <- ""

  for (node in unique(directed$to)) {
    arcs <- directed[which(directed$to == node), ]
    Ps <- write_sum(as.character(arcs$from))
    Cs <- paste(as.character(arcs$to[1]), collapse = "")
    formula <- paste(Cs, "~", Ps)

    modelStr <- paste(modelStr, "\n", formula)
  }

  undirected <- arcs_df[duplicated(arcs_df$id) | duplicated(arcs_df$id, fromLast = TRUE), ]
  undirected <- undirected[duplicated(undirected$id), ] %>% mutate(formula = paste(from, "~~", to))

  for (x in undirected$formula) {
    modelStr <- paste(modelStr, "\n", x)
  }

  return(modelStr)
}

bn_to_sem <- function(bn.net, data) {
  # Takes a bnlearn graph and returns a lavaan SEM fit to data
  # graph: bnlearn bn object
  # data: data.frame of data to fit

  graphStruct <- blearn_to_df(bn.net)

  model <- make_lavaan_paths(graphStruct)
  fit <- lavaan::sem(model = model, data = data)
  return(fit)
}

# arcs_df --> HydeNet

make_hyde_node_string <- function(node, children) {
  output <- as.character(node)

  if (length(children) > 0) {
    output <- paste0(output, "|")
    for (c in children) {
      output <- paste0(output, c, "*")
    }
    output <- substr(output, 1, nchar(output) - 1)
  }

  return(wrap_space(output))
}

df_to_hyde_string <- function(arcs_df) {
  output <- ""

  for (node in unique(unlist(arcs_df))) {
    edges <- arcs_df[which(arcs_df$to == node), "from"]
    nodeString <- make_hyde_node_string(node, edges)
    output <- paste0(output, "+", nodeString)
  }
  output <- paste0("~", substr(output, 2, nchar(output)))

  return(output)
}

bnlearn_to_hyde_string <- function(bn.net) {
  # renamed from bnlearn_to_hyde() 4/24/18, 4:04 PM
  df <- bnlearn_to_df(bn.net)
  return(df_to_hyde_string(df))
}

# arcs_df --> brms formulas

make_brms_node_string <- function(node, children) {
  output <- as.character(node)

  if (length(children) > 0) {
    output <- paste0(output, " ~ ")
    for (c in children) {
      output <- paste0(output, c, "+")
    }
    output <- substr(output, 1, nchar(output) - 1)
  }
  return(trimws(output))
}

df_to_flist <- function(arcs_df, as.strings=FALSE) {
  # take arcs_df and make list of formulas
  arcs_df$to <- as.character(arcs_df$to)
  arcs_df$from <- as.character(arcs_df$from)

  nodes <- unique(unlist(arcs_df))
  output <- lapply(nodes, function(node) {
    edges <- arcs_df[which(arcs_df$to == node), "from"]

    if (length(edges) > 0) {
      nodeString <- make_brms_node_string(node, edges)
    } else {
      nodeString <- paste0(node, " ~ 1")
    }
    node_formula <- as.formula(nodeString)

    if (as.strings) {
      return(nodeString)
    } else {
      return(node_formula)
    }
  })

  return(output)
}


flist_to_bf <- function(flist) {
  lapply(flist, brms::bf)
}


df_to_brms <- function(arcs_df) {
  flist <- lapply(df_to_flist(arcs_df), function(x) {
    brms::bf(x)
  })

  flist <- brms::mvbf(flist = flist, rescor = FALSE)
}


split_equation <- function(equation, side=c("left", "right")) {
  if (side == "left") {
    trimws(strsplit(equation, "~ ")[[1]][1])
  } else if (side == "right") {
    trimws(strsplit(equation, "~ ")[[1]][2])
  }
}


bnlearn_to_brms <- function(bn.net) {
  arcs <- as.data.frame(bnlearn::arcs(bn.net))
  return(df_to_brms(arcs))
}

# brms fit --> jags

### ----------------------------------------
### MCMC inference, prediction, and evaluation

make_jags_eq_brms <- function(to, from, coefs) {
  # to: response node
  # from: list of predictor nodes
  # coefs: broom::tidy(brms_fit) dataframe of coefficients from brms, expects brms formatted labels
  # returns jags code

  lhs <- paste0(to, "  ~ ")
  rhs <- "dnorm("
  from <- as.list(from)

  for (predictor in from) {
    coef_name <- paste0("b_", to, "_", predictor)

    coef_val <- round(coefs[which(coefs$term == coef_name), "estimate"], 5)

    rhs <- paste0(rhs, coef_val, "*", predictor, " + ")
  }
  intercept_label <- paste0("b_", to, "_Intercept")
  intercept_coef <- round(coefs[which(coefs$term == intercept_label), "estimate"], 5)

  rhs <- paste0(rhs, intercept_coef)

  s <- coefs[which(coefs$term == paste0("sigma_", to)), "estimate"]

  eq <- paste0(lhs, rhs, ", ", round(s, 5), ")")
  return(eq)
}

write_jags_model_brms <- function(arcs_df, obs_df, coefs) {
  # bn_df: arcs dataframe specifying a bayesian network
  # obs_df: observed values for variables in bn_df
  # coefs: coefficients from fitting bn in BRMS
  to <- unique(arcs_df$to)
  from <- unique(arcs_df$from)
  exogenous_nodes <- from[!(from %in% to)]

  model <- "model{\n "

  for (node in to) {
    node_parents <- arcs_df[arcs_df$to == node, "from"]
    eq <- make_jags_eq_brms(node, node_parents, coefs)
    model <- paste0(model, eq, "\n ")
  }

  for (node in exogenous_nodes) {
    node_s <- round(sd(obs_df[, node], na.rm = TRUE), 5)
    node_mean <- round(mean(obs_df[, node], na.rm = TRUE), 5)
    eq <- paste0(node, " ~ dnorm(", node_mean, ",", node_s, ")")
    model <- paste0(model, eq, "\n ")
  }
  model <- paste(model, "}")
  return(model)
}

brms_to_jags_model <- function(arcs_df, obs_df, brms_fit) {
  coefs <- broom::tidy(brms_fit)
  write_jags_model_brms(arcs_df, obs_df, coefs)
}

do_jags_inference <- function(model, nodes_to_predict, data = NULL, iter = 1e4, silent=TRUE) {
  # model: jags model string

  if (silent) {
    pbar <- "none"
  } else {
    pbar <- "text"
  }

  jags <- rjags::jags.model(textConnection(model), data = data, quiet = silent)
  model_samples <- rjags::coda.samples(jags, nodes_to_predict, iter * 2, progress.bar = "none")
  return(model_samples)
}

make_predictions <- function(model, orig_data, nodes_to_predict, nodes_held_out=c(), iter=1e4) {
  # takes a model and dataframe and generates model predictions for variables listed in
  # nodes_to_predict. By default, uses all remaining nodes in orig_data. Hold out nodes from
  # prediction with nodes_held_out. iter sets mcmc samples

  pred_df <- orig_data
  pred_df[, nodes_to_predict] <- NA

  for (i in 1:nrow(pred_df)) {
    row <- pred_df[i, ]
    row[nodes_held_out] <- NA
    pred <- do_jags_inference(model, nodes_to_predict, data = row, iter = iter)
    point_preds <- rowMeans(t(pred[[1]][(iter / 2 + 1):iter, nodes_to_predict]))
    pred_df[i, nodes_to_predict] <- point_preds
  }

  return(pred_df)
}

calc_rmse <- function(pred_df, orig_df, predicted_vars) {
  # calculate root mean squared error
  pred_df <- pred_df[, predicted_vars]
  orig_df <- orig_df[, predicted_vars]

  res <- pred_df - orig_df
  ss <- colMeans(res^2)
  rmse <- sqrt(ss)
}

calc_mae <- function(pred_df, orig_df, predicted_vars) {
  # calculate mean absolute error
  pred_df <- pred_df[, predicted_vars]
  orig_df <- orig_df[, predicted_vars]

  res <- abs(pred_df - orig_df)
  mae <- mean(res)
}

truncate_to_range <- function(x, lower, upper) {
  # bound a vector of values within a range by truncation
  if (length(dim(x)) == 2) {
    as.data.frame(apply(x, 2, function(x) {
      ifelse(x < upper, ifelse(x > lower, x, lower), upper)
    }))
  } else {
    ifelse(x < upper, ifelse(x > lower, x, lower), upper)
  }
}

calc_change_scores <- function(df) {
  # calculate change scores for data with variables and variables_t0
  var_names <- names(df)
  var_names <- var_names[!grepl("_t0", var_names)]
  for (v in var_names) {
    df[, v] <- df[, v] - df[, paste0(v, "_t0")]
  }
  df <- df %>% select(-contains("_t0"))
  return(df)
}

### ----------------------------------------
### psuedu-dbn stuff

make_dbn_df <- function(arcs_df) {
  # make a pseudo-dynamic bayesian network arcs_df
  # by adding a t0 time step node causing each variable
  arcs_df$to <- as.character(arcs_df$to)
  arcs_df$from <- as.character(arcs_df$from)
  nodes <- unique(as.character(unlist(arcs_df)))
  output <- data.frame(to = as.character(nodes))
  output$from <- paste0(output$to, "_t0")

  output$to <- as.character(output$to)
  output$from <- as.character(output$from)

  return(bind_rows(arcs_df, output))
}


