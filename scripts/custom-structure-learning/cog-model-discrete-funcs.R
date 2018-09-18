make_custom_cpt <- function(ps, dims, dim_names, data) {
  # makes custom cpt in HydeNet format
  output <- array(ps, dim=dims, dimnames=dim_names)
  attr(output, "model") <- data
  attr(output, "class") <- c("cpt","array")
  
  return(output)
}


result_to_cpt <- function(result, child, predictors, data) {
  # compute cpt for child node given optim result and predictor variables
  # result: optim result
  # child: name of child node (string)
  # data: dataframe
  # return: array
  states <- lapply(c(predictors,child), function(x){c(0,1)})
  names(states) <- c(predictors,child)
  states_labeled <- lapply(c(predictors,child), function(x){c("No","Yes")}) # not working 8/17/18, 9:44 AM
  names(states_labeled) <- c(predictors,child)
  ds <- ifelse(result > 0,1,0)
  names(ds) <- sapply(names(result), function(x){paste0("d_",substr(x,3,100))})
  result <- abs(c(result, ds[-(1:2)]))
  
  z <- expand.grid(states) %>%
    mutate(
      p = with(
        as.data.frame(t(result)),
        eval(parse(text=create_equation_expression(predictors)))
      )
    )
  z[which(z[child]==0),"p"] <- 1-z[which(z[child]==0),"p"] # unnecessary to preserve child name oh well
  
  output <- make_custom_cpt(z$p, rep(2,length(predictors)+1), states_labeled, data)
  return(output)
}


fake_xtabs <- function(node, data, trueval = 1, falseval = 0){
  data <- as.data.frame(data)
  probs <- data[, node]
  p1 <- mean(probs)
  count1 <- round(p1*1e5)
  count0 <- 1e5 - count1
  counts <- matrix(c(rep(trueval,count1), rep(falseval,count0)), ncol=1)
  colnames(counts) <- c(node)
  f <- paste0("~",node)
  output <- xtabs(as.formula(f), data = counts)
  attr(output, "call") <- list(formula = as.formula(f) )
  
  return(output)
}


fit_node_cpt <- function(node, data) {
  child <- node$child # specific to format above ,could be imrpoved
  predictors <- node$parents # specific to format above ,could be imrpoved
  data <- as.data.frame(data)
  if (is.null(predictors)) {
    node_cpt <- fake_xtabs(child, data, "Yes","No")
  } else {
    
    result <- find_node_coefs(child, predictors, data)
    
    node_cpt <- result_to_cpt(result, child, predictors, data)
  }
  
  return(node_cpt)
}


# # Tests and examples
# # ----------------------------------------

# # write pipeline to go from data input to estimated "mental model"
# # 1. need to parse dag into separate models for each node

# # bnlearn::parents(x, node) for each to get list of nodes + parents
# # for now, manually create a list

graph_list <- list( # simpler list for df_cog
  A = list(child = "A", parents = c()),
  B = list(child = "B", parents = c()),
  C = list(child = "C", parents = c("A","B"))
)

# # 2. run optim for each family & 3. spit out cpt object for each family w/ function fit_node()
# bagOfModels <- lapply(graph_list, function(x){fit_node_cpt(x, df_cog)})
# 
# # 4. assemble into hydenet model
# bagNet <- HydeNetwork(bagOfModels)
# writeNetworkModel(bagNet, pretty=TRUE)
