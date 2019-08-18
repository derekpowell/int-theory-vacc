hyde_to_bn_cpt <- function(hyde_cpt) {
  ndims <- length(dim(hyde_cpt))
  
  if ("xtabs" %in% class(hyde_cpt)) {
    output <- hyde_cpt/1e5
  } else {
    
    output <- aperm(hyde_cpt, ndims:1)
  }
  return (output) 
}


create_bnlearn_cog_model <- function(dag, data, evid_cpt = NULL, targetted_belief = NULL) {
  arcs_df <- bnlearn_to_df(dag)
  graph_list <- df_to_list(arcs_df)
  bagOfModels <- lapply(graph_list, function(x){fit_node_cpt(x, data)})
  model_string <- as.character(dag)
  
  if (!(is.null(evid_cpt))) {
    bagOfModels$evid <- evid_cpt
    model_string <- paste0(model_string,"[evid|", targetted_belief, "]")
  }
  model_network <- model2network(model_string)
  bnlearn_models <- lapply(bagOfModels, hyde_to_bn_cpt)
  predModel <- custom.fit(model_network, dist = bnlearn_models)
  
  return(predModel)
}


add_evid_node <- function(data, bnlearn_model, targetted_belief, evid_id = NULL) {
  
  if (!is.null(evid_id)) {
    evid <- paste0("evid",evid_id)
  } else {
    evid <- "evid"
  }
  
  evid_cpt <- make_evid_cpt(data, targetted_belief, evid)
  model_list <- lapply(bnlearn_model, function(x){x$prob})
  model_list[[evid]] <- hyde_to_bn_cpt(evid_cpt)
  
  model_string <- modelstring(bnlearn_model)
  model_string <- paste0(model_string,"[", evid,"|", targetted_belief, "]")
  
  model_network <- model2network(model_string)
  predModel <- custom.fit(model_network, dist = model_list)
  
  return(predModel)
}


add_evid_node_custom <- function(bnlearn_model, parents_string, evid_cpt) {
  
  model_list <- lapply(bnlearn_model, function(x){x$prob})
  model_list$evid <- evid_cpt
  
  model_string <- modelstring(bnlearn_model)
  model_string <- paste0(model_string,"[evid|", parents_string, "]")
  model_network <- model2network(model_string)
  predModel <- custom.fit(model_network, dist = model_list)
  
  return(predModel)
}


cogModel_predictions <- function(model){
  pred0 <- sapply(nodes(model), function(x){
    statement <- paste0("cpquery(model, event = (", x, " == 'Yes'), evidence=TRUE, n = 1e5)")
    eval(parse(text=statement))
  })
  
  pred1 <- sapply(nodes(model), function(x){
    statement <- paste0("cpquery(model, event = (", x, " == 'Yes'), evidence= (evid == 'Yes'), n = 1e5)")
    eval(parse(text=statement))
  })
  
  pred_changes <- pred1 - pred0
  
  pred_changes <- as.data.frame(pred_changes)
  
  pred_changes <- pred_changes %>%
    mutate(scale = nodes(model)) %>%
    rename(Mean = pred_changes) %>%
    mutate(type = "predicted")
  
  return(pred_changes)
}


## ----------------------------------------
## functions predictions for given evidence ratio
## as in s2-modeling.rmd notebook

predict_hypothetical <- function(df, nodes, outcome, evid_ratio){
  ## example:
  # predict_hypothetical(d_bn, nodes, "vaccIntent", 2)
  
  pred_results <- data.frame(scale=nodes, value = NA)
  
  for (i in 1:length(nodes)) {
    if (cor(df[, outcome], df[,nodes[i]]) >= 0) {
      evid_ratio <- evid_ratio
    } else {
      evid_ratio <- 1/evid_ratio
    }
    evid_cpt <- hyde_to_bn_cpt(make_evid_cpt_custom(nodes[i], evid_ratio))
    predNet <- add_evid_node_custom(base_model, nodes[i], evid_cpt)
    
    preds <- cogModel_predictions(predNet) %>%
      filter(scale==outcome)
    
    pred_results[i,"value"] <- preds[1,"Mean"]
    
  }
  
  return(pred_results)
}




filter_dscored <- function(data, scale_name){
  # prep data for estimating evid ratios (rename)
  
  data %>%
    filter(scale == scale_name) %>%
    mutate(mean = rescale_beta(mean, -3, 3)) %>%
    spread(phase, mean) %>% 
    mutate(condition = relevel(condition, ref="noInterv")) %>%
    mutate(evid = ifelse(condition=="noInterv",0,1))
}


## example:

## predModel <- add_evid_node(filter_dscored(d_scored, "diseaseSevere"), base_model, "diseaseSevere")
## pred_changes <- cogModel_predictions(predModel)