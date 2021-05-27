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
  # this is a bit hacky and maybe could be made more elegant but probably doesn't matter
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

to_odds <- function(p){
  p/(1-p)
}

to_prob <- function(x){
  x/(1+x)
}


# get_evid_probs <- function(result){
#   # this is a heuristic/hack, could be made better but probably doesn't matter 
#   # 5/25/21, 3:52 PM thought this would be better but it's bad actually
#   log_evid_ratio <- result$par[1] + result$par[2]
#   p0 <- .5
#   # logit_p1 <- plogis(p0) + log_evid_ratio
#   p1 <- to_prob( exp(logit(p0) + log_evid_ratio))
# 
#   c(p0, p1)
# }


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


cogModel_predictions <- function(model){
  pred_nodes <- nodes(model)
  pred_nodes <- pred_nodes[pred_nodes!="evid"]
  
  pred0 <- sapply(pred_nodes, function(x){
    statement <- paste0("cpquery(model, event = (", x, " == 'Yes'), evidence = TRUE, n = 5e5)")
    # statement <- paste0("cpquery(model, event = (", x, " == 'Yes'), evidence = (evid=='No'), n = 5e5)")
    eval(parse(text=statement))
  })
  
  pred1 <- sapply(pred_nodes, function(x){
    statement <- paste0("cpquery(model, event = (", x, " == 'Yes'), evidence = (evid == 'Yes'), n = 5e5)")
    eval(parse(text=statement))
  })
  
  pred_changes <- pred1 - pred0
  
  pred_changes <- as.data.frame(pred_changes)
  
  pred_changes <- pred_changes %>%
    mutate(scale = pred_nodes) %>%
    rename(Mean = pred_changes) %>%
    mutate(type = "predicted")
  
  return(pred_changes)
}


make_pred_df <- function(data, base_model, target_node){
  predModel <- add_evid_node(filter(data, scale == target_node), base_model, target_node)
  pred_changes <- cogModel_predictions(predModel)
  
  # generate model predictions + create plot
  obs_changes <- data %>%
    mutate(changeScore = post-pre) %>%
    filter(evid==1) %>%
    group_by(scale) %>%
    summarise(
      Mean = mean(changeScore),
      ul = mean(changeScore) + 1.96*sd(changeScore)/sqrt(n()), # replace with bootstrapping
      ll = mean(changeScore) - 1.96*sd(changeScore)/sqrt(n())
    ) %>%
    mutate(type = "observed")
  
  obs_pred <- bind_rows(pred_changes, obs_changes) %>%
    spread(type, Mean) %>%
    group_by(scale) %>%
    summarize_all(mean, na.rm=TRUE)
  
  return(obs_pred)
}
