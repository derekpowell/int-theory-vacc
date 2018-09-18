# author: Derek Powell
# created: 11/3/17, 10:28 AM
# last updated: 11/3/17, 1:44 PM

# script for translating bnlearn arcs to lavaan

library(tidyverse)
library(stringr)

cat_string_list <- function(listOfStrings) {
  paste(listOfStrings, collapse="")
}


cat_plus <- function(x) {
  return(paste(x,"+ "))
}


cat_weight <- function(x,wVal) {
  wVal <- round(wVal,4) # to keep human-readable
  paste0(wVal,"*",x)
}


cat_weight_list <- function(x,wVals) {
  mapply(cat_weight,x,wVals)
}


write_sum <- function(toBeSummed) {
  firstPart <- sapply(toBeSummed[1:length(toBeSummed)-1],cat_plus)
  cat_string_list(c(firstPart,toBeSummed[length(toBeSummed)]))
}


cat_and_sort <- function(x,y) {
  sort(c(x,y)) %>% as.character(collapse="") %>% paste(collapse="")
}


tag_arcs <- function(graphStruct) {
  g <- graphStruct %>% 
    mutate(from =  as.character(from)) %>%
    mutate(to = as.character(to)) %>%
    mutate(id = map2_chr(from, to, cat_and_sort))
  return(as_tibble(g))
}

get_graphStruct <- function(bnObject) {
  # Returns the graph structure ($arcs) from a bn object
  # bnObject: any object of class bn
  
  graphStruct <- bnObject$arcs %>% as_tibble()
  return(graphStruct)
}

make_weighted_lavaan_paths <- function(graphStruct, pathParams, residParams) {
  # Translates bnlearn graph into lavaan syntax with randomly sampled weights
  # for simulating data from graphical models
  # Note: this function requires fully directed graph
  # ---
  # baseGraph: DAG from bnlearn
  # pathParams: c(a,b) for Beta(a, b) for path coefficients 
  # residParams: c(a,b) for Beta(a, b) for variances
  
  graphStruct <- graphStruct %>% mutate(from = as.character(from),
                                        to = as.character(to))
  allNodes <- unique(c(graphStruct$from, graphStruct$to)) %>% as.character()
  
  modelStr <- ""
  for (node in unique(graphStruct$to)) {
    arcs <- graphStruct[which(graphStruct$to==node),]
    weights <- rbeta(n=dim(arcs)[1], shape1=pathParams[1], shape2=pathParams[2])
    Ps <- write_sum(cat_weight_list(arcs$from, weights))
    Cs <- paste(as.character(arcs$to[1]),collapse="")
    formula <- paste(Cs,"~",Ps)
    
    modelStr <- paste(modelStr,"\n",formula)
  }
  
  for (node in allNodes) {
    resWeight <- rbeta(n=1, shape1=residParams[1], shape2=residParams[2])
    formula <- paste(node, " ~~ ", resWeight)
    modelStr <- paste(modelStr,"\n",formula)
  }
  
  return(modelStr)
}


make_lavaan_paths <- function(graphStruct) {
  # Translates bnlearn graph into lavaan syntax
  # Accepts both fully directed and partially directed graphs
  # graphStruct: $arcs object from bn object (get with get_graphStruct())

  graphStruct <- as_tibble(graphStruct)
  graphStruct <- graphStruct %>% tag_arcs()
  directed <- graphStruct[!duplicated(graphStruct$id) & !duplicated(graphStruct$id, fromLast=TRUE),]
  modelStr <- ""

  for (node in unique(directed$to)) {
    arcs <- directed[which(directed$to==node),]
    Ps <- write_sum(as.character(arcs$from))
    Cs <- paste(as.character(arcs$to[1]),collapse="")
    formula <- paste(Cs,"~",Ps)

    modelStr <- paste(modelStr,"\n",formula)
  }
  
  undirected <- graphStruct[duplicated(graphStruct$id) | duplicated(graphStruct$id, fromLast=TRUE),]
  undirected <- undirected[duplicated(undirected$id),] %>% mutate(formula = paste(from, "~~", to))

  for (x in undirected$formula) {

    modelStr <- paste(modelStr,"\n",x)
  }

  return(modelStr)
}


bn_to_sem <- function(graph, data) {
  # Takes a bnlearn graph and returns a lavaan SEM fit to data
  # graph: bnlearn bn object
  # data: data.frame of data to fit

  graphStruct <- get_graphStruct(graph)
  
  model <- make_lavaan_paths(graphStruct)
  fit <- sem(model=model, data=data)
  return(fit)
}

# # ---
# # uncomment below to test

# # A basic example with no other dependencies
# 
# testDF <- data.frame(from=c("A","B","C","B"),
#                      to = c("B","C","A","A"))
# testDF %>% make_lavaan_paths() %>% cat()
# 
# # # Some more complex examples with bnlearn
# 
# library(bnlearn)
# 
# set.seed(11)
# nodes <- LETTERS[1:8]
# rgraph <- random.graph(nodes, num = 1, method = "ordered")
# 
# # # make lavaan model syntax for model fitting 
# rgraph %>% get_graphStruct() %>% make_lavaan_paths() %>% cat
# 
# # # make lavaan model syntax for data simulation (random weights)
# rgraph %>% get_graphStruct() %>% make_weighted_lavaan_paths(c(5,5),c(5,5)) %>% cat
# 
# # # No problem with very large graphs
# rgraph.big <- random.graph(as.character(1:100), num = 1, method = "ordered", prob=.05)
# rgraph.big %>% get_graphStruct() %>% make_lavaan_paths() %>% cat
# 
# # # that graph is pretty big
# graphviz.plot(rgraph.big, highlight = NULL, layout = "dot",
#               shape = "circle", main = NULL, sub = NULL)