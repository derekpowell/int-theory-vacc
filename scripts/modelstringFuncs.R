# Functions for translating network objects across packages via model strings / formulas
# Author: Derek Powell
# Created: 3/8/18, 3:46 PM
# Updated: 3/8/18, 4:24 PM

wrap_brackets <- function(string) {
  return(paste0("[",string,"]"))
}

wrap_space <- function(string) {
  return(paste0(" ",string," "))
}

make_bnlearn_node_string <- function(node, children) {
  output <- as.character(node)
  
  if (length(children) > 0) {
    output <- paste0(output,"|")
    for (c in children) {
      output <- paste0(output,c,":") 
      
    }
    output <- substr(output,1,nchar(output)-1)
  }
  
  return(wrap_brackets(output))
}

make_hyde_node_string <- function(node, children) {
  output <- as.character(node)
  
  if (length(children) > 0) {
    output <- paste0(output,"|")
    for (c in children) {
      output <- paste0(output,c,"*") 
      
    }
    output <- substr(output,1,nchar(output)-1)
  }
  
  return(wrap_space(output))
}

df_to_hyde_string <- function(df) {
  output <- ""
  for (node in unique(unlist(df))) {
    edges <- df[which(df$to==node),"from"]
    edges <- edges$from
    nodeString <- make_hyde_node_string(node, edges)
    output <- paste0(output,"+",nodeString)
    
  }
  output <- paste0("~",substr(output, 2, nchar(output)))
  
  return(output)
}


df_to_bnlearn_string <- function(df) {
  output <- ""
  for (node in unique(unlist(df))) {
    edges <- df[which(df$to==node),"from"]
    edges <- edges$from
    nodeString <- make_node_string(node, edges)
    output <- paste0(output,nodeString)
    
  }
  return(output)
}

bnlearn_to_hyde <- function(bn.net) {
  arcs <- as.data.frame(bnlearn::arcs(bn.net))
  return(df_to_hyde_string(arcs))
}

bn_to_dbn <- function(model) {
  # takes a bnlearn bn, bn.fit, or modelstring and transforms it into a dynamic bayes net
  # with a second mirrored stage (could later be generalized to add a timestep)
  require(dplyr)
  
  if (class(model)=="character") {
    net <- model2network(model)
  }
  else { net <- model}
  
  edges <- arcs(net)
  t1 <- edges %>% as_tibble() %>% mutate(from=paste0(from,"_1"), to=paste0(to,"_1"))
  t2 <- edges %>% as_tibble() %>% mutate(from=paste0(from,"_2"), to=paste0(to,"_2"))
  nodes <- nodes(net)
  bridge <- data.frame(from=nodes,to=nodes) %>% mutate(from=paste0(from,"_1"), to=paste0(to,"_2"))
  
  edges <- bind_rows(t1,t2) %>% bind_rows(bridge)
  
  # return(df_to_bnlearn_string(edges))
  return(edges)
}
