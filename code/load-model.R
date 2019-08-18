# first load in mcmc results
all_results <- readRDS("../local/allresults.rds")
nodes <- colnames(train)

# then grab map dag for theory-based constraints
merged_chains <- merge_chains(all_results[[5]]) # note: this will need to be updated! [fix needed]
map_dag <- extract_mapdag(merged_chains)

arcs_df <- bnlearn_to_df(map_dag)

df_to_list <- function(arcs_df, as.strings=FALSE) {
  # take arcs_df and make list of formulas
  arcs_df$to <- as.character(arcs_df$to)
  arcs_df$from <- as.character(arcs_df$from)
  
  nodes <- unique(unlist(arcs_df))
  output <- lapply(nodes, function(node) {
    edges <- arcs_df[which(arcs_df$to == node), "from"]
    if (length(edges) < 1) {
      edges <- NULL
    }
    list(child=node, parents=edges)
    
  })
  names(output) <- nodes
  return(output)
}

graph_list <- df_to_list(arcs_df)

bagOfModels <- lapply(graph_list, function(x){fit_node_cpt(x, train)})