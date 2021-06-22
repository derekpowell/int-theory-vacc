custom_plot.HydeNetwork <- function(x, 
                             customNodes = NULL,
                             customEdges = NULL,
                             ..., 
                             removeDeterm = FALSE,
                             useHydeDefaults = TRUE)
{
  if (removeDeterm) x <- plot_nondeterm_only(x)
  
  node_df <- 
    DiagrammeR::create_node_df(n = length(x[["nodes"]]),
                               label = x[["nodes"]])
  # 
  # node_df <- data.frame(nodes = x[["nodes"]],
  #                       stringsAsFactors = FALSE)
  if (useHydeDefaults) node_df <- HydeNet:::mergeDefaultPlotOpts(x, node_df)
  
  if (!is.null(customNodes)) node_df <- HydeNet::mergeCustomNodes(node_df, customNodes)
  
  edge_table <- do.call("rbind", 
                        mapply(FUN = HydeNet:::mapEdges, 
                               x[["nodes"]], 
                               x[["parents"]],
                               MoreArgs = list(node_df = node_df)))
  
  edge_df <- DiagrammeR::create_edge_df(from = edge_table[, 2], 
                                        to = edge_table[, 1])
  
  if (!is.null(customEdges)) HydeNet::mergeCustomEdges(edge_df, customEdges)
  
  
  
  DiagrammeR::create_graph(nodes_df = node_df,
                           edges_df = edge_df,
                           ...
                           )
  
}


## ----------------------------------------

# code to demo. For some reason set_edge_attrs() isn't working right, so have to do it manually
# z <- custom_plot.HydeNetwork(h)
# z$edges_df <- z$edges_df %>% mutate(color = "blue", penwidth=rep(c(1,2,3),10))
# render_graph(z)
