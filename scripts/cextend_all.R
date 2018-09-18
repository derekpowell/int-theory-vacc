# Author: Derek Powell
# Created: 12/20/17, 3:34 PM
# Last updated: 12/20/17, 3:34 PM

tag_arcs <- function(graphStruct) {
  g <- graphStruct %>% 
    mutate(from =  as.character(from)) %>%
    mutate(to = as.character(to)) %>%
    mutate(id = map2_chr(from, to, cat_and_sort))
  return(as_tibble(g))
}

cextend_all <- function(cpDAG) {
  # Generate all consistent extensions of equivalence class
  # cpDAG: complete partially directed graph (cpdag(bnobject))
  # returns: a list of bn objects
  library(gtools)
  undirected <- undirected.arcs(cpDAG) %>% as_tibble() %>% tag_arcs()
  uniqueArcs <- unique(undirected$id)
  
  arcDF <- undirected %>% distinct(id, .keep_all=TRUE)
  arcNum <- nrow(arcDF)
  
  # m <- as.data.frame(matrix(rbinom(arcNum*2, 1, 0.5),ncol=arcNum))
  
  grid <- gtools::permutations(2,arcNum,v=c(0,1),repeats.allowed=TRUE)
  # permutations(2,i,v=c(0,1),repeats.allowed=TRUE)
  # grid <- expand.grid(m)
  
  graphList <- list()
  for (i in 1:dim(grid)[1]) {
    # for (i in 1:50) {
    gDF <- arcDF %>% mutate(direction = as.numeric(grid[i,]))
    g <- cpDAG
    for (j in 1:nrow(gDF)) {
      arc <- gDF[j,]
      if (arc$direction==1){
        g <- set.arc(g,arc$from,arc$to, check.cycles=FALSE, check.illegal=FALSE)
        
      } else if (arc$direction==0) {
        g <- set.arc(g, arc$to, arc$from, check.cycles=FALSE, check.illegal=FALSE)
      }
    }
    checkDag <- acyclic(g)
    
    if (checkDag==TRUE) {
      # graphList <- c(graphList, g)
      graphList[[length(graphList)+1]] <- g
    }
    
  }
  return(graphList)
}