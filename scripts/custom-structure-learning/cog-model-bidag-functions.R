
## ----------------------------------------
## Generating BiDAG scoretable 
## -------------------------------------------
## Code below is modified from BiDAG package v1.2 under GPL 2/3 license 
## 8/21/18, 1:11 PM
library(parallel)

c_definestartspace <- function(alpha, param, cpdag = FALSE, n, algo = "pc") {
  # condense, just error out for other score types 
  if (is.null(alpha)) {
    if (n < 50) {
      alpha <- 0.4
    }
    else {
      alpha <- max(20/n, 0.01)
    }
  }
  if (param$type == "bde") {
    if (cpdag) {
      pc.skel <- pc(suffStat = list(d1 = param$d1, d0 = param$d0, 
                                    data = param$data), indepTest = weightedbinCItest, 
                    alpha = alpha, labels = colnames(param$data), 
                    verbose = FALSE)
    }
    else {
      pc.skel <- pcalg::skeleton(suffStat = list(d1 = param$d1, 
                                                 d0 = param$d0, data = param$data), indepTest = weightedbinCItest, 
                                 alpha = alpha, labels = colnames(param$data), 
                                 verbose = FALSE)
    }
  }
  else if (param$type != "bde") {
    if (is.null(param$weightvector)) {
      cormat <- cor(param$data)
      N <- nrow(param$data)
    }
    else {
      N <- sum(param$weightvector)
      cormat <- cov.wt(param$data, wt = param$weightvector, 
                       cor = TRUE)$cor
    }
    if (cpdag) {
      pc.skel <- pcalg::pc(suffStat = list(C = cormat, 
                                           n = N), indepTest = pcalg::gaussCItest, alpha = alpha, 
                           labels = colnames(param$data), skel.method = "stable", 
                           verbose = FALSE)
    }
    else {
      pc.skel <- pcalg::skeleton(suffStat = list(C = cormat, 
                                                 n = N), indepTest = pcalg::gaussCItest, alpha = alpha, 
                                 labels = colnames(param$data), method = "stable", 
                                 verbose = FALSE)
    }
  }
  g <- pc.skel@graph
  startspace <- 1 * (dag2adjacencymatrix(g))
}


c_listpossibleparents.PC.aliases<-function(skeletonedges,isgraphNEL=FALSE,n,updatenodes=c(1:n)){
  if(isgraphNEL==FALSE){
    l1<-ncol(skeletonedges)
  } else {l1<-length(skeletonedges)}
  listy<-vector("list",l1)
  aliases<-vector("list",l1)
  numparents<-vector("numeric",l1)
  
  #we keep record of which parent table lengths we already constructed
  table.with.k.parents<-matrix(rep(0,l1*2),nrow=2,ncol=l1)
  
  for (i in updatenodes){
    if (isgraphNEL==TRUE) {possparents<-skeletonedges[[i]]$edges
    } else{possparents<-which(skeletonedges[,i]==1)}
    aliases[[i]]<-possparents
    l<-length(possparents)
    numparents[i]<-l
    possparents<-c(1:l)
    if (l==0){
      matrixofparents<-matrix(rep(NA,1),1,1)
    } else if (table.with.k.parents[1,l]>0){
      matrixofparents<-listy[[table.with.k.parents[2,l]]]
    } else {
      matrixofparents<-rep(NA,l)
      for (r in 1:l){
        combpossparents<-combinations(l,r,possparents)
        if(r<l){
          for (j in 1:(l-r)){
            combpossparents <- cbind(combpossparents, NA)
          }
        }
        matrixofparents<-rbind(matrixofparents,combpossparents,deparse.level=0)
      }
    }
    listy[[i]] <- matrixofparents
    table.with.k.parents[1,l]<-1
    table.with.k.parents[2,l]<-i
  }
  listz<-list()
  listz$parenttable<-listy
  listz$aliases<-aliases
  listz$numparents<-numparents
  listz$numberofparentsvec<-lapply(numparents,function(x)rep(c(0:x),choose(x,c(0:x))))
  
  return(listz)
}


c_scorepossibleparents.alias <- function(parenttable, aliases, n, param, updatenodes=c(1:n), parentmaps=NULL, numparents=NULL, numberofparentsvec=NULL, cores = NULL) {
  if (is.null(cores)) {
    cores <- detectCores()
  }

  listz <- mclapply(updatenodes, function(i) {
    scoretemp <- c_TableDAGscore.alias(parenttable[[i]], i, n, aliases[[i]], param, parentmaps[[i]], numparents[i], numberofparentsvec[[i]])
    return(as.matrix(scoretemp))
    # listz[[i]] <- as.matrix(scoretemp)
  },
  mc.preschedule = FALSE,
  mc.cores = cores
  )

  return(listz)
}


c_TableDAGscore.alias <- function(parentrows, j, n, alias, param, parentmaps = NULL, numparents = NULL, numberofparentsvec = NULL){
  # condense, just error out for other score types 
  if (param$type == "bge" | param$type == "custom") { 
    nrows <- nrow(parentrows)
    P_local <- numeric(nrows)
    for (i in 1:nrows) {
      parentnodes <<- alias[parentrows[i, !is.na(parentrows[i, # modified 8/21/18, 2:05 PM
                                                            ])]]
      P_local[i] <- c_DAGcorescore(j, parentnodes, n, param)
    }
  }
  else {
    nrows <- nrow(parentrows)
    parentnodes <- alias[parentrows[nrows, !is.na(parentrows[nrows, 
                                                             ])]]
    P_local <- DAGbinarytablescore(j, parentnodes, n, param, 
                                   parentrows, parentmaps, numparents, numberofparentsvec)
  }
  return(P_local)
}


c_DAGcorescore <- function(j, parentnodes, n, param) {
  # j_export <<- j
  # parentnodes_export <<- parentnodes
  # n_export <<- n
  
  d <- param$data
  node_names <- colnames(d)
  child_node <- node_names[j]
  parent_nodes <- node_names[parentnodes]
  # mychild <<- child_node
  # myparents <<- parent_nodes
  
  # f <- make_family_formula(child_node, parent_nodes)
  # m <- glm(f, data = d)
  # corescore <- as.numeric(LogLik(m))
  
  # corescore <- tryCatch(
  #   score_family(child_node, parent_nodes, d, method="BFGS"),
  #   error = function(e){ score_family(child_node, parent_nodes, d) }
  # )
  
  corescore <- score_family(child_node, parent_nodes, d) # stick with nelder-mead
  
  return(corescore)
}

source("cog-model-iterativeMCMC.R")


reinit_BiDAG <- function(){
  
  if("BiDAG" %in% (.packages())){
    detach("package:BiDAG", unload=TRUE)
  }
  
  library(BiDAG)
}


customize_BiDAG <- function(){
  if( !("BiDAG" %in% (.packages())) ){
    library(BiDAG)
  }
  # assignInNamespace("DAGcorescore", c_DAGcorescore, "BiDAG")
  # assignInNamespace("TableDAGscore.alias", c_TableDAGscore.alias, "BiDAG") 
  assignInNamespace("iterativeMCMCplus1", c_iterativeMCMCplus1, "BiDAG") 
  
}


find_startspace <- function(data, alpha = .4, blacklist = NULL, prune=TRUE){
  nodes <- colnames(data)
  n <- length(nodes)
  
  if (prune){
    startspace_obj <- BiDAG::iterativeMCMCsearch(
      ncol(data),
      scoreparam = BiDAG::scoreparameters(n, scoretype="bge", data),
      blacklist = blacklist,
      scoreout = TRUE,
      alpha = .4
    )
    
    startspace <- startspace_obj$space$adjacency
    
  } else {

    # create starting search space adjacency matrix
    startspace <- matrix(rep(1,n^2), ncol=n)
    rownames(startspace) <- nodes
    colnames(startspace) <- nodes
    diag(startspace) <- 0
  }

  return(startspace)

}

generate_scoretable <- function(startspace, data, cores=NULL){
  # for simple graph, allow all possible edges in space
  library(gtools) # needed for combinations()
  
  n <- ncol(startspace)
  ptab <- BiDAG:::listpossibleparents.PC.aliases(startspace,isgraphNEL=FALSE,n)
  param <- scoreparameters(n, "bge", data)
  
  parenttable<-ptab$parenttable
  aliases<-ptab$aliases
  numberofparentsvec<-ptab$numberofparentsvec
  numparents<-ptab$numparents
  updatenodes<-c(1:n)

  rowmapsneeded<-BiDAG:::parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
  
  scoretable <- c_scorepossibleparents.alias(parenttable,aliases,n,param,updatenodes,rowmaps,numparents,numberofparentsvec, cores)
}


check_scoretable_size <- function(startspace){
  
  ptab <- c_listpossibleparents.PC.aliases(startspace, isgraphNEL=F,nrow(startspace))
  length(unlist(ptab$parenttable))
}

parallel_orderMCMC <- function(data, startspace, scoretable, iterations = 1e4, stepsave=NULL, blacklist = NULL, chains = 2, cores = NULL, ...) {
  
  if (is.null(cores)) {
    cores <- detectCores()
  }
  
  n <- ncol(data)
  # orders <- rep(list(sample(1:n)), chains)
  orders <- purrr::map(1:chains, function(x){sample(1:n)})
  
  chains <- parallel::mclapply(orders, function(order) {
    orderMCMC(
      n,
      startspace = startspace,
      startorder = order,
      plus1 = FALSE, # orderMCMC only
      scoreparam = scoreparameters(n, scoretype = "bge", data),
      scoretable = scoretable,
      blacklist = blacklist, # works on orderMCMC
      chainout = TRUE, # only orderMCMC
      scoreout = TRUE, # only orderMCMC
      # verbose = TRUE
      iterations = iterations,
      stepsave = stepsave,
      ...
    )
  },
  mc.cores = cores,
  mc.set.seed = TRUE
  )

  
  return(chains)
}

merge_chains <- function(chains_list, burn=.20) {
  n_iter <- length(chains_list[[1]]$chain$DAGscores)
  n_burnin <- round(n_iter*burn) + 1
  merged_chain <- list(DAGscores = unlist( map(chains_list, function(x){x$chain$DAGscores[n_burnin:n_iter]}) ),
                       incidence = unlist( map(chains_list, function(x){x$chain$incidence[n_burnin:n_iter]}), recursive=FALSE)
  )
}


extract_mapdag <- function(chain){
  map_score <- max(chain$DAGscores)
  map_index <- match(map_score, chain$DAGscores)
  
  map_adj <- chain$incidence[map_index]
  
  map_dag <- adjacency2dag(map_adj[[1]], nodes = nodes)
  map_dag <- as.bn(map_dag)
}

BiDAG_traceplot <- function(chains, hide = 200, show_chains = NULL){
  
  n <- length(chains)
  
  if (is.null(show_chains)){
    show_chains = 1:n
  }
  
  scoreChains <- lapply(chains, function(x){unlist(x$chain$DAGscores)})
  
  chaindf <- data.frame( 
    # value = c(scoreChains[[1]], scoreChains[[2]]),
    value = unlist(scoreChains),
    chain = rep(1:n, each = length(scoreChains[[1]])),
    iter = rep(seq(1:length(scoreChains[[1]])),2)
  )
  
  
  chaindf %>%
    filter(iter > hide, chain %in% show_chains) %>%
    ggplot(aes(x=iter, y = value, color=as.factor(chain))) +
    # ggplot(aes(x=iter, y = value)) +
    geom_line(alpha=.9) +
    labs(color = "chain", x = "iteration", y = "score") +
    # facet_wrap(~chain, ncol=1) +
    NULL
}
