
# load data
# -----------------------------------------
# source("../../study1/vacc_import_data.R", chd)


# load packages
# ----------------------------------------
library(tidyverse)
library(BiDAG)
# BiDAG depends on 'graph' package, source: http://www.bioconductor.org/packages/release/bioc/html/graph.html
# BiDAG depends on 'RBGL' package, source: https://bioconductor.org/packages/release/bioc/html/RBGL.html
library(bnlearn)
library(HydeNet)
library(gtools)

# source("cog-model-graph-tools.R")
# source("cog-model-main.R")
# source("../../Scripts/gmod_tools.R")
# source("../../Scripts/bnTheoryBlacklist.R")


# transform data
# ----------------------------------------

rescale_beta <- function(x, lower, upper) {
  # rescales onto the open interval (0,1)
  # rescales over theoretical bounds of measurement, specified by "upper" and "lower"
  
  N <- length(x)
  res <- (x - lower) / (upper - lower)
  res <- (res * (N - 1) + .5) / N
  
  return(as.vector(res))
}

d_bn_scaled <- d_bn %>%
  mutate_all(function(x){rescale_beta(x,-3,3)})

## set the seed to make your partition reproductible
set.seed(123)
trainInd <- sample(seq_len(nrow(d_bn_scaled)), size = floor(nrow(d_bn_scaled)*.80))

train <- d_bn_scaled[trainInd, ]
test <- d_bn_scaled[-trainInd, ]


# generate blacklists
# ----------------------------------------
nodes <- c(
  "diseaseRare",
  "diseaseSevere",
  "hb",
  "infantImmLimCap",
  "infantImmWeak",
  "medSkept",
  "nat",
  "overpar",
  "parentExpert",
  "vaccDanger",
  "vaccEff",
  "vaccIntent",
  "vaccStrain",
  "vaccTox"
)

theoryBasedHierarchy <- data.frame(
  node = nodes,
  order = c(
    3,
    3,
    1,
    3,
    3,
    2,
    1,
    2,
    2,
    3,
    3,
    4,
    3,
    3
  )
)

intent_is_dv <- data.frame(
  node = nodes,
  order = c(
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    1,
    1
  )
)

abstract_are_parents <- data.frame(
  node = nodes,
  order = c(
    2,
    2,
    1,
    2,
    2,
    2,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2
  )
)

abstract_are_parents_intent_dv <- data.frame(
  node = nodes,
  order = c(
    2,
    2,
    1,
    2,
    2,
    2,
    1,
    2,
    2,
    2,
    2,
    3,
    2,
    2
  )
)


make_adjacency_blacklist <- function(ordering_df){
  library(dplyr)
  
  blacklist <- make_theory_blacklist(ordering_df)
  blacklist <- suppressWarnings(bind_rows(blacklist, data.frame(from = nodes, to = nodes))) %>% arrange(from)
  adjBlacklist <- igraph::graph.data.frame(blacklist)
  adjBlacklist <- igraph::get.adjacency(adjBlacklist, sparse=FALSE)
  
  # adjBlacklist <- abs(1-adjBlacklist)
  
  return(adjBlacklist)
}

bl_intent_is_dv <- make_adjacency_blacklist(intent_is_dv)
bl_theoryBasedHierarchy <- make_adjacency_blacklist(theoryBasedHierarchy)
bl_abstract_are_parents <- make_adjacency_blacklist(abstract_are_parents)
bl_abstract_are_parents_intent_dv <- make_adjacency_blacklist(abstract_are_parents_intent_dv)


# identify seachspace
# ----------------------------------------
set.seed(123)
ss_intent_is_dv <- find_startspace(train, alpha =.4, blacklist = bl_intent_is_dv)
ss_theoryBasedHierarchy <- find_startspace(train, alpha =.4, blacklist = bl_theoryBasedHierarchy)
ss_abstract_are_parents <- find_startspace(train, alpha =.4, blacklist = bl_abstract_are_parents)
ss_abstract_are_parents_intent_dv <- find_startspace(train, alpha =.4, blacklist = bl_abstract_are_parents_intent_dv)
ss_unconstrained <- find_startspace(train, alpha =.4)

check_scoretable_size(ss_intent_is_dv)
check_scoretable_size(ss_theoryBasedHierarchy)
check_scoretable_size(ss_abstract_are_parents)
check_scoretable_size(ss_abstract_are_parents_intent_dv)
check_scoretable_size(ss_unconstrained)

# calculate scoretables and do MCMC
# ----------------------------------------
nCores <- parallel::detectCores(logical=TRUE)

startspaces <- list(ss_unconstrained,
                    ss_intent_is_dv, 
                    ss_abstract_are_parents, 
                    ss_abstract_are_parents_intent_dv, 
                    ss_theoryBasedHierarchy)

blacklists <- list(NULL,
                   bl_intent_is_dv, 
                   bl_abstract_are_parents, 
                   bl_abstract_are_parents_intent_dv, 
                   bl_theoryBasedHierarchy)

fsuffixes <- c("unconstrained",
               "intent_is_dv", 
               "abstract_are_parents",
               "abstract_are_parents_intent_dv", 
               "theoryBasedHierarchy")


if (!exists("compute_scoretables")) {
  compute_scoretables <- TRUE
}

if (compute_scoretables) {
  for (i in 1:length(startspaces)) {
    t1 <- Sys.time()
    scoretable <- generate_scoretable(startspaces[[i]], train, cores = nCores)
    print(Sys.time() - t1)
    saveRDS(scoretable, paste0("scoretable-", fsuffixes[i], ".rds"))
  }
  
}


# do this separately, it's comparativly fast and less parallelized
# all_results <- c()
# 
# for (i in 1:length(startspaces)) {
#   all_results[i] <- parallel_orderMCMC(train, 
#                                        startspaces[[i]], 
#                                        scoretables[[i]], 
#                                        1e6, 
#                                        100, 
#                                        blacklist = blacklists[[i]], 
#                                        chains = 4, 
#                                        cores = 4
#   )
#   
# }
# 
# saveRDS(all_results, "allresults.rds")