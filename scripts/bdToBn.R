# author: Kara Weisman
# created: 11/10/17, 2:28 PM
# last updated: 11/10/17, 2:28 PM

# script for translating BDgraph posterior probabilities to bnlearn

library(tidyverse)
library(bnlearn)
library(BDgraph)
library(tidygraph)

# make model with BDgraph -----
fit.bd.ggm <- bdgraph(d_bn, method = "ggm")
plot(fit.bd.ggm)

# approach #1: from last graph -----

## get links
links.bd.ggm <- fit.bd.ggm$last_graph

# ## format as arcs for bnlearn
# arcs.whitelist <- links.bd.ggm %>%
#   data.frame() %>%
#   rownames_to_column("from") %>%
#   gather(to, yes_no, -from) %>%
#   filter(yes_no == 1) %>%
#   data.frame() %>%
#   dplyr::select(from, to)

## make the blacklist for bnlearn
arcs.blacklist <- links.bd.ggm %>%
  data.frame() %>%
  rownames_to_column("from") %>%
  gather(to, yes_no, -from) %>%
  filter(yes_no == 0) %>%
  data.frame() %>%
  dplyr::select(from, to)

## make model with bnlearn
fit.mmhc.blacklist <- mmhc(d_bn, blacklist = arcs.blacklist)
plot(fit.mmhc.blacklist)

# approach #2: from posterior probabilities -----

## format as arcs for bnlearn
arcs.whitelist995 <- fit.bd.ggm$p_links %>%
  data.frame() %>%
  rownames_to_column("from") %>%
  gather(to, posterior, -from) %>%
  filter(posterior >= 0.995) %>%
  data.frame() %>%
  dplyr::select(from, to)

## make the blacklist for bnlearn
arcs.blacklist995 <- fit.bd.ggm$p_links %>%
  data.frame() %>%
  rownames_to_column("from") %>%
  gather(to, posterior, -from) %>%
  filter(posterior < 0.995) %>%
  data.frame() %>%
  dplyr::select(from, to)

## make model with bnlearn
fit.mmhc.blacklist995 <- mmhc(d_bn, blacklist = arcs.blacklist995)
plot(fit.mmhc.blacklist995)

# convert to tidy -----
tidy_graph <- as_tbl_graph(fit.mmhc.blacklist995$arcs)
