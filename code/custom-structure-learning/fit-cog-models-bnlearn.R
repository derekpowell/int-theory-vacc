
# load data
# -----------------------------------------
source("../load-data-s1.R", chdir=TRUE)

# load packages and functions
# ----------------------------------------
library(tidyverse)
library(bnlearn)
library(HydeNet)
library(parallel)

source("cog-model-main.R")
source("../graph-model-tools.R")

# transform data
# ----------------------------------------
d_bn_scaled <- d_bn %>%
  mutate_all(function(x){rescale_beta(x,-3,3)})

## set the seed to make your partition reproducible
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

theory_bl <- tiers2blacklist(
  list(
    c("hb","nat"),
    c("medSkept", "overpar", "parentExpert"),  # idk if overparenting makes sense here
    c("diseaseRare","diseaseSevere","vaccEff","vaccTox","vaccStrain",
      "vaccDanger","infantImmLimCap","infantImmWeak"), 
    c("vaccIntent") 
    )
  )

intent_dv_bl <- tiers2blacklist(
  list(
    c(  "diseaseRare",
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
        "vaccStrain",
        "vaccTox"),
    c("vaccIntent")
  )
)


intent_dv_abstract_parents_bl <- tiers2blacklist(
  list(
    c("nat","hb"),
    c(  "diseaseRare",
        "diseaseSevere",
        "infantImmLimCap",
        "infantImmWeak",
        "medSkept",
        "overpar",
        "parentExpert",
        "vaccDanger",
        "vaccEff",
        "vaccStrain",
        "vaccTox"),
    c("vaccIntent")
  )
)

## ------- cross validation
## set up different model algos to compare

hc_args <- list(score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000)

models <- list(
  hc = partial(hc, score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000),
  hc_intentdv = partial(hc, blacklist=intent_dv_bl, score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000),
  hc_idv_abspar = partial(hc, blacklist=intent_dv_abstract_parents_bl, score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000),
  hc_theory = partial(hc, blacklist=theory_bl, score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000),
  mmhc_05 = partial(mmhc, restrict.args=list(alpha=.05), maximize.args = hc_args),
  mmhc_01 = partial(mmhc, restrict.args=list(alpha=.01), maximize.args = hc_args),
  mmhc_intentdv_05 = partial(mmhc, blacklist = intent_dv_bl, restrict.args=list(alpha=.05), maximize.args = hc_args),
  mmhc_intentdv_01 = partial(mmhc, blacklist = intent_dv_bl, restrict.args=list(alpha=.01), maximize.args = hc_args),
  mmhc_idv_abspar_05 = partial(mmhc, blacklist = intent_dv_abstract_parents_bl, restrict.args=list(alpha=.05), maximize.args = hc_args),
  mmhc_idv_abspar_01 = partial(mmhc, blacklist = intent_dv_abstract_parents_bl, restrict.args=list(alpha=.01), maximize.args = hc_args),
  mmhc_theory_05 = partial(mmhc, blacklist = theory_bl, restrict.args=list(alpha=.05), maximize.args = hc_args),
  mmhc_theory_01 = partial(mmhc, blacklist = theory_bl, restrict.args=list(alpha=.01), maximize.args = hc_args)
)

do_cross_validation <- function(fold, data=NULL, model=NULL){
  train_fold <- data[fold,]
  test_fold <- data[-fold,]
  network <- model(train_fold)
  oos_score <- score(network, test_fold, type="custom", fun=cog_model_score_func, args=list(train_data=train_fold))
  
  return(
    oos_score
  )
}

## split training data into k-folds
set.seed(12345)
cv <- caret::createMultiFolds(train$vaccIntent, k=10, times=10)

print("beginning cross-validation ...")

res_df_all <- tibble() # init table to hold results

## loop through models, doing fitting over runs with parallel processing
for (i in 1:length(models)) {
  
  model_name <- names(models)[i]
  model <- models[[model_name]]
  
  cv_res <- mclapply(
    cv, 
    do_cross_validation, 
    data = train, 
    model = model, 
    mc.cores = detectCores(logical=TRUE)
  )
  
  res_df_temp <- as_tibble(cv_res) %>% 
    gather(run, loss) %>% 
    mutate(algo = model_name)
  
  res_df_all <- bind_rows(res_df_all, res_df_temp)
  
  print(paste(model_name, "completed"))
}

saveRDS(res_df_all, file="../../local/cv-res.rds")

print("all completed, output saved")


## as a first pass just do mmhc to see if blacklist is being respected




# cog_model_score_func <- function(node, parents, data, args=NULL){
#   return(score_family(node, parents, data))
# }

# mmhc_cv <- purrr::partial(mmhc, 
#                           maximize.args = hc_args
# )
# 
# mmhc_theory_cv <- purrr::partial(mmhc, 
#                                  blacklist = theory_bl,
#                                  maximize.args = hc_args
# )
# 
# mmhc_intent_dv_cv <- partial(mmhc, blacklist = intent_dv_bl, maximize.args = hc_args)



## going to have to implement my own CV functions and loss functions for held-out data
## the jags thing is gonna be too slow, so gotta do it family-by-family ala scoring

# already have the start of this with find_node_coefs(), just need to create a score from that
## -------
# 
# bntest3 <- mmhc(
#   train, 
#   blacklist = theory_bl,
#   maximize.args = list(score="custom", fun=cog_model_score_func, maxp=5, max.iter=1000)
#   )
# 
# # can get score but it is finding new parameters given the structure
# # could pass in the training data in args and estimate w/ that rather than
# # the original data being passed
# # score(bntest3, test, type="custom", fun=cog_model_score_func, args=NULL)
# 
# # cog_model_score_func("vaccIntent",c("diseaseRare","vaccEff"), test)
# 
# 
# # testing .. it works!
# h <- HydeNetwork(as.formula(bnlearn_to_hyde_string(bntest3)))
# plot(h)
# 
# model_string <- suppressWarnings(write_jags_model(bntest3, train))
# cat(model_string)
