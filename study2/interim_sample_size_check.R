setwd("~/Projects/Vaccines/Repos/vaccineBeliefs/study3/")

read_qualtrics_csv <- function(fname) {
  df <- read.csv(fname, skip = 3, header = F)
  headers <- as.matrix(read.csv(fname, skip = 0, header = F, nrows = 1, as.is = T))
  colnames(df) <- headers
  df <- df[which(df[, "DistributionChannel", ] == "anonymous"), ] # remove survey previews
  total_n <- nrow(df)
  
  remove_cols <- c(
    "RecipientLastdata",
    "RecipientFirstdata",
    "RecipientEmail",
    "Finished",
    "ResponseId",
    "ExternalReference",
    "DistributionChannel",
    "UserLanguage",
    "Status"
  )
  
  df <- df[, -which(colnames(df) %in% remove_cols)]
  
  return(df)
}

pre <- read_qualtrics_csv("./data/vaccine+many+beliefs+-+interv+1+pre+-+APRIL_April+24%2C+2018_10.54.csv")
post1 <- read_qualtrics_csv("./data/vaccine+many+beliefs+-+interv+1+post+g1+-+APRIL_April+24%2C+2018_10.55.csv") %>%
  select(Progress, hb_check_5, vaccEff_check_1,vaccTox_check_6, vaccDanger_check_3)
post2 <- read_qualtrics_csv("./data/vaccine+many+beliefs+-+interv+1+post+g2+-+APRIL_April+24%2C+2018_10.55.csv") %>%
  select(Progress, hb_check_5, vaccEff_check_1,vaccTox_check_6, vaccDanger_check_3)
post3 <- read_qualtrics_csv("./data/vaccine+many+beliefs+-+interv+1+post+g3+-+APRIL_April+24%2C+2018_10.56.csv") %>%
  select(Progress, hb_check_5, vaccEff_check_1,vaccTox_check_6, vaccDanger_check_3)

df_post <- post1 %>%
  bind_rows(post2) %>%
  bind_rows(post3) %>%
  as_tibble() %>%
  filter(Progress==100) %>%
  filter(hb_check_5==5,
         vaccEff_check_1==1,
         vaccTox_check_6==6,
         vaccDanger_check_3==3) 

# started with 667, invited 574 to return, got 450 to return and complete, 405 of them passed all checks
# so lost 14%, 22%, and 10% at each step.