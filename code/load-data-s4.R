library(tidyverse)
library(lubridate)

read_qualtrics <- function(fname){
  # load raw qualtrics data
  header <- colnames(read.csv(fname, header = TRUE))
  df <- read.csv(fname, skip = 3, header = FALSE, col.names = header)
  
  # do some preliminary tidying
  df <- df %>%
    as_tibble()
  
  return(df)
}

df4 <- read_qualtrics("../data/vaccine+many+beliefs+-+may2019_May+21,+2019_16.12.csv")

df4 <- df4 %>%
  filter(!duplicated(workerId)) %>%
  filter(check_3==3, check_5==5) %>%
  mutate(StartDate = ymd_hms(StartDate)) %>%
  rename(Duration = Duration..in.seconds.) %>%
  select(workerId, matches(".*(_[0-9])"), vacc_news, vacc_news_freq, 
         news_time, news_favorable, vacc_conv, conv_favorable, outbreaks, outbreak_favorable, other,
         sex, age, race, religion, educ, income, 
         pol_party, pol_ind_leaning, parent, expecting, children,
         youngest_child, refuse, comments, StartDate,
         -check_3, -check_5
         ) %>%
  mutate_at(vars(matches(".*(_[0-9]*r)"), -pol_party_4_TEXT, -pol_rating_1), funs(8 - .)) %>% # reverse code
  gather(item, resp, matches(".*(_[0-9])"), -pol_party_4_TEXT, -pol_rating_1) %>%
  mutate(scale = gsub("(_([0-9])+(r)*)","",item)) %>%
  mutate(resp = as.numeric(resp))

