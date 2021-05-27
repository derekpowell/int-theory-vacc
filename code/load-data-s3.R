library(tidyverse)

read_qualtrics <- function(fname){
  # load raw qualtrics data
  header <- colnames(read.csv(fname, header = TRUE))
  df <- read.csv(fname, skip = 3, header = FALSE, col.names = header)
  
  # do some preliminary tidying
  df <- df %>%
    as_tibble()
  
  return(df)
}





# pre <- read_qualtrics("../data/vacc+tox+interv+full+-+pre_April+3,+2019_15.38.csv")
# post <- read_qualtrics("../data/vacc+tox+interv+full+-+post_April+5,+2019_13.59.csv")

pre_a <- read_qualtrics("../data/vacc+tox+interv+full+-+pre_April+3,+2019_15.38.csv") %>%
  mutate(study="A")
post_a <- read_qualtrics("../data/vacc+tox+interv+full+-+post_April+5,+2019_13.59.csv") %>%
  mutate(study="A")

pre_b <- read_qualtrics("../data/vacc+tox+interv+full2+-+pre_June+6,+2019_09.35.csv") %>%
  mutate(study="B")
post_b <- read_qualtrics("../data/vacc+tox+interv+full2+-+post_June+7,+2019_18.00.csv") %>%
  mutate(study="B")

pre <- bind_rows(pre_a, pre_b)
post <- bind_rows(post_a, post_b)

summarize <- dplyr::summarise

df_pre <- pre %>%
  filter(workerId!="e082cf59a57d4562") %>% # hashed manually
  filter(Progress==100, check_3==3, check_5==5) %>%
  mutate(parent = recode(parent, `23`=1, `24`=0)) %>%
  # select(-check_3, -check_5) %>% 
  select(workerId, matches(".*(_[0-9])"), -check_3, -check_5, parent,
         sex, age, race, religion, educ, income, expecting, children, 
         youngest_child, refuse, return, StartDate, study) %>%
  rename(StartDate_pre = StartDate) %>% 
  gather(item, resp, matches(".*(_[0-9])")) %>%
  mutate(item = paste0("pre_",item)) %>%
  spread(item, resp) %>% 
  mutate(
    sex = case_when(
      sex==1 ~ "male",
      sex==0 ~ "female",
      sex==-1 ~ "other or prefer not to say"
    ),
    race = case_when(
      race == 1 ~ "White",
      race == 2 ~ "Hispanic/Latino",
      race == 3 ~ "Black or African American",
      race == 4 ~ "Hawaiian or PI",
      race == 5 ~ "Asian",
      race == 6 ~ "Native American",
      race == 7 ~ "Other or prefer not to say"
    ),
    religion = case_when(
      religion == 1 ~ "Christian",
      religion == 2 ~ "Muslim",
      religion == 3 ~ "Jewish",
      religion == 4 ~ "Hindu",
      religion == 5 ~ "Buddhist",
      religion == 6 ~ "Non-religious",
      religion == 7 ~ "Other or prefer not to say"
    ),
    educ = case_when(
      educ == 1 ~ "No Diploma",
      educ == 2 ~ "High School",
      educ == 3 ~ "Some Undergraduate",
      educ == 4 ~ "Undergraduate Degree",
      educ == 5 ~ "Some Graduate",
      educ == 6 ~ "Graduate Degree",
      educ == 7 ~ "Doctorate",
      educ == 8 ~ "Prefer not to say"
    ),
    income = case_when(
      income == 1 ~ "Less than $20,000",
      income == 2 ~ "$20,000 - $30,000",
      income == 3 ~ "$30,001 - $50,000",
      income == 4 ~ "$50,001 - $70,000",
      income == 5 ~ "$70,001 - $100,000",
      income == 6 ~ "More than $100,000",
      income == 0 ~ "Prefer not to say",
    ),
    expecting = case_when(
      expecting == 23 ~ "yes", #F@!&ing qualtrics
      expecting == 24 ~ "no"
    ),
    parent = ifelse(parent==1, "Yes","No")
  )

df_post <- post %>%
  filter(Progress==100, check_3==3, check_5==5) %>%
  select(workerId, matches(".*(_[0-9])"), -check_3, -check_5,
         condition, interv_t_Page.Submit, StartDate) %>%
  rename(read_time = interv_t_Page.Submit) %>%
  gather(item, resp, -workerId, -condition, -read_time, -StartDate) %>%
  mutate(item = paste0("post_",item)) %>%
  spread(item, resp) %>% 
  rename(StartDate_post = StartDate)

df <- left_join(filter(df_pre, workerId %in% df_post$workerId), df_post, by="workerId") %>%
  mutate_at(vars(matches(".*(_[0-9]*r)")), funs(8 - .)) %>% # reverse code
  gather(item, resp, matches(".*(_[0-9])")) %>%
  mutate(phase = ifelse(grepl("pre_",item),"pre","post")) %>%
  mutate(
    item = gsub("pre_","",item),
    item = gsub("post_","",item)
  ) %>%
  mutate(scale = gsub("(_([0-9])+(r)*)","",item))

df_scored <- df %>%
  group_by(workerId, condition, scale, phase, parent) %>%
  summarize(mean_resp = mean(resp)) %>%
  spread(phase, mean_resp) %>%
  mutate(change = post - pre) %>%
  ungroup()

df_demo <- df %>% 
  select(-item, -resp, -phase, -scale) %>% 
  distinct(workerId, .keep_all = TRUE) %>% 
  mutate_at(vars(contains("StartDate")), lubridate::as_datetime)  %>% 
  mutate(days_between = StartDate_post - StartDate_pre)

s3_pre <- df_pre
s3_long <- df
s3_scored <- df_scored
s3_demo <- df_demo


s3_mdf <- s3_scored %>% 
  select(-change) %>% 
  gather(phase, mean, pre, post) %>% 
  mutate(condition = factor(condition)) %>% 
  # filter(scale == "vaccTox") %>%
  mutate(mean = rescale_beta(mean, 1, 7)) %>%
  spread(phase, mean) %>%
  mutate(condition = relevel(condition, ref="noInterv")) %>%
  mutate(evid = ifelse(condition=="noInterv",0,1))


s3_meta <- list(
  n_recruited = pre %>% filter(Progress==100) %>% distinct(workerId) %>% nrow(),
  n_eligible = df_pre %>% filter(return==26) %>% distinct(workerId) %>% nrow(),
  n_recruited_post = post %>% filter(Progress==100) %>%  distinct(workerId) %>%  nrow(),
  n_failed_post = post %>% filter(Progress==100, (check_3!=3|check_5!=5))  %>%  distinct(workerId) %>%  nrow(),
  n_complete = nrow(df_post),
  n_duplicate = 0
)

# clean up unused datasets
rm(pre, post, df_pre, df_post, df, df_demo)

# in the future refactor and clean up df_scored
