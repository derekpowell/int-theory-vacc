
# set wd (user-specific)
# setwd("/Users/kweisman/Documents/Research (Stanford)/Projects/Vaccines/vaccineBeliefs/study1/")
library(tidyverse)
library(stringr)

# working_dir <- getwd()
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

# load in raw data
d_raw <- read.csv("../data/study1_data.csv", header = T)

# filter rows and select columns
d1 <- d_raw %>%
  # filter(!grepl("Import", StartDate), !grepl("Start", StartDate)) %>% # remove extra headers
  filter(Status == "IP Address") %>% # remove survey previews (?)
  filter(Progress == 100) %>% # only get people who made it all the way through (?)
  filter(as.numeric(as.character(age)) >= 18) %>% # exclude people who said they were <18 yo (could be more stringent)
  select(StartDate, Duration..in.seconds., workerId, hb_1:comments) %>% # include only selected columns
  filter(!duplicated(workerId)) # we had a worker take it twice somehow! caught 5/1/19, 2:33 PM

# reformat variables
d2 <- d1 %>% 
  # mutate_at(vars(StartDate), funs(parse_datetime)) %>% # parse dates # 1/14/19, 3:12 PM this is throwing an error
  # mutate_at(vars(Duration..in.seconds., age, children), funs(as.character)) %>% # parse numbers step 1
  mutate_at(vars(Duration..in.seconds., age, children), funs(as.numeric)) %>% # parse numbers step 2
  mutate(educ = factor(educ, # recode ordinal variables
                       levels = c("No Diploma", "High School or Equivalent",
                                  "Some Undergraduate Education", "Undergraduate Degree",
                                  "Some Graduate or Professional Education",
                                  "Graduate or Professional Degree", "Doctorate")),
         income = factor(income, 
                         levels = c("Less than $20,000", "$20,000 - $30,0000",
                                    "$30,001 - $50,000", "$50,001 - $70,000",
                                    "$70,001 - $100,000", "More than $100,000")),
         parent = factor(parent, levels = c("No", "Yes")),
         expecting = factor(expecting, levels =c ("No", "Yes")),
         youngest_child = factor(youngest_child, 
                                 levels = c("Less than 1 year old", "1 year old", "2 years old", 
                                            "3 years old", "4 years old", "5 years old", "6 years old",
                                            "7 years old", "8 years old", "9 years old", "10 years old",
                                            "11 years old", "12 years old", "13 years old", "14 years old",
                                            "15 years old", "16 years old", "17 years old", "18 years old",
                                            "Over 18 years old"))) %>%
  gather(question, response, c(starts_with("hb_"),
                               starts_with("nat_"),
                               starts_with("medSkept_"),
                               starts_with("disease"),
                               starts_with("vacc"),
                               starts_with("parentE"),
                               starts_with("infant"),
                               starts_with("overpar_"))) %>% # convert to long form
  mutate(response_num = as.numeric(recode(response, # recode likert scales as numeric
                                          "Strongly disagree" = -3,
                                          "Disagree" = -2,
                                          "Somewhat disagree" = -1,
                                          "Neither agree nor disagree" = 0,
                                          "Somewhat agree" = 1,
                                          "Agree" = 2,
                                          "Strongly agree" = 3))) %>%
  mutate(response_num_rev = ifelse(str_sub(as.character(question), start = -1) == "r", # reverse-code as appropriate
                                   -1 * response_num, response_num),
         question_block = factor(gsub("_.*$", "", as.character(question))), # label question blocks
         attention_check = as.numeric(ifelse(grepl("check", as.character(question)), 
                                             str_sub(as.character(question), start = -1),
                                             NA))) %>% # label attention checks with correct answer
  rename(Duration = Duration..in.seconds.)

# find participants who failed any attention check
exclude_ids <- d2 %>%
  filter(!is.na(attention_check),
         response_num != attention_check -4) %>% # -4 to center
  distinct(workerId) %>%
  mutate(workerId = as.character(workerId))

# finish cleaning
d3 <- d2 %>%
  filter(!workerId %in% exclude_ids$workerId, # exclude participants
         grepl("check", question) == F) # omit check questions

# save dataframes for analysis

s1_all <- left_join(d1, exclude_ids %>% mutate(exclude=1), by="workerId") %>% 
  mutate(exclude = ifelse(is.na(exclude), 0 , 1)) %>% 
  as_tibble()

# demographics info (1 row per participant, demographics only)
d_demo <- d3 %>%
  select(StartDate:comments) %>%
  distinct()

# long-form (many rows per participant, all data)
d_long <- d3 %>% select(StartDate:comments, question_block, question:response_num_rev)

# wide-form (1 row per participant, all data)
d_wide <- d_long %>% 
  select(-question_block, -response, -response_num) %>%
  spread(question, response_num_rev)

# summary scores (1 row per participant, all data in score form)
d_sum <- d_long %>%
  group_by(workerId, question_block) %>%
  summarise(score = mean(response_num_rev, na.rm = T)) %>%
  ungroup() %>%
  full_join(d_demo)

s1_long <- d_long
s1_wide <- d_wide


# make dataframe for heatmpa, efa, etc.
# d_efa <- d_wide %>% select(-c(StartDate, Duration, sex:comments)) %>%
#   remove_rownames() %>%
#   column_to_rownames("workerId")

# make dataframe for graphical model fitting
d_bn <- d_sum %>%
  select(workerId, question_block, score) %>%
  spread(question_block, score) %>%
  remove_rownames() %>%
  column_to_rownames("workerId")

# remove extraneous dataframes
rm(d1, d2, d3, exclude_ids, d_long, d_sum)

# setwd(working_dir) # reset working directory