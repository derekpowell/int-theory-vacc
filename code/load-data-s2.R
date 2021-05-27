# load data for study 2
# Kara's code

library(tidyverse)

# make refactoring functions
as.numchar <- function(x){return(as.numeric(as.character(x)))}
as.likert7 <- function(x){
  x1 <- factor(tolower(as.character(x)),
               levels = c("1", "2", "3", "4", "5", "6", "7"),
               labels = c("strongly disagree", "disagree",
                          "somewhat disagree",
                          "neither agree nor disagree",
                          "somewhat agree",
                          "agree", "strongly agree"))
  x2 <- recode(x1,
               "strongly disagree" = -3,
               "disagree" = -2,
               "somewhat disagree" = -1,
               "neither agree nor disagree" = 0,
               "somewhat agree" = 1,
               "agree" = 2,
               "strongly agree" = 3)
  return(x2)
}

d_raw_pre <- read.csv("../data/vaccine+many+beliefs+-+interv+1+pre+-+APRIL_May+1%2C+2018_13.38.csv", encoding = "latin1") %>%
  rename(gender = Q49,
         age = Q51,
         ethnicity = Q53,
         education = Q55,
         job = Q57,
         income = Q59,
         political_party = Q61,
         political_beliefs = Q63_1)

# get question text
question_text_pre <- d_raw_pre[c(1,2),] %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("question") %>%
  rename("question_text" = "X1", question_importId = "X2") %>%
  mutate(question_importId = gsub('.*:\\"', '', question_importId),
         question_importId = gsub('\\"}', '', question_importId), 
         question_text = gsub(".* - ", "", question_text),
         question_text = gsub("Ûª", "'", question_text),
         question_text = gsub("â", "'", question_text))

# tidy dataset 
d_tidy_pre <- d_raw_pre[-c(1,2),] %>%
  filter(DistributionChannel == "anonymous") %>% 
  # mutate(StartDate = parse_datetime(StartDate, format = "%m/%d/%y %H:%M"),
  #        EndDate = parse_datetime(EndDate, format = "%m/%d/%y %H:%M"),
  #        RecordedDate = parse_datetime(RecordedDate, format = "%m/%d/%y %H:%M")) %>%
  mutate_at(vars(Progress, Duration..in.seconds., starts_with("Location"),
                 starts_with("Q"), age, study_time, payment),
            funs(as.numchar)) %>%
  mutate_at(vars(starts_with("hb"), starts_with("nat"), starts_with("medSkept"),
                 starts_with("disease"), starts_with("vacc"), 
                 starts_with("parent"),
                 starts_with("infant"), starts_with("overpar")),
            funs(as.likert7))

# reverse code
d_reverse_pre <- d_tidy_pre %>%
  mutate_at(vars(matches("_[0-9]r$")), # this is risky but it works
            funs(. * -1))


d_raw_post <- read.csv("../data/vaccine+many+beliefs+-+interv+1+post+g1+-+APRIL_May+1%2C+2018_13.40.csv", encoding = "latin1") %>%
  full_join(read.csv("../data/vaccine+many+beliefs+-+interv+1+post+g2+-+APRIL_May+1%2C+2018_13.40.csv", encoding = "latin1")) %>%
  full_join(read.csv("../data/vaccine+many+beliefs+-+interv+1+post+g3+-+APRIL_May+1%2C+2018_13.41.csv", encoding = "latin1")) %>%
  rename(is_parent = Q48,
         children_num = Q50,
         children_oldest = Q52,
         children_youngest = Q54,
         plan_parent = Q56,
         flushot_self_this = Q58,
         flushot_self_next = Q60,
         flushot_child_this = Q62,
         flushot_child_next = Q64,
         vax_delay_spread = Q66,
         vax_refuse = Q68,
         vax_exemption = Q70,
         attention = Q72,
         comments1 = Q74)

# get question text
question_text_post <- d_raw_post[c(1,2),] %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("question") %>%
  rename("question_text" = "X1", question_importId = "X2") %>%
  mutate(question_importId = gsub('.*:\\"', '', question_importId),
         question_importId = gsub('\\"}', '', question_importId), 
         question_text = gsub(".* - ", "", question_text),
         question_text = gsub("Ûª", "'", question_text),
         question_text = gsub("â", "'", question_text))

# tidy dataset 
d_tidy_post <- d_raw_post[-c(1,2),] %>%
  filter(DistributionChannel == "anonymous") %>% 
  # mutate(StartDate = parse_datetime(StartDate, format = "%m/%d/%y %H:%M"),
  #        EndDate = parse_datetime(EndDate, format = "%m/%d/%y %H:%M"),
  #        RecordedDate = parse_datetime(RecordedDate, format = "%m/%d/%y %H:%M")) %>%
  mutate_at(vars(Progress, Duration..in.seconds., starts_with("Location"),
                 starts_with("Q"), children_num, children_oldest, 
                 children_youngest, study_time, payment),
            funs(as.numchar)) %>%
  mutate_at(vars(starts_with("hb"), starts_with("nat"), starts_with("medSkept"),
                 starts_with("disease"), starts_with("vacc"), 
                 starts_with("parent"),
                 starts_with("infant"), starts_with("overpar")),
            funs(as.likert7)) %>%
  mutate_at(vars(ends_with("_parent"), starts_with("flushot_"),
                 starts_with("vax_"), attention),
            funs(. %>% as.character %>% as.numeric))

# reverse code
d_reverse_post <- d_tidy_post %>%
  mutate_at(vars(matches("_[0-9]r$")), # this is risky but it works
            funs(. * -1))

d_merge <- d_reverse_pre %>% 
  select(workerId, gender:political_beliefs, 
         comments, eligible, Duration..in.seconds., StartDate,
         starts_with("hb_"), starts_with("nat_"), starts_with("medSkept_"), 
         starts_with("disease"), starts_with("vacc"), starts_with("parent"), 
         starts_with("infant"), starts_with("overpar_")) %>%
  rename(duration = Duration..in.seconds.) %>%
  # mutate(phase = "pretest") %>%
  # gather(question, response, -c(workerId, phase)) %>%
  rename_at(vars(-c(workerId, gender:political_beliefs)), 
            funs(paste(., "pretest", sep = "_"))) %>%
  # distinct() %>%
  full_join(d_reverse_post %>%
              select(workerId, is_parent:comments, 
                     condition, Duration..in.seconds., StartDate,
                     starts_with("hb_"), starts_with("nat_"), 
                     starts_with("medSkept_"), starts_with("disease"),
                     starts_with("vacc"), starts_with("parent"),
                     starts_with("infant"), starts_with("overpar_")) %>%
              rename(duration = Duration..in.seconds.) %>%
              # mutate(phase = "post_test") %>%
              # gather(question, response, -c(workerId, phase)) %>%
              rename_at(vars(-workerId, -condition), 
                        funs(paste(., "posttest", sep = "_"))) %>%
              # distinct()) %>%
  filter(!is.na(workerId), !workerId %in% c("", " ", "NA")) %>%
  mutate(condition = factor(condition, levels = c("noInterv", "diseaseRisk")))
  )

d_tidy <- d_merge %>%
  mutate(duplicate = ifelse(workerId %in% data.frame(d_merge %>%
                                                       count(workerId) %>%
                                                       filter(n > 1))$workerId,
                            TRUE, FALSE),
         attn_fail_pretest = ifelse(eligible_pretest == "false",
                                    TRUE, FALSE),
         attn_fail_posttest = ifelse(hb_check_5_posttest != 5 - 4 |
                                       is.na(hb_check_5_posttest) |
                                       vaccEff_check_1_posttest != 1 - 4 |
                                       is.na(vaccEff_check_1_posttest) |
                                       vaccDanger_check_3_posttest != 3 - 4 |
                                       is.na(vaccDanger_check_3_posttest) |
                                       vaccTox_check_6_posttest != 6 - 4 |
                                       is.na(vaccTox_check_6_posttest),
                                     TRUE, FALSE),
         incomplete = ifelse(
           !workerId %in% data.frame(d_merge %>%
                                       drop_na(workerId, condition,
                                               starts_with("hb_"),
                                               starts_with("nat_"),
                                               starts_with("medSkept_"),
                                               starts_with("disease"),
                                               starts_with("vacc"),
                                               starts_with("parent"),
                                               starts_with("infant"),
                                               starts_with("overpar_")))$workerId,
           TRUE, FALSE)
         ) %>% 
  as_tibble()

# d_tidy %>% 
#   distinct(workerId, attn_fail_pretest, attn_fail_posttest, 
#            incomplete, duplicate) %>% 
#   count(attn_fail_pretest, attn_fail_posttest, incomplete, duplicate) 

d_all <- d_tidy %>%
  # filter(!duplicate) %>%
  # filter(attn_fail_pretest != TRUE, attn_fail_posttest != TRUE, 
  #        incomplete != TRUE, duplicate != TRUE) %>%
  # select(-c(attn_fail_pretest, attn_fail_posttest, incomplete, duplicate)) %>%
  # select(-incomplete, -duplicate) %>% 
  mutate(
    StartDate_pretest = lubridate::as_datetime(StartDate_pretest),
    StartDate_posttest = lubridate::as_datetime(StartDate_posttest)
  ) %>% 
  as_tibble()


d_keep <- d_all %>% 
  filter(!attn_fail_pretest, !attn_fail_posttest, !incomplete, !duplicate)
# write_csv(d, "data/study2_data.csv")

# d_raw <- read.csv("../data/vaccine+many+beliefs+-+interv+1+pre+-+APRIL_May+1%2C+2018_13.38.csv", encoding = "latin1")

# reformat to match previous analyses (i.e., 2 rows per participant)
d_demo <- d_keep %>% 
  select(workerId, condition, gender, age, ethnicity, education, job, income,
         political_party, political_beliefs, eligible_pretest, 
         is_parent_posttest, children_num_posttest, children_oldest_posttest, 
         children_youngest_posttest, plan_parent_posttest,
         starts_with("flushot_"), starts_with("vax_"), starts_with("attention_"),
         starts_with("comments"), starts_with("duration"), starts_with("attn_")
         )


d_pre <- d_keep %>% 
  select(workerId, ends_with("_pretest")) %>%
  select(-c(eligible_pretest, starts_with("flushot_"), starts_with("vax_"), 
            starts_with("attention_"), starts_with("comments"), 
            starts_with("duration"))) %>%
  rename_all(funs(gsub("_pretest", "", .))) %>%
  mutate(phase = "pre")

d_post <- d_keep %>% 
  select(-c(is_parent_posttest, children_num_posttest, children_oldest_posttest,
            children_youngest_posttest, plan_parent_posttest,
            starts_with("flushot_"), starts_with("vax_"), starts_with("attention_"),
            starts_with("comments"), starts_with("duration"))) %>%
  rename_all(funs(gsub("_posttest", "", .))) %>%
  mutate(phase = "post")

## problems here, redoing this below

# d <- bind_rows(d_pre, d_post) %>% 
#   gather(question, response, -c(workerId, phase)) %>%
#   mutate(phase = factor(phase,
#                         levels = c("pre", "post")),
#          reverse_cat = ifelse(grepl("_[1-9]r$", question), TRUE, FALSE),
#          # NOTE: "response" has already been reverse coded!
#          question = factor(question),
#          scale = factor(gsub("_.*$", "", question),
#                         levels = c("vaccIntent", "vaccDanger", "vaccEff", 
#                                    "vaccStrain", "vaccTox", 
#                                    "diseaseSevere", "diseaseRare", 
#                                    "infantImmLimCap", "infantImmWeak", 
#                                    "medSkept", "hb", "nat", 
#                                    "overpar", "parentExpert"))) %>%
#   full_join(d_demo) %>%
#   mutate(condition = factor(condition, levels = c("noInterv", "diseaseRisk"))) %>%
#   filter(!is.na(response), !is.na(workerId), !is.na(condition)) #%>%
#   # distinct()
# 
# # score all scales
# d_scored <- d %>%
#   # select(workerId, condition, phase, scale, response,
#   #        gender:duration_posttest) %>%
#   group_by(workerId, condition, phase, scale) %>%
#   mutate(response = as.numeric(response)) %>%
#   summarise(mean = mean(response, na.rm = TRUE)) %>%
#   ungroup()
#   # distinct()

# rename

s2_all <- d_all
s2_demo <- d_demo %>% as_tibble()

s2_long <- d_keep %>% 
  # filter(!attn_fail_pretest, !attn_fail_posttest, !incomplete, !duplicate) %>% 
  select(-attn_fail_pretest, -attn_fail_posttest, -incomplete, -duplicate, -attention_posttest) %>% 
  gather(item, response, matches("_[0-9]_"), matches("_[0-9]r_")) %>% 
  mutate(
    phase = ifelse(grepl("_pretest", item), "pre", "post"),
    scale = str_match(item, "([A-z]*)_")[,2],
    item = gsub("(_pretest|_posttest)", "", item)
    ) %>% 
  select(-comments_pretest, -comments1_posttest) %>% 
  rename(
    is_parent = is_parent_posttest,
    children_num = children_num_posttest,
    children_oldest = children_oldest_posttest,
    children_youngest = children_youngest_posttest,
    plan_parent = plan_parent_posttest,
    flushot_self = flushot_self_this_posttest,
    flushot_self_next = flushot_self_next_posttest,
    flushot_child = flushot_child_this_posttest,
    flushot_child_next = flushot_child_next_posttest,
    vax_delay_spread = vax_delay_spread_posttest,
    vax_refuse = vax_refuse_posttest,
    vax_exemption = vax_exemption_posttest
  ) %>% 
  relocate(workerId, condition, phase, scale, item, response)

s2_scored <- s2_long %>% 
  filter(!grepl("check",item)) %>% 
  group_by(workerId, condition, phase, scale) %>% 
  summarize(mean = mean(response, na.rm=TRUE)) %>% 
  ungroup()


s2_mdf <- s2_scored %>% 
  mutate(mean = rescale_beta(mean, -3, 3)) %>%
  spread(phase, mean) %>% 
  mutate(condition = relevel(condition, ref="noInterv")) %>%
  mutate(evid = ifelse(condition=="noInterv",0,1))

s2_meta <- list(
  n_recruited = d_all %>% distinct(workerId) %>% nrow(),
  n_eligible = d_all %>%  filter(eligible_pretest=="true") %>% distinct(workerId) %>%  nrow(),
  n_recruited_post = d_all %>% filter(!is.na(StartDate_posttest), !incomplete) %>% distinct(workerId) %>%  nrow(),
  n_failed_post = d_all %>% filter(!is.na(StartDate_posttest), !incomplete, !duplicate, attn_fail_posttest) %>% nrow(),
  n_complete = d_keep %>% filter(!is.na(StartDate_posttest)) %>% nrow(),
  n_duplicate = d_all %>% filter(duplicate, !attn_fail_posttest) %>% distinct(workerId) %>% nrow()
)

rm(d_all, d_keep, d_demo, d_pre, d_post, d_tidy, d_tidy_post, d_tidy_pre, d_raw_post, d_raw_pre )