# load data for study 2
# Kara's code

library(tidyverse)
d_all <- read.csv("study2_data.csv")[-1] %>%
  # filter out anyone who failed attention checks
  filter(hb_check_5_pretest == 1, hb_check_5_posttest == 1,
         vaccEff_check_1_pretest == -3, vaccEff_check_1_posttest == -3, 
         vaccDanger_check_3_pretest == -1, vaccDanger_check_3_posttest == -1, 
         vaccTox_check_6_pretest == 2, vaccTox_check_6_posttest ==2) %>%
  # remove attention checks from dataset
  select(-contains("check"))

# reformat to match previous analyses (i.e., 2 rows per participant)
d_demo <- d_all %>% 
  select(workerId, condition, gender, age, ethnicity, education, job, income,
         political_party, political_beliefs, eligible_pretest, 
         is_parent_posttest, children_num_posttest, children_oldest_posttest, 
         children_youngest_posttest, plan_parent_posttest,
         starts_with("flushot_"), starts_with("vax_"), starts_with("attention_"),
         starts_with("comments"), starts_with("duration"))

d_pre <- d_all %>% 
  select(workerId, ends_with("_pretest")) %>%
  select(-c(eligible_pretest, starts_with("flushot_"), starts_with("vax_"), 
            starts_with("attention_"), starts_with("comments"), 
            starts_with("duration"))) %>%
  rename_all(funs(gsub("_pretest", "", .))) %>%
  mutate(phase = "pre")

d_post <- d_all %>% 
  select(workerId, ends_with("_posttest")) %>%
  select(-c(is_parent_posttest, children_num_posttest, children_oldest_posttest,
            children_youngest_posttest, plan_parent_posttest,
            starts_with("flushot_"), starts_with("vax_"), starts_with("attention_"),
            starts_with("comments"), starts_with("duration"))) %>%
  rename_all(funs(gsub("_posttest", "", .))) %>%
  mutate(phase = "post")

d <- bind_rows(d_pre, d_post) %>% 
  gather(question, response, -c(workerId, phase)) %>%
  mutate(phase = factor(phase,
                        levels = c("pre", "post")),
         reverse_cat = ifelse(grepl("_[1-9]r$", question), TRUE, FALSE),
         # NOTE: "response" has already been reverse coded!
         question = factor(question),
         scale = factor(gsub("_.*$", "", question),
                        levels = c("vaccIntent", "vaccDanger", "vaccEff", 
                                   "vaccStrain", "vaccTox", 
                                   "diseaseSevere", "diseaseRare", 
                                   "infantImmLimCap", "infantImmWeak", 
                                   "medSkept", "hb", "nat", 
                                   "overpar", "parentExpert"))) %>%
  full_join(d_demo) %>%
  mutate(condition = factor(condition, levels = c("noInterv", "diseaseRisk"))) %>%
  filter(!is.na(response), !is.na(workerId), !is.na(condition)) %>%
  distinct()

# score all scales
d_scored <- d %>%
  select(workerId, condition, phase, scale, response,
         gender:duration_posttest) %>%
  group_by(workerId, condition, phase, scale) %>%
  mutate(response = as.numeric(response)) %>%
  summarise(mean = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct()
