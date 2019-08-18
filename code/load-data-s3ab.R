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
  distinct(workerId, .keep_all = TRUE) %>%
  filter(Progress==100, check_3==3, check_5==5) %>%
  mutate(parent = recode(parent, `23`=1, `24`=0)) %>%
  select(workerId, matches(".*(_[0-9])"), -check_3, -check_5, parent, study) %>%
  gather(item, resp, -workerId, -parent, -study) %>%
  mutate(item = paste0("pre_",item)) %>%
  spread(item, resp)

df_post <- post %>%
  filter(Progress==100, check_3==3, check_5==5) %>%
  select(workerId, matches(".*(_[0-9])"), -check_3, -check_5,
         condition, interv_t_Page.Submit) %>%
  rename(read_time = interv_t_Page.Submit) %>%
  gather(item, resp, -workerId, -condition, -read_time) %>%
  mutate(item = paste0("post_",item)) %>%
  spread(item, resp)


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

s3_pre <- df_pre
s3 <- df
s3_scored <- df_scored

# clean up unused datasets
rm(pre, post, df_pre, df_post, df)
