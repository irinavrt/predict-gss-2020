library(tidyverse)


# Argument data -----------------------------------------------------------

issue_list <- read_csv("data/issue-list.csv")

# mf_data <- read_rds("data_processed/full_mf_data_with_predicted_year_slope.rds")
# mf_data_issp <- read_rds("data_processed/mf-data-or-issp-issues.rds")
# mf_data_1time <- read_rds("data_processed/mf-data-1time-asked.rds")
# 
# 
# mf_data_comb <- full_join(mf_data %>% 
#                             select(issue, hfvl = hfvl_advantage),
#                           mf_data_issp %>% 
#                             select(issue, hfvl_issp = hflv)) %>% 
#   bind_rows(mf_data_1time %>% select(issue, hfvl = hflv))
# 
# mf_data_comb <- mf_data_comb %>%
#   mutate(hfvl = ifelse(is.na(hfvl), hfvl_issp, hfvl)) %>%
#   select(-hfvl_issp)
# 
# mf_data_comb <- left_join(issue_list, mf_data_comb) %>% select(-asked_after_2010)
# 
# write_rds(mf_data_comb, "data/gss-issp-hfvl-measures.rds")


# Old code for 98 items -------------------------------------

amt <- read_csv("data/gss-arguments-raw-data.csv")

amt <- amt %>%
  select(code,
         polviews,
         intel,
         WorkerId,
         eval, time, time_took,
         quest = Input.question,
         answer = Answer.answer,
         arg = Answer.arg,
         counter_arg = Answer.counterarg,
         exposure = Answer.heard,
         belief_before = Answer.opinionbefore,
         belief_now = Answer.opinionnow) 


# Clean arguments data

full_arguments <- amt %>% 
  separate(arg, str_c("arg.", 1:9), fill = "right") %>% 
  separate(counter_arg, str_c("counter_arg.", 1:9), fill = "right") %>% 
  gather(order, mf, arg.1:counter_arg.9) %>% 
  drop_na(mf) %>% 
  separate(order, c("type", "order"), sep = "\\.") %>% 
  select(-order) %>% 
  mutate(mf = factor(mf, labels = c("harm", 
                                    "fair", 
                                    "ingr", 
                                    "auth", 
                                    "pure", 
                                    "lib", 
                                    "viol",
                                    "govern",
                                    "other")),
         value = 1) %>% 
  spread(mf, value) %>% 
  mutate_at(vars(harm:other), funs(replace_na(., 0))) %>% 
  gather(mf, value, harm:other) %>% 
  spread(type, value) %>% 
  mutate(pro = ifelse(answer == 1, arg, counter_arg),
         against = ifelse(answer == 0, arg, counter_arg)) %>% 
  select(-arg, -counter_arg) %>% 
  gather(type, value, pro, against) %>% 
  mutate(agree_position = ifelse(type == "pro", answer, 1 - answer))

full_arguments <- full_arguments %>% 
  mutate(position = str_c(issue, type, sep = "_"))

full_arguments <- full_arguments %>% 
  left_join(workers %>% 
              select(WorkerId,
                     group, aff, strength, polviews_cat,
                     round, age, sex, edu,
                     wordsum, wordsum_cat,
                     sci_knowledge, sci_kn_cat, 
                     oth_att, oth_att_cat))

full_arguments <- full_arguments %>% 
  left_join(time_trends %>% 
              select(variable, lib_side),
            by = c("issue" = "variable"))

write_rds(full_arguments, "data_processed/full_arg_positions_98v.rds")

full_arguments <- read_rds("data_processed/full_arg_positions_98v.rds")


mf_data <- full_arguments %>% 
  drop_na(value) %>% 
  group_by(issue, mf, type) %>% 
  summarise(value = mean(value)) %>% 
  spread(type, value) %>% 
  transmute(aa = pro - against)

mf_data <- mf_data %>% 
  spread(mf, aa)

mf_data <- mf_data %>% 
  mutate(hf = (harm + fair)/2, 
         hflv = (harm + fair + lib + viol)/4, 
         aip = (auth + ingr + pure)/3,
         aipg = (auth + ingr + pure + govern )/4)

write_rds(mf_data, "data_processed/mf-measures-98v-upd.rds")
mf_data <- read_rds("data_processed/mf-measures-98v-upd.rds")

harm_viol_upd <- full_arguments %>% 
  drop_na(value) %>% 
  filter(mf %in% c("harm", "viol")) %>% 
  group_by(WorkerId, issue, type) %>% 
  mutate(both_selected = sum(value) == 2)

harm_viol_upd <- harm_viol_upd %>% 
  ungroup() %>% 
  mutate(upd_value = ifelse(both_selected, value/2, value)) 

hv_data <- harm_viol_upd %>% 
  group_by(issue, mf, type) %>% 
  summarise(value = mean(upd_value)) %>% 
  spread(type, value) %>% 
  transmute(aa = pro - against)

hv_data <- hv_data %>% 
  spread(mf, aa) %>% 
  rename(harm_upd = harm, viol_upd = viol)

harm_or_viol <-  full_arguments %>% 
  drop_na(value) %>% 
  filter(mf %in% c("harm", "viol")) %>% 
  group_by(WorkerId, issue, type) %>% 
  summarise(sum_value = sum(value))

harm_or_viol <- harm_or_viol %>% 
  mutate(value = ifelse(sum_value > 0, 1, 0)) %>% 
  group_by(issue, type) %>% 
  summarise(value = mean(value)) %>% 
  spread(type, value) %>% 
  mutate(hv = pro - against) %>% 
  rename(hv_pro =  pro, hv_against = against)


write_rds(harm_or_viol, "data_processed/harm-viol-98v.rds")

# calculate mf measures ---------------------------------------------------

full_measures_with_pro_con <- amt %>% 
  separate(arg, str_c("arg.", 1:9)) %>% 
  separate(counter_arg, str_c("counter_arg.", 1:9)) %>% 
  gather(number, mf, arg.1:counter_arg.9) %>% 
  drop_na(mf) %>% 
  separate(number, c("type", "number"), sep = "\\.") %>% 
  select(-number) %>% 
  mutate(mf = factor(mf, labels = c("harm", 
                                    "fair", 
                                    "ingr", 
                                    "auth", 
                                    "pure", 
                                    "lib", 
                                    "viol",
                                    "govern",
                                    "other" )),
         value = 1) %>% 
  spread(mf, value) %>% 
  mutate_at(vars(harm:other), funs(ifelse(is.na(.), 0, 1))) %>% 
  gather(mf, value, harm:other) %>% 
  spread(type, value) %>% 
  mutate(pro = ifelse(answer == 1, arg, counter_arg),
         against = ifelse(answer == 0, arg, counter_arg)) %>% 
  select(-arg, -counter_arg) 

write_rds(full_measures_with_pro_con, "data_processed/full_arguments_new_data.rds")

mf_measures <- full_measures_with_pro_con %>% 
  filter(mf != "other") %>% 
  mutate(mf = str_c(mf, "_pro")) %>% 
  group_by(mf, issue) %>% 
  summarise(pro = mean(pro, na.rm = TRUE)) %>% 
  spread(mf, pro) %>% 
  left_join(full_measures_with_pro_con %>% 
              filter(mf != "other") %>% 
              mutate(mf = str_c(mf, "_against")) %>% 
              group_by(mf, issue) %>% 
              summarise(against = mean(against, na.rm = TRUE)) %>% 
              spread(mf, against)) %>% 
  left_join(full_measures_with_pro_con %>% 
              select(issue, WorkerId, exposure) %>% 
              distinct() %>% 
              group_by(issue) %>% 
              summarise(exposure = mean(exposure)))

mf_measures <- mf_measures %>%   
  mutate(hf_advantage = (harm_pro + fair_pro - harm_against - fair_against)/2, 
         hf_mean_argued = (harm_pro + fair_pro + harm_against + fair_against)/2,
         hfvl_advantage = (harm_pro + fair_pro + lib_pro + viol_pro - harm_against - fair_against - lib_against - viol_against)/4, 
         hfvl_mean_argued = (harm_pro + fair_pro + lib_pro + viol_pro - harm_against + fair_against + lib_against + viol_against)/4,
         lap_advantage = (auth_pro + ingr_pro + pure_pro - auth_against - ingr_against - pure_against)/3,
         lap_mean_argued = (auth_pro + ingr_pro + pure_pro + auth_against + ingr_against + pure_against)/3,
         lapg_advantage = (auth_pro + ingr_pro + pure_pro + govern_pro - auth_against - ingr_against - pure_against - govern_against)/4,
         lapg_mean_argued = (auth_pro + ingr_pro + pure_pro + govern_pro + auth_against + ingr_against + pure_against + govern_against)/4)
