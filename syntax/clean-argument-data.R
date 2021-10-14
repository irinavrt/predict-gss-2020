library(tidyverse)

# Argument data -----------------------------------------------------------

issue_list <- read_csv("data/issue-list.csv") 

# Issues from a previous project

arg_data_previous <- read_csv("data/argument-data-previous-project-clean.csv") 

# New data
arg_data_new <- read_csv("data/arg_data_new_clean.csv") 

# Save participants' demographic information for sample description  
arg_data_new %>% 
  distinct(id, issue, age, sex, polviews, edu) %>% 
  write_rds("data/for-sample-description.rds")

# Combine the two data sources
arg_data_comb <- bind_rows(arg_data_previous %>% 
                            select(-polviews), 
                           arg_data_new %>% 
                             select(-polviews))

# Aggregate arguments per issue and calculate the advantage scores.

mf_data <- arg_data_comb %>% 
  filter(mf != "other") %>%
  rename(position = type) %>% 
  group_by(issue, position, mf) %>% 
  summarise(value = mean(value)) # value = the proportion of respondents who chose a given arg. as supporting a given issue-pos. 

mf_data <- mf_data %>%
  spread(position, value) %>% 
  transmute(mf, arg_adv = pro - against) %>% 
  spread(mf, arg_adv) %>% 
  ungroup()

mf_data <- mf_data %>% 
  mutate(aa = (Harm + Fairness + Liberty + Violence)/4) %>%
  select(issue, aa) %>% 
  # Add issues' wording, short labels, and version:
  left_join(issue_list %>% 
              select(issue, wording, short_label, version)) %>% 
  select(issue, wording, version, aa, short_label) 

write_rds(mf_data, "data/gss-issp-aa-measures.rds")
