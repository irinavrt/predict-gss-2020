
gss_aggr <- gss_aggr %>%
  rename(aa = hfvl)

# Transform public opinion variable.

gss_aggr <- gss_aggr %>%
  mutate(log_odds = log(mean_opin/(1 - mean_opin)))

# Select issues for training set
# that were measured in 2018 and at least twice before that.

issues_18 <- gss_aggr %>% 
  group_by(issue) %>% 
  filter(n() > 2, year == 2018) %>% 
  pull(issue)

# Create training and test sets.

gss_ts_expand <- gss_aggr %>%
  filter(issue %in% issues_18) %>% 
  expand_grid(last_year = seq(2010, 2016, by = 2))

train_ts <- gss_ts_expand %>% 
  filter(year <= last_year) 
test_ts <- gss_ts_expand %>% 
  filter(year == 2018) 

# Estimate drift parameter c_trend (by two year shift).

train_drift <- train_ts %>%
  group_by(last_year, issue, aa) %>% 
  summarise(last_mean_opin = last(mean_opin),
            last_log_odds = last(log_odds),
            c_trend = (last(log_odds) - first(log_odds))/((last(year) - first(year))/2))

# Estimate c_aa based on the relationship between c_trend and aa.

train_fit_pr <- train_drift %>%
  group_by(last_year) %>%
  mutate(c_aa = predict(lm(c_trend ~ aa, cur_data())))

# Make the forecast.

test_fc <- train_fit_pr %>% 
  inner_join(test_ts) %>% 
  group_by(last_year, issue) %>% 
  mutate(steps = (year - last_year)/2,
         fc_benchmark = last_log_odds,
         fc_trends = last_log_odds + c_trend*steps,
         fc_aa = last_log_odds + c_aa*steps) %>% 
  ungroup()

inv_logit <- inv_logit <- function(log_odds){
  exp(log_odds)/(1 + exp(log_odds))
}

test_fc %>% 
  mutate_at(vars(starts_with("fc_")), inv_logit) %>% 
  mutate(fc_resp_trend = inv_logit(fc_trend),
         resid_rw =  log_odds - last_log_odds,
         resid_trend =  log_odds - fc_trend, 
         resid_args = log_odds - fc_args) %>% 
  group_by(last_year) %>% 
  summarise(mean(abs(resid_trend)))

args_fc_acc <- args_fc %>% 
  filter(year2 == 2018) %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median) %>% 
  group_by(.model, last_year, issue, hfvl) %>% 
  summarise(MAE = mean(abs(error)),
            RMSE = sqrt(mean(error^2)))


accuracy_combined <- bind_rows(
  fc_acc %>% 
    mutate(.model = str_replace_all(.model, c("mod_rw" = "Benchmark", "mod_drift" = "Trends"))),
  args_fc_acc %>% 
    mutate(.model = "AA")
) %>% 
  as_tibble() %>% 
  mutate(.model = fct_relevel(.model, "Benchmark", "Trends"))

