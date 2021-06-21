library(tidyverse)
library(fable)
library(tsibble)
# library(feasts)
# library(tsibbledata)
# library(distributional)

gss_aggr <- read_rds("data/gss-issp-trends-to-predict.rds")
mf_data <- read_rds("data/gss-issp-hfvl-measures.rds")

gss_aggr <- left_join(gss_aggr, mf_data)

# select issues for training set
issues_18 <- gss_aggr %>% 
  group_by(issue) %>% 
  filter(n() > 2,
         year == 2018) %>% 
  pull(issue)

gss_aggr %>% 
  filter(issue %in% issues_18) %>% 
  group_by(issue) %>% 
  summarise(min = min(year), max = max(year), n = n()) %>% 
  arrange(n)

gss_ts <- gss_aggr %>%
  mutate(year2 = ifelse(year %% 2 == 1, year - 1, year)) %>% 
  group_by(issue, hfvl, year2) %>%
  summarise(mean_opin = mean(mean_opin), .groups = "drop") %>%
  as_tsibble(key = c(issue, hfvl), index = year2)

gss_ts <- fill_gaps(gss_ts)


logit <- function(p){
  log(p/(1 - p))
}

inv_logit <- function(log_odds){
  exp(log_odds)/(1 + exp(log_odds))
}

logit_transf <- new_transformation(logit, inv_logit)

gss_ts_fit <- gss_ts %>%
  model(mod = ARIMA(logit_transf(mean_opin) ~ 1 + pdq(0, 1, 0)),
        drift = RW(logit_transf(mean_opin) ~ drift())) 

gss_ts_fit$mod[[5]]$fit$model
gss_ts_fit$drift[[5]]$fit[1:4]

gt <- logit(gss_ts$mean_opin[gss_ts$issue == "abpoor"])
param_est(gt)

gss_ts_expand <- gss_ts %>%
  filter(issue %in% issues_18) %>% 
  expand_grid(last_year = seq(2010, 2016, by = 2)) %>%
  as_tsibble(key = c(last_year, issue, hfvl), index = year2)

train_ts <- gss_ts_expand %>% 
  filter(year2 <= last_year) 

test_ts <- gss_ts_expand %>% 
  filter(year2 > last_year) 

train_fit <- train_ts %>%
  model(mod_rw = ARIMA(logit_transf(mean_opin)  ~ 0 + pdq(0, 1, 0)),
        mod_drift = ARIMA(logit_transf(mean_opin)  ~ 1 + pdq(0, 1, 0))) 

# fix sigma for issue with t = 2 for 2012 sample, estimated at Inf
# train_fit$mod_drift[train_fit$issue == "spkmslm" & train_fit$last_year == 2012][[1]]$fit$model$sigma2["year"] <- 0
# train_fit$mod_drift[train_fit$issue == "spkmslm" & train_fit$last_year == 2012][[1]]$fit$fit$sigma2 <- 0
# train_fit$mod_drift[train_fit$issue == "colmslm" & train_fit$last_year == 2012][[1]]$fit$model$sigma2["year"] <- 0
# train_fit$mod_drift[train_fit$issue == "colmslm" & train_fit$last_year == 2012][[1]]$fit$fit$sigma2 <- 0


gss_fc <- train_fit %>%
  filter(map_lgl(mod_drift, ~class(.$fit) != "null_mdl")) %>% 
  forecast(test_ts, point_forecast = list(.median = median))

fc_acc <- gss_fc %>% 
  filter(year2 == 2018) %>% 
  mutate(y = inv_logit(log_odds), 
         error = y - .median) %>% 
  group_by(.model, last_year, issue, hfvl) %>% 
  summarise(MAE = mean(abs(error)),
            RMSE = sqrt(mean(error^2)))
  
# Argument predictions ----------------------------------------------------

# new_data <- gss_aggr_rscl %>% 
#   group_by(issue, hfvl) %>% 
#   summarise(within_10 = max(year) > 2008,
#             after_test = max(year) > 2010,
#             .groups = "drop") 

predicted <- train_fit %>%
  tidy() %>%
  group_nest(last_year) %>%
  mutate(mod = map(data, ~lm(estimate ~ hfvl, .x)),
         data = map2(data, mod, 
                     ~.x %>% 
                       mutate(hfvl_estimate = predict(.y, newdata = .x),
                              hfvl_se = predict(.y, newdata = .x, se.fit = TRUE)$se.fit,
                              hfvl_se = sqrt(hfvl_se^2 + summary(.y)$sigma^2),
                              diff = estimate - hfvl_estimate))) %>% 
  unnest(data)

# estimates_115 <- gss_ts %>%
#   as_tibble() %>%
#   group_by(issue, hfvl) %>% 
#   mutate(mean_opin_fill = zoo::na.approx(mean_opin), 
#          log_odds_fill = log(mean_opin_fill/(1 - mean_opin_fill)), 
#          last_10yrs = max(year2) >= 2008) %>% 
#   group_by(issue, hfvl, last_10yrs) %>% 
#   filter(year2 < 2012) %>% 
#   summarise(estimate = mean(log_odds_fill - lag(log_odds_fill), na.rm = TRUE), 
#             .groups = "drop") 

# hfvl_mod_115 <- lm(estimate ~ hfvl, estimates_115) 

# predicted <- predicted %>% 
#   mutate(hfvl_est_115 = predict(hfvl_mod_115, newdata = cur_data()),
#          hfvl_se_115 = predict(hfvl_mod_115, newdata = cur_data(), se.fit = TRUE)$se.fit,
#          hfvl_se_115 = sqrt(hfvl_se_115^2 + summary(hfvl_mod_115)$sigma^2),
#          diff_115 = estimate - hfvl_est_115) 

predicted_c <- left_join(tidy(train_fit),
                         predicted %>% 
                           select(issue, last_year,
                                  hfvl_estimate, hfvl_se, diff))
#                                  hfvl_est_115, hfvl_se_115, diff_115)) 
train_fit_augm <- augment(train_fit)
train_res <- left_join(train_fit_augm,
                       predicted_c  %>% 
                         select(-hfvl, -.model, -term, -statistic, -p.value))

args_spec <- train_res %>% 
  as_tibble() %>% 
  drop_na(estimate) %>%
  group_by(last_year, issue) %>% 
  mutate(hfvl_innov = .innov + diff) %>% 
#         hfvl_innov_115 = .innov + diff_115) %>% 
  summarise(hfvl_estimate = unique(hfvl_estimate),
            #hfvl_est_115 = unique(hfvl_est_115),
            estimate = unique(estimate),
            future = last(mean_opin),
            future_log_odds = logit_transf(future),
            mse_orig = mean(.innov^2, na.rm = TRUE),
            mse = mean(hfvl_innov^2, na.rm = TRUE),
            #mse_115 = mean(hfvl_innov_115^2, na.rm = TRUE),
            se_orig = unique(std.error),
            se_hfvl = unique(hfvl_se))
            #se_hfvl_115 = unique(hfvl_se_115))

test_param <- test_ts %>% 
  inner_join(args_spec) %>% 
  group_by(last_year, issue) %>% 
  mutate(steps = (year2 - last_year + 2)/2,
         fc_orig = future_log_odds + estimate*steps,
         fc = future_log_odds + hfvl_estimate*steps,
#         fc_115 = future_log_odds + hfvl_est_115*steps,
         se_orig = sqrt(mse_orig * steps + (steps * se_orig)^2),
         se = sqrt(mse * steps + (steps * se_hfvl)^2)) %>% 
#         se_115 = sqrt(mse_115 * steps + (steps * se_hfvl_115)^2)) %>% 
  ungroup()

# reproduce_res <- test_param %>% 
#   mutate(mean_opin = dist_normal(fc_orig, se_orig),
#          .model = "Trends")
# 
# reproduce_fc <- reproduce_res %>% 
#   as_fable(response = "mean_opin",
#            distribution = "mean_opin")
# 
# reproduce_fc <- reproduce_fc %>% 
#   mutate(mean_opin = dist_transformed(mean_opin, 
#                                        invert_transformation(structure(function (p) 
#                                        {
#                                          log(p/(1 - p))
#                                        }, class = "transformation", inverse = function (log_odds) 
#                                        {
#                                          exp(log_odds)/(1 + exp(log_odds))
#                                        })), 
#                                        logit_transf)) 
# 
# reproduce_acc <- reproduce_fc %>% 
#   filter(year2 == 2018) %>% 
#   accuracy(test_ts)
# 
# args_res_115 <- test_param %>% 
#   mutate(mean_opin = dist_normal(fc_115, se_115),
#          .model = "HFVL")
# 
# attr(args_res_115$mean_opin, "vars") <- "mean_opin"
# 
# args_fc_115 <- args_res_115 %>% 
#   as_fable(response = "mean_opin",
#            distribution = "mean_opin")
# 
# args_fc_115 <- args_fc_115 %>% 
#   mutate(mean_opin = dist_transformed(mean_opin, 
#                                       invert_transformation(structure(function (p) 
#                                       {
#                                         log(p/(1 - p))
#                                       }, class = "transformation", inverse = function (log_odds) 
#                                       {
#                                         exp(log_odds)/(1 + exp(log_odds))
#                                       })), 
#                                       logit_transf)) 
# 
# args_fc_acc_115 <- args_fc_115 %>% 
#   filter(year2 == 2018) %>% 
#   accuracy(test_ts)


args_res <- test_param %>% 
  mutate(mean_opin = dist_normal(fc, se),
         .model = "HFVL")

attr(args_res$mean_opin, "vars") <- "mean_opin"

args_fc <- args_res %>% 
  as_fable(response = "mean_opin",
           distribution = "mean_opin")

args_fc <- args_fc %>% 
  mutate(mean_opin = dist_transformed(mean_opin, 
                                      invert_transformation(structure(function (p) 
                                      {
                                        log(p/(1 - p))
                                      }, class = "transformation", inverse = function (log_odds) 
                                      {
                                        exp(log_odds)/(1 + exp(log_odds))
                                      })), 
                                      logit_transf)) 

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

accuracy_combined %>% 
  mutate(time_span = 2018 - last_year + 2, 
         time_span = paste(time_span, "years")) %>% 
  group_by(.model, time_span) %>% 
  summarise(MAE = round(mean(MAE)*100, 1)) %>% 
  spread(.model, MAE) %>% 
  mutate(ratio = AA*100/Benchmark, 
         ratio = paste0(round(ratio), "%")) %>% 
  mutate_at(vars(Benchmark:AA), ~sprintf("%.1f", .)) %>% 
  arrange(desc(time_span)) %>% 
  kableExtra::kbl(col.names = c("Time span of forecast", 
                                "Benchmark", 
                                "Trends", 
                                "AA", 
                                "Ratio AA: benchmark"))

accuracy_combined %>% 
  mutate(group = ifelse(abs(hfvl) > .16, "high adv", "low adv")) %>% 
  group_by(.model, last_year, group) %>% 
  summarise(RMSE = sqrt(mean(RMSE^2)), MAE = mean(MAE)) %>% 
  pivot_wider(names_from = group, values_from = c(RMSE, MAE)) %>% 
  mutate(last_year = last_year - 2) %>% 
  arrange(last_year) %>% 
  rename(Model = .model, `Last year` = last_year) %>% 
  kableExtra::kbl(digits = 3)

args_fc_115 %>% 
  filter(year2 == 2018) %>% 
  mutate(actual_value = inv_logit(log_odds),
         forecast = mean(mean_opin),
         error_resp = actual_value - forecast) %>% 
  as_tibble() %>% 
  drop_na(error_resp) %>% 
  group_by(last_year) %>% 
  summarise(ME = mean(error_resp),
            MAE = mean(abs(error_resp)),
            RMSE = sqrt(mean(error_resp^2)))


# sampling error ----------------------------------------------------------

ch <- gss_aggr_18 %>% 
  filter(year == 2018) %>% 
  semi_join(pl_data %>% select(issue)) %>% 
  group_by(issue) %>% 
  mutate(mean_opin_draws = list(map_df(1:1000, ~tibble(i = .x, p = mean(rbinom(n = n_sample, 1, prob = mean_opin))))))

ch %>% 
  select(issue, mean_opin, mean_opin_draws) %>% 
  unnest(mean_opin_draws) %>% 
  group_by(i) %>% 
  mutate(err = mean_opin - p) %>% 
  summarise(MAE = mean(abs(err))) %>% 
  summarise(MAE = mean(MAE))


# plot issue example ------------------------------------------------------

pl_spkhomo <- gss_ts %>% 
  filter(issue == "spkhomo") %>% 
  mutate(mean_opin_fill = zoo::na.approx(mean_opin), 
         log_odds_fill = log(mean_opin_fill/(1 - mean_opin_fill)),
         estimate = mean(log_odds_fill - lag(log_odds_fill), na.rm = TRUE),
         pred_log_odds = first(log_odds) + estimate*(0:(n()-1)),
         pred = inv_logit(pred_log_odds))

pl_spkhomo %>% 
  as_tibble() %>% 
  drop_na(log_odds) %>% 
  summarise(unique(estimate),
            sigma = sqrt(sum((log_odds - pred_log_odds)^2)/(n()-1)))

pl_spkhomo %>% 
  as_tibble() %>% 
  ggplot(aes(year2, mean_opin_fill))+
  geom_line() +
  geom_line(aes(y = pred), linetype = 2) +
  theme_classic() +
  scale_color_viridis_d() +
  scale_x_continuous(expand = expansion(mult = c(0.05, .07)),
                     breaks = seq(1970, 2020, by = 10)) +
  labs(x = NULL, y = "Public opinion") +
  theme(legend.position = c(.8, .2))


pl_colhomo <- pl_colhomo %>% 
  mutate(pred = ifelse(year2 > 2010, NA_real_, pred))

opin_2010 <- pl_colhomo$log_odds_fill[pl_colhomo$year2 == 2010]
hfvl_est <- 0.0880

pl_colhomo_adj <- pl_colhomo %>% 
  mutate(Benchmark = ifelse(year2 < 2010, NA_real_, opin_2010),
         Trends = Benchmark + ifelse(year2 <= 2010, 0, 1:4)*estimate,
         AA = Benchmark + ifelse(year2 <= 2010, 0, 1:4)*hfvl_est) %>% 
  mutate_at(vars(Benchmark:AA), inv_logit) %>% 
  select(Benchmark, Trends, AA) 


pl_colhomo_adj <- pl_colhomo_adj %>% 
  select(year2, Benchmark, Trends, AA) %>% 
  gather(Method, pred, -year2) %>% 
  drop_na() %>% 
  mutate(Method = fct_relevel(Method, "Benchmark"))

pl_colhomo %>% 
  as_tibble() %>% 
  ggplot(aes(year2, mean_opin_fill))+
  geom_line() +
  geom_line(aes(y = pred), color = "grey60") +
  geom_line(data = pl_colhomo_adj, aes(y = pred, color = Method)) +
  theme_classic() +
  scale_color_viridis_d() +
  scale_x_continuous(expand = expansion(mult = c(0.05, .07)),
                     breaks = seq(1970, 2020, by = 10)) +
  labs(x = NULL, y = "Public opinion") +
  theme(legend.position = c(.8, .2))


# plot forecasts ----------------------------------------------------------

pl_data <- args_fc %>% 
    filter(last_year == 2012, year2 %in% c(2018)) %>% 
    select(hfvl, last_value = future, mean_opin, log_odds) %>% 
    mutate(actual_value = inv_logit(log_odds),
           forecast = mean(mean_opin), 
           # lwr = quantile(mean_opin, 0.025), 
           # upr = quantile(mean_opin, 0.975),
           model = "HFVL") %>% 
    as_tibble() %>% 
    drop_na(forecast) 

  
pl_data <- pl_data %>% 
  gather(type, value, actual_value, forecast) %>% 
  mutate(issue = fct_reorder(issue, last_value), 
         type = ifelse(type == "actual_value", "True change", "Predicted change")) 

pl_data %>% 
  ggplot(aes(issue))+ 
  geom_segment(aes(xend = issue, y = last_value, yend = value, 
                   color = type), 
               data = filter(pl_data, type == "True change"),
               arrow = arrow(length = unit(0.1, "cm")), 
               lineend = "round",  linejoin = "mitre",
               size = .6, 
               position = position_nudge(y = 0, x = 0.2)) +
  geom_segment(aes(xend = issue, y = last_value, yend = value, 
                   color = type), 
               data = filter(pl_data, type == "Predicted change"),
               arrow = arrow(length = unit(0.1, "cm")), 
               lineend = "round",  linejoin = "mitre",
               size = .6, 
               position = position_nudge(y = 0, x = -0.2)) +
    coord_flip() +
  theme_bw() +
#  theme(legend.position = c(.8, .2)) +
  scale_color_viridis_d(end = 0.85) +
  labs(x = NULL, y = "Forecast", color = NULL) 

ggsave("forecast.jpeg", width = 7, height = 9, dpi = 350)

pl_data %>% 
  drop_na(forecast) %>% 
  mutate_at(vars(last_value, actual_value, forecast),
            ~ifelse(hfvl < 0, 1 - ., .)) %>% 
  mutate(issue = fct_reorder(issue, last_value, mean),
         model = fct_rev(model)) %>%
  ggplot(aes(issue))+ 
  geom_point(aes(y = forecast), color = "#CC6677", size = 1.5) +
  geom_segment(aes(xend = issue, y = last_value, yend = actual_value), 
               arrow = arrow(length = unit(0.1, "cm")), 
               lineend = "round",  linejoin = "mitre",
               size = 1, color = "grey30", alpha = .8) +
  coord_flip() + 
  facet_wrap(~year2, ncol = 4) +
  theme_bw() +
  labs(x = NULL, y = "Forecast") 

pl_diff_data <- pl_data %>% 
  drop_na(forecast) %>% 
  mutate_at(vars(last_value, actual_value, forecast),
            ~ifelse(hfvl < 0, 1 - ., .)) %>% 
  mutate(hfvl = abs(hfvl), 
         change = actual_value - last_value, 
         predicted_change = forecast - last_value) 

pl_diff_data %>% 
  group_by(year2) %>% 
  summarise(min_change = min(change),
            max_change = max(change),
            min_predicted_change = min(predicted_change),
            max_predicted_change = max(predicted_change))

panel_data <- tibble(year2 = c(2012, 2012, 2018, 2018),
                     change = c(-0.05, 0.04, -0.07, 0.222),
                     predicted_change = c(-0.05, 0.04, -0.07, 0.222))

pl_diff_data %>% 
  ggplot(aes(change, predicted_change))+
  geom_blank(data = panel_data) + 
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm") +
  facet_wrap(~year2, scales = "free") +
  labs(x = "Actual change", y = "Predicted change") +
  theme_bw()
  



# Forecast for 2020 -------------------------------------------------------

gss_res <- gss_ts_fit %>% 
  augment()

mse_median <- gss_res %>%
  as_tibble() %>%
  summarise(mse = median(.innov^2, na.rm = TRUE))  %>%
  pull()

estimates_full <- gss_ts %>%
  as_tibble() %>% 
  group_by(issue, hfvl) %>%
  mutate(mean_opin_fill = zoo::na.approx(mean_opin),
         log_odds_fill = log(mean_opin_fill/(1 - mean_opin_fill)),
         max_year = max(year2)) %>%
  group_by(issue, hfvl, max_year) %>%
  summarise(estimate = mean(log_odds_fill - lag(log_odds_fill), na.rm = TRUE),
            .groups = "drop")

hfvl_mod_full <- lm(estimate ~ hfvl, estimates_full)

fcst_2020 <- estimates_full %>% 
  filter(max_year >= 2008) %>% 
  mutate(hfvl_est = predict(hfvl_mod_full, newdata = .),
         hfvl_se = predict(hfvl_mod_full, newdata = ., se.fit = TRUE)$se.fit,
         hfvl_se = sqrt(hfvl_se^2 + summary(hfvl_mod_full)$sigma^2)) 

# Forecat sample

new_issues <- fcst_2020 %>% 
  filter(!issue %in% args_fc$issue[!is.na(args_fc$fc)],
         !issue %in% mf_data$issue[mf_data$version == "GSS"])

sample <- read_rds("data_processed/for-sample-description.rds")
sample <- sample %>% 
  mutate()

sample_sub <- sample %>% 
 filter(issue %in% new_issues$issue)

sample_sub %>%
  count(issue) %>% 
  summarise(mean(n), sum(n))

sample_sub %>% 
  mutate(age = as.numeric(age)) %>% 
  summarise(mean(age), sd(age), 
            mean(sex == "Female", na.rm = TRUE))

sample_sub %>% 
  mutate(age = as.numeric(age),
         polviews = as.numeric(polviews),
         polviews = cut(polviews, c(-1, 3, 6, 10))) %>% 
  count(polviews) %>% 
  mutate(n/sum(n))

sample_sub %>% 
  group_by(id, sex, age, polviews) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  count(n)

sample_sub %>% 
  group_by(id, sex, age, polviews) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  mutate(age = as.numeric(age)) %>% 
  summarise(n_distinct(id), mean(age), sd(age), 
            mean(sex == "Female", na.rm = TRUE))

fcst_2020_param <- left_join(fcst_2020,
                             gss_aggr %>% 
                               group_by(issue) %>% 
                               filter(year == last(year)) %>% 
                               select(issue, year, mean_opin) %>% 
                               expand_grid(year2 = seq(2020, 2030, by = 2))) %>% 
  group_by(issue) %>% 
  mutate(steps = (year2 - year)/2,
         future_log_odds = logit_transf(mean_opin),
         fc = future_log_odds + hfvl_est*steps,
         se = sqrt(mse_median * steps + (steps * hfvl_se)^2)) %>% 
  ungroup()


fcst_2020 <- fcst_2020_param %>% 
  mutate(mean_opin = dist_normal(fc, se))

attr(fcst_2020$mean_opin, "vars") <- "mean_opin"

fcst_2020 <- fcst_2020 %>% 
  as_fable(key = "issue",
           index = "year2",
           response = "mean_opin",
           distribution = "mean_opin")

fcst_2020 <- fcst_2020 %>% 
  mutate(mean_opin = dist_transformed(mean_opin, 
                                      invert_transformation(structure(function (p) 
                                      {
                                        log(p/(1 - p))
                                      }, class = "transformation", inverse = function (log_odds) 
                                      {
                                        exp(log_odds)/(1 + exp(log_odds))
                                      })), 
                                      logit_transf)) 

fcst_2020 <- fcst_2020 %>% 
  mutate(fc = median(mean_opin)) 

fcst_2020 %>% 
  as_tibble() %>% 
  left_join(mf_data %>% select(issue, wording)) %>% 
  mutate(wording = ifelse(issue %in% new_issues$issue, 
                          paste0(wording, "*"),
                          wording)) %>% 
  transmute(wording, issue, hfvl, year2, 
            last_opin = sprintf("%.2f (%.0f)",
                                inv_logit(future_log_odds),
                                year) ,
            fc = sprintf("%.2f", fc)) %>% 
  spread(year2, fc) %>% 
  arrange(desc(hfvl)) %>% 
  kableExtra::kbl(col.names = c("Item", "GSS code", "AA",
                                "Latest public opinion",
                                "Forecast 2020",
                                "Forecast 2022",
                                "Forecast 2024",
                                "Forecast 2026",
                                "Forecast 2028",
                                "Forecast 2030"),
                  digits = 2)


# Error in predictions for two different periods --------------------------


gss_ts_miss_yr <- gss_ts %>%
  filter(issue %in% sub_iss) %>% 
  expand_grid(last_year = c(1998, 2008, 2018)) %>%
  as_tsibble(key = c(last_year, issue, hfvl), index = year2)

train_ts <- gss_ts_miss_yr %>% 
  filter(year2 < last_year, year2 >= (last_year - 12)) 

test_ts <- gss_ts_miss_yr %>% 
  filter(year2 == last_year) 

train_fit <- train_ts %>%
  model(mod_rw = ARIMA(logit_transf(mean_opin)  ~ 0 + pdq(0, 1, 0)),
        mod_drift = ARIMA(logit_transf(mean_opin)  ~ 1 + pdq(0, 1, 0))) 

gss_fc <- train_fit %>%
  filter(map_lgl(mod_drift, ~class(.$fit) != "null_mdl")) %>% 
  forecast(test_ts, point_forecast = list(.median = median))

fc_acc <- gss_fc %>% 
  mutate(y = inv_logit(log_odds), 
         error = y - .median) %>% 
  group_by(.model, last_year, issue, hfvl) %>% 
  summarise(MAE = mean(abs(error)),
            RMSE = sqrt(mean(error^2)))

predicted <- train_fit %>%
  tidy() %>%
  group_nest(last_year) %>%
  mutate(mod = map(data, ~lm(estimate ~ hfvl, .x)),
         data = map2(data, mod, 
                     ~.x %>% 
                       mutate(hfvl_estimate = predict(.y, newdata = .x),
                              hfvl_se = predict(.y, newdata = .x, se.fit = TRUE)$se.fit,
                              hfvl_se = sqrt(hfvl_se^2 + summary(.y)$sigma^2),
                              diff = estimate - hfvl_estimate))) %>% 
  unnest(data)

predicted_c <- left_join(tidy(train_fit),
                         predicted %>% 
                           select(issue, last_year,
                                  hfvl_estimate, hfvl_se, diff))
train_fit_augm <- augment(train_fit)
train_res <- left_join(train_fit_augm,
                       predicted_c  %>% 
                         select(-hfvl, -.model, -term, -statistic, -p.value))

args_spec <- train_res %>% 
  as_tibble() %>% 
  drop_na(estimate) %>%
  group_by(last_year, issue) %>% 
  mutate(hfvl_innov = .innov + diff) %>% 
  summarise(hfvl_estimate = unique(hfvl_estimate),
            estimate = unique(estimate),
            future = last(mean_opin),
            future_log_odds = logit_transf(future),
            mse_orig = mean(.innov^2, na.rm = TRUE),
            mse = mean(hfvl_innov^2, na.rm = TRUE),
            se_orig = unique(std.error),
            se_hfvl = unique(hfvl_se))

test_param <- test_ts %>% 
  inner_join(args_spec) %>% 
  group_by(last_year, issue) %>% 
  mutate(steps = (year2 - last_year + 2)/2,
         fc_orig = future_log_odds + estimate*steps,
         fc = future_log_odds + hfvl_estimate*steps,
         se_orig = sqrt(mse_orig * steps + (steps * se_orig)^2),
         se = sqrt(mse * steps + (steps * se_hfvl)^2)) %>% 
  ungroup()

args_res <- test_param %>% 
  mutate(mean_opin = dist_normal(fc, se),
         .model = "HFVL")

attr(args_res$mean_opin, "vars") <- "mean_opin"

args_fc <- args_res %>% 
  as_fable(response = "mean_opin",
           distribution = "mean_opin")

args_fc <- args_fc %>% 
  mutate(mean_opin = dist_transformed(mean_opin, 
                                      invert_transformation(structure(function (p) 
                                      {
                                        log(p/(1 - p))
                                      }, class = "transformation", inverse = function (log_odds) 
                                      {
                                        exp(log_odds)/(1 + exp(log_odds))
                                      })), 
                                      logit_transf)) 

args_fc_acc <- args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median) %>% 
  group_by(.model, last_year, issue, hfvl) %>% 
  summarise(MAE = mean(abs(error)),
            RMSE = sqrt(mean(error^2)))

args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  select(issue, hfvl, year2, error) %>% 
  spread(year2, error) %>%
  drop_na() %>% 
  ggplot(aes(`1998`, `2018`, label = issue, color = hfvl)) +
  geom_point() +
  geom_text() +
  scale_color_viridis_c() +
  theme_classic()

ggsave("plot.jpeg")



args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  group_by(issue) %>% 
  select(issue, hfvl, year2, error) %>% 
  spread(year2, error) %>%
  filter(!str_detect(issue, "^ab")) %>% 
  select(-issue) %>% 
  cor(use = "pair") %>% 
  round(2)

args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median, 
         dist_half = abs(future - .5)) %>% 
  as_tibble() %>% 
  group_by(last_year) %>% 
  summarise(cor(dist_half, abs(error), use = "pair"))

args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median, 
         dist_half = abs(future - .5)) %>% 
  as_tibble() %>% 
  ggplot(aes(dist_half, abs(error))) +
  geom_point() +
  facet_wrap(~last_year) + 
  geom_smooth(method = "lm")

args_fc %>% 
    mutate(.median = median(mean_opin),
           y = inv_logit(log_odds), 
           error = y - .median) %>% 
    as_tibble() %>% 
    ungroup() %>% 
    select(issue, hfvl, year2, error) %>% 
    spread(year2, error) %>%
    drop_na() %>% 
  summarise(cor(`2008`, `2018`))


args_fc %>% 
  mutate(.median = median(mean_opin),
         y = inv_logit(log_odds), 
         error = y - .median) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  select(issue, hfvl, year2, error) %>% 
  spread(year2, error) %>%
  drop_na() %>% 
  mutate(mean = (`2008` + `2018`)/2) %>% 
  summarise(cor(mean, hfvl))


accuracy_combined <- bind_rows(
  fc_acc %>% 
    mutate(.model = str_replace_all(.model, c("mod_rw" = "Benchmark", "mod_drift" = "Trends"))),
  args_fc_acc %>% 
    mutate(.model = "AA")
) %>% 
  as_tibble() %>% 
  mutate(.model = fct_relevel(.model, "Benchmark", "Trends"))

accuracy_combined %>% 
  group_by(.model, last_year) %>% 
  summarise(MAE = round(mean(MAE)*100, 1)) %>% 
  spread(.model, MAE) %>% 
  mutate(ratio = AA*100/Benchmark, 
         ratio = paste0(round(ratio), "%")) %>% 
  mutate_at(vars(Benchmark:AA), ~sprintf("%.1f", .)) %>% 
  arrange(desc(time_span)) %>% 
  kableExtra::kbl(col.names = c("Time span of forecast", 
                                "Benchmark", 
                                "Trends", 
                                "AA", 
                                "Ratio AA: benchmark"))



train_ts %>% 
  filter(issue == "colhomo", last_year == 2018) %>% 
  model(mod = RW(logit_transf(mean_opin)  ~ drift())) %>% 
  forecast(test_ts)


object <- train_ts %>% 
  filter(issue == "colhomo", last_year == 2018) %>% 
  model(mod = RW(logit_transf(mean_opin)  ~ drift())) 
bind_new_data(object, test_ts)

ch_fc <- train_ts %>% 
  filter(issue == "colhomo", last_year == 2018) %>% 
  model(mod = ARIMA(logit_transf(mean_opin)  ~ 1 + pdq(0, 1, 0)),
        mod_lin = ARIMA(mean_opin  ~ 1 + pdq(0, 1, 0))) %>% 
  forecast(h = 2, point_forecast = lst(mean, median)) %>% 
  hilo()

ch_fc$mean_opin[[3]] %>% class()
ch_fc$mean_opin[[3]] %>% str()

new_data <- train_fit %>% 
  filter(issue == "colhomo", last_year == 2018) %>% 
  fabletools:::bind_new_data(train_ts) %>% 
  pull(new_data) %>% 
  .[[1]]

b.se <- 0.0286
fc <- rep(object$future, fullperiods)[1:h] + steps * b
mse <- mean(residuals(ch)$.resid^2, na.rm = TRUE)
se <- sqrt(mse * steps + (steps * b.se)^2)
distributional::dist_normal(fc, se)


ch <- gss_ts_fit$mod[gss_ts_fit$issue == "marhomo"][[1]] 


ch_nomiss <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0))) %>%
  pull(mod) %>% 
  .[[1]]


gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0)), 
        drift = RW(log(odds) ~ drift())) %>% 
  augment() %>% 
  as_tibble() %>% 
  pivot_wider(names_from = .model, values_from = c(.fitted:.innov))

gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0)), 
        drift = RW(log(odds) ~ drift())) %>% 
  forecast(h = 5, point_forecast = lst(mean, median)) %>% 
  hilo()


gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0)), 
        drift = RW(log(odds) ~ drift())) %>% 
  augment() %>% 
  as_tibble() %>% 
  group_by(.model) %>% 
  slice(-1) %>% 
  summarise(n = n(), 
            se = sqrt(sum(.innov^2, na.rm = TRUE)/(n-2)),
            se2 = sqrt(sum(.innov^2, na.rm = TRUE)/(n-1)))

gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log(odds)  ~ 1 + pdq(0, 1, 0)), 
        drift = RW(log(odds) ~ drift())) %>% 
  tidy()


homo_ts <- gss_aggr %>% 
  as_tsibble(key = c(issue, hfvl), index = year) %>% 
  filter(issue == "marhomo") %>% 
  fill_gaps() %>% 
  pull(log_odds) %>% 
  ts()
homo_ts2 <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  pull(log_odds) %>% 
  ts()

xreg <- matrix(seq(1:length(homo_ts)), dimnames = list(NULL, "constant"))
ch <- arima(homo_ts, order = c(0, 1, 0),
             xreg = xreg)

xreg2 <- matrix(seq(1:length(homo_ts2)), dimnames = list(NULL, "constant"))
ch2 <- arima(homo_ts2, order = c(0, 1, 0),
            xreg = xreg2)

ch$coef 



ch %>%
  residuals()

ch2

ci_boot <- gss_ts_fit %>% 
  filter(issue == "colhomo") %>%
  forecast(h = 1, bootstrap = TRUE, times = 1000)
  
ci_boot %>% 
  summarise(mean = mean(log_odds), median = median(log_odds), lower = quantile(log_odds, 0.025), higher = quantile(log_odds, 0.975)) %>% 
  mutate(lower - mean, higher - mean)

1.96*sqrt(0.01725)
ci_boot %>% hilo()
  
drift_test <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = RW_adj(log(odds) ~ drift()))
ar_test <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log(odds) ~ 1 + pdq(0, 1, 0)))

gss_ts %>% 
  filter(issue == "marhomo") %>% 
  left_join(augment(drift_test)) %>% 
  left_join(tidy(drift_test) %>% select(issue, estimate)) %>% 
  mutate(hfvl_estimate = 0.000855 + 0.267*hfvl, 
         diff = estimate - hfvl_estimate,
         hfvl_fitted = exp(log(.fitted) - diff), 
         hfvl_resid = exp(.innov + diff)) %>% 
  select(-estimate:-diff)
  
  
  
gss_ts %>% 
  filter(issue == "marhomo")


drift_test %>% 
  tidy()
drift_test %>% 
  model_sum()

drift_test %>% 
  forecast() %>% 
  hilo()

ch <- gss_ts_fit %>% 
  filter(issue == "marhomo") %>% 
  pull(mod) %>% 
  .[[1]]

ch$fit$fit$sigma2 %>% sqrt()

sqrt(sum(ch$residuals^2, na.rm = TRUE)/(length(homo_ts) - 2))



mod_test <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log_odds  ~ 0 + pdq(0, 1, 0)))

ch <- mod_test$mod[[1]]
class(ch)

ch$fit$model$residuals

ch$fit %>% str()

ch_fc <- forecast(ch, h = 1, bootstrap = TRUE, times = 100)
ch_fc2 <- forecast(ch, h = 1, bootstrap = TRUE, times = 100)

hilo(ch_fc)
table(ch_fc$log_odds[[1]]$x - ch_fc$.mean[[1]] )
table(ch_fc2$log_odds[[1]]$x - ch_fc2$.mean[[1]] )

hilo(ch_fc2)

ch_fc %>% 
  hilo()

sqrt(sum(ch$fit$model$residuals^2, na.rm = TRUE)/(length(homo_ts) - 2))

ch_adj <- ch
new_const <- 0.1119851
ch_adj$fit$model$coef["constant"] <- new_const
ch_adj$fit$model$residuals <- ch_adj$fit$model$residuals + (ch$fit$model$coef - new_const)

sqrt(sum(ch_adj$fit$model$residuals^2, na.rm = TRUE)/(length(homo_ts) - 2))

x <- ch_adj$fit
h <- 1
new_data <- newdata

intercept <-  fable:::arima_constant(NROW(x$est) + h, x$spec$d, 
                            x$spec$D, x$spec$period)[NROW(x$est) + seq_len(h)]
xreg <- matrix(intercept, dimnames = list(NULL, "constant"))

narma <- sum(x$model$arma[1L:4L])
coef <- x$model$coef
coef <- coef[(narma + 1):length(coef)]

off <- 
xm <- drop(xreg %*% coef)

res <- residuals(x)
new_data$.innov <- sample(na.omit(res) - mean(res, 
                                              na.rm = TRUE), nrow(new_data), replace = TRUE)
stats:::predict.Arima()
stats:::arima.sim()

newdata <- fabletools:::make_future_data(ch$data, 1)
ch_adj_fc <- forecast(ch_adj$fit, new_data = newdata, h = 1, bootstrap = TRUE, times = 100) %>% 
  hilo() 
ch_adj$fit %>% str()


ch_adj_fc$log_odds[[1]]$x %>% mean()


fable:::generate.ARIMA(mod_test$mod[[1]]$fit, times = 100)

fable:::forecast.ARIMA(mod_test$mod[[1]]$fit)

resid_const <- stats::lag(homo_ts) - homo_ts - ch2$coef
resid_args <- stats::lag(homo_ts) - homo_ts - 0.1119851

ch2_arg <- arima(resid_args, order = c(0, 1, 0))
ch2_arg
ch2_resid <- arima(resid_const, order = c(0, 1, 0))
ch2_resid

resid(ch2_arg)
resid(ch2_resid)

ch3 <- arima(homo_ts, order = c(0, 1, 0))

ch2$residuals
ch2_arg$residuals
ch2_resid$residuals
ch3$residuals

xreg
ch2_no_miss <- arima(homo_ts_no_miss, order = c(0, 1, 0),
             xreg = xreg)

ch2_no_miss %>% str()
stats:::predict.Arima()

predicted_c %>% 
  filter(issue == "marhomo") %>% 
  select(estimate, hfvl_estimate)

newdata <- matrix(seq(1,length(homo_ts)+5), dimnames = list(NULL, "constant"))

fitted <- predict(ch2, n.ahead = 5, newxreg = newdata)$pred
stats::predict.Arima()

predict(ch3)

arima(homo_ts, order = c(0, 1, 0), xreg = time(homo_ts))



nstar <- length(homo_ts_no_miss) - 1
nstar <- length(homo_ts[!is.na(homo_ts)]) - 1
npar <- length(ch2_no_miss$coef[ch2_no_miss$mask]) + 1
new_sigma2 <- sum(ch2_no_miss$residuals^2, na.rm = TRUE) / (nstar - npar + 1)
ch2$sigma2
mean(ch2$residuals^2, na.rm = TRUE) 
sum(ch2$residuals^2, na.rm = TRUE) / (nstar - npar + 1)

ch2$var.coef %>% sqrt()
ch2 %>% tidy()
ch2$sigma2
ch$fit$model$sigma2

new_data <- fabletools:::make_future_data(ch$data, h = 5)
forecast(ch$fit, new_data)

# 0.03584:  log likelihood = 5.05,  aic = -8.11

gss_ts %>% 
  filter(issue == "marhomo") %>% 
   model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0)),
        drift = RW(log_odds ~ drift())) %>% 
  tidy()

gss_ts %>% 
  filter(issue == "homosex") %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0)),
        drift = RW(log_odds ~ drift())) %>% 
  tidy()

zoo::na.approx()

gss_ts %>% 
  filter(issue == "homosex") %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0))) %>% 
  interpolate(gss_ts) %>% 
  as_tibble() %>% 
  mutate(diff = log_odds - lag(log_odds)) %>% 
  drop_na() %>% 
  summarise(mean(diff))

  # augment() %>% 
  # pivot_wider(names_from = .model, values_from = c(.fitted, .resid, .innov)) %>% View


gss_ts %>% 
  filter(issue == "marhomo") %>% 
  model(mod = ARIMA(log_odds  ~ pdq(0, 1, 0)),
        drift = RW(log_odds ~ drift())) %>% 
  tidy()

ch_diff <- gss_ts %>% 
  filter(issue == "marhomo") %>% 
  as_tibble() %>% 
  drop_na() %>% 
  mutate(diff = log_odds - lag(log_odds),
         lag = (year2 - lag(year2))/2, 
         diff_adj = diff/lag) %>% 
  drop_na(diff_adj) %>% 
  pull(diff_adj)

lm(ch_diff ~ 1) %>% tidy()

gss_ts_fit %>% 
  filter(issue == "homosex") %>% 
  tidy()



train_fit %>% 
  filter(issue == "homosex", last_year == 2018) %>% 
  gg_tsresiduals()

train_fit %>% 
  filter(issue == "homosex", last_year == 2018) %>% 
  forecast(h = 1) %>% 
  autoplot(gss_ts_miss_yr)

pred_int <- trends_fc_boot %>% 
  mutate(median = median(log_odds),
         lower = quantile(log_odds, 0.025),
         higher = quantile(log_odds, 0.925)) %>% 
  select(issue, log_odds:higher) %>% 
  left_join(gss_ts %>% 
              filter(year2 == 2018) %>% 
              select(issue, actual = log_odds)) 

pred_int %>% 
  mutate_at(vars(median:higher, actual), 
            ~ exp(.)/(exp(.) + 1)) %>% 
  drop_na(actual) %>% 
  ggplot(aes(fct_reorder(issue, actual), median, 
             ymin = lower, ymax = higher))+
  geom_pointrange() +
  geom_point(aes(y = actual), color = "red", size = 3)+
  coord_flip()

trends_fc <- train_fit %>%
  filter(last_year == 2018) %>% 
  forecast(h = 1)

trends_fc %>% 
  filter(issue == "abhlth") %>% 
  hilo(level = c(95))

trends_fc_boot %>% 
  filter(issue == "abhlth") %>% 
  hilo(level = c(95))

trends_fc_sim %>% 
  filter(issue == "abhlth") %>% 
  hilo(level = c(95))

train_fit %>%
  filter(last_year == 2018, issue == "marhomo") %>% 
  glance()

trends_fc %>% 
  filter(issue == "abhlth") %>% 
  autoplot(gss_ts_miss_yr)

accuracy(gss_fc, gss_ts_miss_yr) %>% 
  group_by(last_year, .model) %>% 
  summarise_at(vars(MAE, MAPE:ACF1), mean, na.rm = TRUE)


ch <- accuracy(gss_fc, gss_ts_miss_yr, list(qs=quantile_score), probs=0.10) 

gss_fc %>% 
  filter(issue == "abhlth") %>% 
  autoplot(gss_ts_miss_yr)


# Check n years threshold -------------------------------------------------

yrs <- 3

thr_year <- gss_aggr %>% 
  group_by(issue) %>% 
  filter(year < max(year) - yrs) %>% 
  summarise(thr_year = max(year))


gss_aggr <- left_join(gss_aggr, thr_year)

by_issue <- gss_aggr %>% 
  group_by(issue) %>% 
  filter(n() > 5) %>% 
  nest()

by_issue <- by_issue %>% 
  mutate(data_sub = map(data, ~filter(., year <= thr_year)),
         data_last = map(data, ~filter(., year %in% c(thr_year, max(year)))))

fit_lin <- function(data){
  lm(mean_opin ~ year_r, 
     data, weights = n_sample)
}

by_issue <- by_issue %>% 
  mutate(mod_lin = map(data_sub, fit_lin))

by_issue <- by_issue %>% 
  mutate(mod_lin_last = map(data_last, fit_lin))

by_issue <- by_issue %>% 
  mutate(map_df(mod_lin, coef)) %>% 
  rename(int_lin = `(Intercept)`, slope_lin = year_r) %>% 
  mutate(map_df(mod_lin_last, coef)) %>% 
  rename(int_lin_last = `(Intercept)`, slope_lin_last = year_r)

by_issue <- left_join(by_issue, mf_data %>% select(issue = code, hfvl = hfvl_advantage))

by_issue <- by_issue %>% drop_na(hfvl)
by_issue %>% 
  mutate(map_df(data_last, ~summarise(., act_year_diff = max(year) - min(year)))) %>% 
  select(issue, act_year_diff) %>% 
  ungroup() %>% 
  count(act_year_diff)


m_tr <- lm(slope_lin_last ~ slope_lin, by_issue)
m_slope <- lm(slope_lin_last ~ hfvl, by_issue)
m_both <- lm(slope_lin_last ~ slope_lin +  hfvl, by_issue)

modelsummary::msummary(list(m_tr, m_slope, m_both), 
                       statistic = "conf.int")



# subset for 2020 ---------------------------------------------------------

issues2018 <- gss_bin %>% 
  drop_na(opinion) %>% 
  group_by(issue) %>% 
  filter(2018 %in% year, 2016 %in% year) %>% 
  distinct(issue) %>% 
  pull(issue)

gss_aggr <- gss_bin %>% 
  drop_na(opinion) %>% 
  filter(issue %in% issues2018) %>% 
  group_by(issue, year) %>%
  summarise(mean_opin = weighted.mean(opinion, wgt),
            n_agree = sum(opinion*wgt),
            n_sample = sum(wgt),
            n_disagree = n_sample - n_agree,
            .groups = "drop")

gss_aggr <- left_join(gss_aggr, mf_data %>% select(issue = code, hfvl = hfvl_advantage))

gss_aggr <- expand_grid(gss_aggr, last_year = seq(2010, 2018, by = 2))

gss_aggr <- gss_aggr %>% 
  mutate(year_r = (year - (last_year - 2))/2)

by_issue_year <- gss_aggr %>% 
  group_nest(issue, hfvl, last_year) %>% 
  mutate(train_data = map2(data, last_year, ~filter(.x, year < .y)),
         test_data = map2(data, last_year, ~filter(.x, year >= .y - 2)), 
         test_2yrs = map2(data, last_year, ~filter(.x, year %in% c(.y - 2, .y))))
  
fit_mod <- function(data){
  glm(cbind(n_agree, n_disagree) ~ year_r,
      data, family = quasibinomial(link = "logit"))
}

fit_lin <- function(data){
  lm(mean_opin ~ year_r, 
     data, weights = n_sample)
}

coef_by_issue <- by_issue_year %>% 
  mutate(mod_log = map(train_data, fit_mod),
         mod_lin = map(train_data, fit_lin),
         mod_lin_test = map(test_data, fit_lin),
         mod_lin_test_2yr = map(test_2yrs, fit_lin)) %>% 
  mutate(map_df(mod_log, coef)) %>% 
  rename(int_log = `(Intercept)`, slope_log = year_r) %>% 
  mutate(map_df(mod_lin, coef)) %>% 
  rename(int_lin = `(Intercept)`, slope_lin = year_r) %>% 
  mutate(map_df(mod_lin_test, coef)) %>% 
  rename(int_lin_test = `(Intercept)`, slope_lin_test = year_r) %>% 
  mutate(map_df(mod_lin_test_2yr, coef)) %>% 
  rename(int_lin_2yrs = `(Intercept)`, slope_lin_2yrs = year_r) 

coef_by_issue %>% 
  group_by(last_year) %>% 
  summarise(cor(hfvl, slope_lin, use = "pair"), 
            cor(hfvl, slope_lin_test),
            cor(hfvl, slope_lin_2yrs), 
            cor(slope_lin, slope_lin_test, use = "pair"),
            cor(slope_lin, slope_lin_2yrs, use = "pair")) 

coef_by_issue %>% 
  ggplot(aes(hfvl, slope_lin))+
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~last_year)

coef_by_issue %>% 
  ggplot(aes(hfvl, slope_lin))+
  geom_point() +
  geom_smooth() +
  facet_wrap(~last_year)

coef_by_issue %>% 
  select(last_year, hfvl, slope_lin_2yrs) %>% 
  spread(last_year, slope_lin_2yrs) %>% 
  cor()

coef_by_issue %>% 
  filter(last_year == 2012) %>% 
  ggplot(aes(hfvl, slope_lin_2yrs))+
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text(aes(label = issue), size = 3)



# Multilevel version ------------------------------------------------------

gss_aggr %>% 
  group_by(issue) %>% 
  filter(n_distinct(year) == 2, 
         last_year == 2018)

gss_aggr <- gss_aggr %>% 
  group_by(issue) %>% 
  filter(n_distinct(year) > 2) %>% 
  ungroup()

gss_aggr <- gss_aggr %>% 
  mutate(year = ifelse(year %in% seq(2015, 1995, by = -2), year + 1, year)) 

by_year <- gss_aggr %>% 
  group_nest(last_year)


library(lme4)

fit_mlm <- function(data){
  glmer(cbind(n_agree, n_disagree) ~ year_r*hfvl + (year_r|issue), 
        family = binomial(),
        data)
}

predict_future <- function(mod, data){
  data %>% 
    mutate(pred = predict(mod, newdata = data, type = "response", re.form = ~ (1|issue)))
}

by_year <- by_year %>% 
  mutate(data_train = map2(data, last_year, ~filter(.x, year < .y)),
         mod = map(data_train, fit_mlm)) 
         
by_year <- by_year %>% 
  mutate(pred_data = map2(mod, data, predict_future))


by_year_prediction <- by_year %>% 
  select(last_year, pred_data) %>% 
  unnest(pred_data)

by_year_prediction %>% 
  filter(issue == "marhomo") %>% 
  ggplot(aes(year, pred)) +
  geom_line(aes(color = factor(last_year))) +
  geom_line(data = filter(gss_aggr, issue == "marhomo"), aes(y = mean_opin))


by_year_prediction %>%
  filter(year >= last_year) %>% 
  group_by(last_year, year) %>% 
  summarise(mae = mean(abs(pred - mean_opin)), 
            mape = mean(abs(pred - mean_opin)/mean_opin))

pl <- by_year_prediction %>% 
  filter(year >= last_year) %>% 
  group_by(last_year, issue) %>% 
  filter(year == max(year)) %>% 
  ggplot(aes(pred, mean_opin, label = issue, color = hfvl)) +
  geom_point() +
  facet_wrap(~last_year) +
  geom_abline()

plotly::ggplotly(pl)

extract_ranef <- function(mod){
  fixed_coef <- fixef(mod)
  ranef(mod)$issue %>% 
    as_tibble(rownames = "issue") %>% 
    left_join(mf_data %>% select(issue = code, hfvl = hfvl_advantage)) %>% 
    mutate(fitted = fixed_coef["year_r"] + fixed_coef["year_r:hfvl"]*hfvl) 
}

by_year_randef <- by_year %>% 
  mutate(rand_coef = map(mod, extract_ranef)) %>% 
  select(last_year, rand_coef) %>% 
  unnest(rand_coef)

by_year_randef <- by_year_randef %>% 
  mutate(trend = fitted + year_r)


by_year_randef %>% 
  ggplot(aes(year_r, fitted)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~last_year)
   

by_year_randef %>% 
  ggplot(aes(hfvl, trend)) +
  geom_point() +
  geom_smooth(method = "lm", color = "grey50") +
  geom_smooth() +
  facet_wrap(~last_year)

  


by_issue %>% 
  mutate(k = sd(2*slope_lin)/sd(hfvl),
         ak = k * hfvl) %>% 
  ggplot(aes(hfvl, 2*slope_lin)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_line(aes(y = ak))


m <- lm(2*slope_lin ~ hfvl, by_issue)
m2 <- lm(2*slope_lin_last ~ hfvl, by_issue)

sd(2*by_issue$slope_lin)/sd(by_issue$hfvl)
sd(2*by_issue$slope_lin_last)/sd(by_issue$hfvl)

by_issue %>% 
  ggplot(aes(hfvl, slope_lin))+
  geom_point() +
  geom_smooth(method = "lm")

by_issue %>% 
  ggplot(aes(hfvl, slope_lin_last))+
  geom_point() +
  geom_smooth(method = "lm")

by_issue <- by_issue %>% 
  mutate(tr_pred_18 = opin_2016 + 2*slope_lin, 
         k = sd(2*slope_lin)/sd(hfvl), 
         k_l = sd(2*slope_lin_last)/sd(hfvl),
         k_lm = coef(m)[2],
         k_lm2 = coef(m2)[2],
         arg_pred_18 = opin_2016 + k*hfvl,
         arg_pred_18lm = opin_2016 + k_lm*hfvl,
         arg_pred_18l = opin_2016 + k_l*hfvl,
         arg_pred_18lml = opin_2016 + k_lm2*hfvl)


by_issue <- gss_aggr %>%
  filter(year < 2018) %>% 
  group_nest(issue, hfvl)


by_issue <- by_issue %>% 
  mutate(data_sub = map(data, ~slice(., 1:(n() - 1))),
         data_last = map(data, ~slice(., (n() - 1):n())))

fit_mod <- function(data){
  glm(cbind(n_agree, n_disagree) ~ year_r2, 
      data, family = quasibinomial(link = "logit"))
}

fit_lin <- function(data){
  lm(mean_opin ~ year_r2, 
     data, weights = n_sample)
}

by_issue <- by_issue %>% 
  mutate(mod = map(data, fit_mod))

by_issue <- by_issue %>% 
  mutate(mod_lin = map(data, fit_lin))

by_issue <- by_issue %>% 
  mutate(mod_last = map(data_last, fit_mod))

by_issue <- by_issue %>% 
  mutate(mod_lin_last = map(data_last, fit_lin))


by_issue <- by_issue %>% 
  mutate(map_df(mod, coef)) %>% 
  rename(int = `(Intercept)`, slope = year_r2) %>% 
  mutate(map_df(mod_last, coef)) %>% 
  rename(int_last = `(Intercept)`, slope_last = year_r2) 

by_issue <- by_issue %>% 
  mutate(map_df(mod_lin, coef)) %>% 
  rename(int_lin = `(Intercept)`, slope_lin = year_r2) %>% 
  mutate(map_df(mod_lin_last, coef)) %>% 
  rename(int_lin_last = `(Intercept)`, slope_lin_last = year_r2)


opin2016 <- gss_aggr %>% 
  filter(year %in% c(2014, 2016, 2018)) %>% 
  mutate(measure = str_c("opin_", year)) %>% 
  select(issue, measure, mean_opin) %>% 
  spread(measure, mean_opin)

by_issue <- left_join(by_issue, opin2016)


by_issue %>% 
  mutate(k = sd(2*slope_lin)/sd(hfvl),
         ak = k * hfvl) %>% 
  ggplot(aes(hfvl, 2*slope_lin)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  geom_line(aes(y = ak))


m <- lm(2*slope_lin ~ hfvl, by_issue)
m2 <- lm(2*slope_lin_last ~ hfvl, by_issue)

sd(2*by_issue$slope_lin)/sd(by_issue$hfvl)
sd(2*by_issue$slope_lin_last)/sd(by_issue$hfvl)

by_issue %>% 
  ggplot(aes(hfvl, slope_lin))+
  geom_point() +
  geom_smooth(method = "lm")

by_issue %>% 
  ggplot(aes(hfvl, slope_lin_last))+
  geom_point() +
  geom_smooth(method = "lm")

by_issue <- by_issue %>% 
  mutate(tr_pred_18 = opin_2016 + 2*slope_lin, 
         k = sd(2*slope_lin)/sd(hfvl), 
         k_l = sd(2*slope_lin_last)/sd(hfvl),
         k_lm = coef(m)[2],
         k_lm2 = coef(m2)[2],
         arg_pred_18 = opin_2016 + k*hfvl,
         arg_pred_18lm = opin_2016 + k_lm*hfvl,
         arg_pred_18l = opin_2016 + k_l*hfvl,
         arg_pred_18lml = opin_2016 + k_lm2*hfvl)

by_issue %>% 
  mutate_at(vars(opin_2016, matches("pred_18")), ~ . - opin_2018) %>% 
#  select(matches("pred_18")) %>% 
  summarise_at(vars(opin_2016, matches("pred")), ~sqrt(mean(.^2))) %>% 
  gather() %>% 
  kableExtra::kbl(col.names = c(" ", "RMSE"), digits = 4)


by_issue %>% 
  ggplot(aes(fct_reorder(issue, opin_2018), opin_2018))+
  geom_point() +
  coord_flip() +
  geom_point(aes(y = tr_pred), color = "grey50") +
  geom_point(aes(y = arg_pred), color = "green")

by_issue %>% 
  select(issue, opin_2018, opin_2016,
         Trends  = tr_pred_18, 
         HFVL = arg_pred_18) %>%  
  gather("Predictor", "Prediction", Trends, HFVL) %>% 
  ggplot(aes(Prediction, opin_2018, color = Predictor))+
  geom_point(size = 3, alpha= .5) + 
#  geom_point(aes(x = opin_2016), color = "grey10") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline() +
  theme_classic() +
  ggthemes::scale_color_ptol()


by_issue %>% 
  ggplot(aes(opin_2016, opin_2018))+
  geom_point() +
  geom_abline()

rsq_tr <- cor(by_issue$slope_lin_last, by_issue$slope_lin)^2
rsq_arg <- cor(by_issue$slope_lin_last, by_issue$hfvl)^2
rsq_arg_old <- cor(by_issue$hfvl, by_issue$slope_lin)^2

pl_tr <-  by_issue %>% 
  ggplot(aes(slope_lin, slope_lin_last, label = issue)) +
  geom_point() +
  geom_text(size = 3) +
  geom_smooth(method = "lm") +
  ylim(c(-0.274, 0.6)) +
  annotate("text", x = -0.05, y = .4, label = sprintf("R sq. = %.2f", rsq_tr)) +
  theme_classic()+
  labs(x = "Trends estimated up to 2016", y = "Trends estimates between 2014 and 2016")

pl_arg <- by_issue %>% 
  ggplot(aes(hfvl, slope_lin_last, label = issue)) +
  geom_point() +
  geom_text(size = 3) +
  geom_smooth(method = "lm") +
  ylim(c(-0.274, 0.6)) +
  annotate("text", x = -0.2, y = .4, label = sprintf("R sq. = %.2f", rsq_arg)) +
  theme_classic()+
  labs(x = "HFVL argument advantage", y = "Trends estimates between 2014 and 2016")


pl_arg_old <- by_issue %>% 
  ggplot(aes(hfvl, slope_lin, label = issue)) +
  geom_point() +
  geom_text(size = 3) +
  geom_smooth(method = "lm") +
  ylim(c(-0.274, 0.6)) +
  annotate("text", x = -0.2, y = .2, label = sprintf("R sq. = %.2f", rsq_arg_old)) +
  theme_classic()+
  labs(x = "HFVL argument advantage", y = "Trends estimated up to 2016")


cowplot::plot_grid(pl_tr, pl_arg, pl_arg_old, nrow = 1)


gss_change_wide <- gss_aggr %>% 
  group_by(issue) %>% 
  mutate(diff = (mean_opin - lag(mean_opin)), 
         yrs = (year - lag(year)), 
         diff = diff*10/yrs) %>% 
  drop_na(diff) %>% 
  select(issue, year, hfvl, diff) %>% 
  spread(year, diff)
  



ggsave("figures/trends_for_2012-2018.jpeg", dpi = 350)

m_tr <- lm(slope_last ~ slope, by_issue)
m_slope <- lm(slope_last ~ hfvl, by_issue)
m_both <- lm(slope_last ~ slope +  hfvl, by_issue)


m_tr_lin <- lm(slope_lin_last ~ slope_lin, by_issue)
m_slope_lin <- lm(slope_lin_last ~ hfvl, by_issue)
m_both_lin <- lm(slope_lin_last ~ slope_lin +  hfvl, by_issue)

modelsummary::msummary(list(m_tr, m_slope, m_both, m_tr_lin, m_slope_lin, m_both_lin), 
                       statistic = "conf.int")

gss_aggr %>% 
  filter(issue %in% c("helpblk", "wrkwayup")) %>% 
  ggplot(aes(year, mean_opin, color = fct_reorder(issue, hfvl))) +
  geom_line()

gss_aggr %>% 
  filter(issue %in% c("colhomo", "spkhomo")) %>% 
  ggplot(aes(year, mean_opin, color = fct_reorder(issue, hfvl))) +
  geom_line()

# keep two last points
# by_issue_2 <- by_issue

m_tr2 <- lm(slope_last ~ slope, filter(by_issue_2, issue != "marwht"))
m_slope2 <- lm(slope_last ~ hfvl, filter(by_issue_2, issue != "marwht"))
m_both2 <- lm(slope_last ~ slope +  hfvl, filter(by_issue_2, issue != "marwht"))

m_tr_lin2 <- lm(slope_lin_last ~ slope_lin, filter(by_issue_2, issue != "marwht"))
m_slope_lin2 <- lm(slope_lin_last ~ hfvl, filter(by_issue_2, issue != "marwht"))
m_both_lin2 <- lm(slope_lin_last ~ slope_lin +  hfvl, filter(by_issue_2, issue != "marwht"))

m_tr2 <- lm(slope_last ~ slope, by_issue_2)
m_slope2 <- lm(slope_last ~ hfvl, by_issue_2)
m_both2 <- lm(slope_last ~ slope +  hfvl, by_issue_2)

m_tr_lin2 <- lm(slope_lin_last ~ slope_lin, by_issue_2)
m_slope_lin2 <- lm(slope_lin_last ~ hfvl, by_issue_2)
m_both_lin2 <- lm(slope_lin_last ~ slope_lin +  hfvl, by_issue_2)


modelsummary::msummary(list(m_tr2, m_slope2, m_both2, m_tr_lin2, m_slope_lin2, m_both_lin2), 
                       statistic = "conf.int")

modelsummary::msummary(list(m_tr_lin2, m_slope_lin2, m_both_lin2, m_tr_lin, m_slope_lin, m_both_lin), 
                       statistic = "conf.int")

modelsummary::msummary(list(m_tr2, m_slope2, m_both2, m_tr, m_slope, m_both), 
                       statistic = "conf.int")

rsq_arg <- cor(by_issue_2$slope_lin_last, by_issue_2$hfvl)^2
rsq_tr <- cor(by_issue_2$slope_lin_last, by_issue_2$slope_lin)^2

pl_arg <- by_issue_2 %>% 
  ggplot(aes(hfvl, slope_lin_last, label = issue)) +
  geom_point() +
  geom_text(size = 3) +
  geom_smooth(method = "lm") +
  annotate("text", x = -0.2, y = .4, label = sprintf("R sq. = %.2f", rsq_arg)) +
  theme_classic()+
  labs(x = "HFVL argument advantage", y = "Trends estimates between 2016 and 2018")


pl_tr <-  by_issue_2 %>% 
  ggplot(aes(slope_lin, slope_lin_last, label = issue)) +
  geom_point() +
  geom_text(size = 3) +
  geom_smooth(method = "lm") +
  annotate("text", x = -0.05, y = .4, label = sprintf("R sq. = %.2f", rsq_tr)) +
  theme_classic()+
  labs(x = "Trends estimated up to 2016", y = "Trends estimates between 2016 and 2018")

cowplot::plot_grid(pl_tr, pl_arg)
ggsave("figures/trends_for_2018.jpeg", dpi = 350)

by_issue_2 %>% 
  ggplot(aes(hfvl, slope_last, label = issue)) +
  geom_point() +
  geom_text(data = filter(by_issue_2, (hfvl > 0 & slope_last < 0)|(hfvl < 0 & slope_last > 0))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


other_drctn <- filter(by_issue_2, (hfvl > 0 & slope_last < 0)|(hfvl < 0 & slope_last > 0))

gss_aggr <- gss_aggr %>% 
   left_join( mf_data %>% select(issue = code, hfvl = hfvl_advantage))

gss_aggr %>% 
  filter(issue %in% other_drctn$issue, hfvl < -0.05) %>% 
  ggplot(aes(year, mean_opin, color = fct_reorder(issue, hfvl))) +
  geom_line() + 
  ylim(c(0, 1)) +
  geom_text(aes(label = round(hfvl, 2)), 
            data = filter(gss_aggr, 
                          year == 2018, 
                          issue %in% other_drctn$issue, 
                          hfvl < -0.05))

gss_aggr %>% 
  filter(issue %in% other_drctn$issue, hfvl > 00.05) %>% 
  ggplot(aes(year, mean_opin, color = fct_reorder(issue, hfvl))) +
  geom_line() + 
  ylim(c(0, 1)) +
  geom_text(aes(label = round(hfvl, 2)), 
            data = filter(gss_aggr, 
                          year == 2018, 
                          issue %in% other_drctn$issue, 
                          hfvl > 0.05))
  
  

summary(lm(slope_lin ~ hfvl, by_issue))
summary(m_slope)

modelsummary::msummary(list(m_tr, m_slope, m_both), 
                       statistic = "conf.int")



# With predictions --------------------------------------------------------

gss_aggr <- gss_bin %>% 
  drop_na(opinion) %>% 
  group_by(issue) %>% 
  mutate(min_year = min(year), 
         year_r = (year - min_year)/10) %>% 
  group_by(issue, min_year, year, year_r) %>%
  summarise(mean_opin = weighted.mean(opinion, wgt),
            n_agree = sum(opinion*wgt),
            n_sample = sum(wgt),
            n_disagree = n_sample - n_agree,
            .groups = "drop")

by_issue <- gss_aggr %>% 
  filter(year <= 2010) %>% 
  group_nest(issue, min_year) 


fit_mod <- function(data){
  glm(cbind(n_agree, n_disagree) ~ year_r, 
      data, family = quasibinomial(link = "logit"))
}

fit_lin <- function(data){
  lm(mean_opin ~ year_r, 
     data, weights = n_sample)
}


by_issue <- by_issue %>% 
  mutate(mod = map(data, fit_mod),
         mod_lin = map(data, fit_lin))


by_issue <- by_issue %>% 
  mutate(map_df(mod, coef)) %>% 
  rename(int = `(Intercept)`, slope = year_r) %>% 
  mutate(map_df(mod_lin, coef)) %>% 
  rename(int_lin = `(Intercept)`, slope_lin = year_r) 

by_issue <- left_join(by_issue, mf_data %>% select(issue = code, hfvl = hfvl_advantage))

by_issue %>% 
  ungroup() %>% 
  select(hfvl, slope, slope_lin) %>% 
  cor()

make_preditions <- function(data, int, slope) {
  data %>%
    mutate(pred = int + slope*year_r,
           pred = exp(pred)/(1 + exp(pred)))
}

newdata <- tibble(year = seq(2010, 2018, by = 2))
newdata <- newdata %>% 
  mutate(year_r = (year - 1973)/10)

summary(by_issue$mod_lin[[22]])
predict(by_issue$mod[[22]], newdata, type = "response", se.fit = TRUE)

by_issue <- by_issue %>% 
  mutate(newdata = list(newdata),
         newdata = map2(newdata, min_year, 
                        ~.x %>% mutate(year_r = (year - .y)/10)))

lin_pred_tr <- by_issue %>% 
  mutate(pred = map2(mod_lin, newdata, 
                     ~bind_cols(.y, 
                                predict(.x, .y, interval = "confidence") %>% 
                                  as_tibble()))) %>% 
  select(issue, hfvl, pred) %>% 
  unnest(pred)

m <- lm(slope_lin ~ hfvl, by_issue)
by_issue$exp_slope_lin <- predict(m)

change_slope <- function(model, exp_slope) {
  model$coefficients["year_r"] <- exp_slope
  model
}

lin_pred_arg <- by_issue %>% 
  mutate(mod_lin = map2(mod_lin, exp_slope_lin, change_slope),
         pred = map2(mod_lin, newdata, 
                     ~bind_cols(.y, 
                                predict(.x, .y, interval = "confidence") %>% 
                                  as_tibble()))) %>% 
  select(issue, hfvl, pred) %>% 
  unnest(pred)

comb_pred <- bind_rows(
  lin_pred_arg %>% mutate(predictor = "hfvl"),
  lin_pred_tr %>% mutate(predictor = "trends")
)

comb_pred <- comb_pred %>% 
  left_join(gss_aggr %>% select(issue, year, mean_opin)) %>% 
  drop_na(mean_opin) %>% 
  mutate(sq_er = (fit - mean_opin)^2)

comb_pred_ly <- comb_pred  %>% 
  group_by(issue) %>% 
  filter(year == max(year)) %>% 
  ungroup() 

comb_pred_ly %>% 
  group_by(predictor) %>% 
  summarise(cor(fit, mean_opin))

comb_pred %>% 
  group_by(predictor) %>% 
  summarise(cor(fit, mean_opin))

comb_pred_ly %>% 
  group_by(predictor) %>% 
  summarise(mean(sq_er))

comb_pred %>% 
  group_by(predictor) %>% 
  summarise(mean(sq_er))

comb_pred_ly %>% 
  select(issue, predictor, hfvl_adv = hfvl, sq_er) %>% 
  spread(predictor, sq_er) %>% 
  count(hfvl < trends)


comb_pred_ly %>% 
  ggplot(aes(fct_reorder(issue, hfvl), fit, 
             ymin = lwr, ymax = upr, 
             color = predictor))+
  geom_pointrange() + 
  geom_point(aes(y = mean_opin), color = "red") +
  coord_flip() +
  scale_color_grey()

by_issue <- by_issue %>% 
  mutate(data = map(data, rename, pred_tr = pred),
         data = pmap(list(data, int, pred_slope), make_preditions),
         data = map(data, rename, pred_hfvl = pred))

gss_pred <- by_issue %>% 
  select(issue, hfvl, data) %>% 
  unnest(data) %>% 
  ungroup() 

gss_pred %>% 
  mutate(pred_mean = (pred_tr + pred_hfvl)/2) %>% 
  group_by(issue) %>% 
  summarise(trends_r = cor(mean_opin, pred_tr), 
            hfvl_r = cor(mean_opin, pred_hfvl),
            mean_r = cor(mean_opin, pred_mean)) %>% 
  summarise_at(vars(-issue), mean)


gss_pred <- gss_pred %>% 
  mutate(pred_tr_sq = (mean_opin - pred_tr)^2,
         pred_hfvl_sq = (mean_opin - pred_hfvl)^2, 
         pred_mean = (pred_tr + pred_hfvl)/2,
         pred_mean_sq = (mean_opin - pred_mean)^2)

gss_pred %>% 
  summarise(sqrt(mean(pred_tr_sq)), 
            sqrt(mean(pred_hfvl_sq)), 
            sqrt(mean(pred_mean_sq)))

gss_pred %>% 
  group_by(issue, hfvl) %>% 
  summarise(rmse_tr = sqrt(mean(pred_tr_sq)), 
            rmse_arg = sqrt(mean(pred_hfvl_sq)), 
            rmse_mean = sqrt(mean(pred_mean_sq)),
            .groups = "drop") %>% 
  mutate(arg_better = rmse_tr - rmse_arg, 
         mean_better = rmse_tr - rmse_mean) %>%
  summarise(cor(abs(hfvl), arg_better), 
            cor(abs(hfvl), mean_better))


gss_pred %>% 
  group_by(issue, hfvl) %>% 
  summarise(rmse_tr = sqrt(mean(pred_tr_sq)), 
            rmse_arg = sqrt(mean(pred_hfvl_sq)), 
            rmse_mean = sqrt(mean(pred_mean_sq)),
            .groups = "drop") %>% 
  mutate(arg_better = rmse_tr - rmse_arg, 
         mean_better = rmse_tr - rmse_mean) %>%
  ggplot(aes(hfvl, arg_better)) +
  geom_point() +
  geom_text(aes(label = issue), size = 3) +
  geom_hline(yintercept = 0)

gss_pred %>% 
  group_by(issue, hfvl) %>% 
  summarise(rmse_tr = sqrt(mean(pred_tr_sq)), 
            rmse_arg = sqrt(mean(pred_hfvl_sq)), 
            rmse_mean = sqrt(mean(pred_mean_sq)),
            .groups = "drop") %>% 
  mutate(arg_better = rmse_tr - rmse_arg, 
         mean_better = rmse_tr - rmse_mean) %>%
  ggplot(aes(hfvl, mean_better)) +
  geom_point() +
  geom_text(aes(label = issue), size = 3) +
  geom_hline(yintercept = 0)


gss_pred %>% 
  filter(issue %in% c("homosex","grass", "cappun", "suicide1", "fefam")) %>% 
  pivot_longer(c(mean_opin, pred_tr, pred_hfvl, pred_mean), 
               names_to = "var", 
               values_to = "value") %>% 
  ggplot(aes(year, value, color = var)) +
  geom_line() + 
  facet_wrap(~issue)


gss_pred %>% 
  filter(issue %in% c("marhomo","fehelp", "marasian", "marhisp", "pilok", "racseg")) %>% 
  pivot_longer(c(mean_opin, pred_tr, pred_hfvl, pred_mean), 
               names_to = "var", 
               values_to = "value") %>% 
  ggplot(aes(year, value, color = var)) +
  geom_line() + 
  facet_wrap(~issue)

