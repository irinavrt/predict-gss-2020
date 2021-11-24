library(tidyverse)
library(haven)
library(labelled)

# The file GSS7218_R3.sav is available at ...

# GSS ---------------------------------------------------------------------

gss_full <- read_sav("data/GSS7218_R3.sav")
issue_list <- read_csv("data/issue-list.csv")

# Dropped unused issues from the issue list (not asked in 2010)

gss <- gss_full %>%
  rename_all(tolower) 

gss <- gss %>% 
  mutate_at(vars(hubbywk1, twoincs1, homosex1, premars1, xmarsex1), 
            remove_val_labels) %>% 
  # combine two different wordings of the same issues
  mutate(hubbywrk = coalesce(hubbywrk, hubbywk1),
         twoincs = coalesce(twoincs, twoincs1),
         homosex = coalesce(homosex, homosex1),
         premars = coalesce(premarsx, premars1),
         xmarsex = coalesce(xmarsex, xmarsex1)) %>% 
  select(id, year, wtssall, formwt, oversamp,
         one_of(issue_list$issue))

gss <- gss %>% 
  mutate_at(vars(id, wtssall, oversamp), zap_labels) %>% 
  mutate_if(is.labelled, ~fct_relabel(as_factor(.), tolower)) %>% 
  mutate(wgt = wtssall*oversamp) 


# Recoding levels that are not used
gss <- gss %>% 
  mutate(homosex = fct_recode(homosex, NULL = "other"),
         racopen = fct_recode(racopen, NULL = "neither"),
         sexeduc = fct_recode(sexeduc, NULL = "depends"))

#Excluding years where racial questions were asked of non-blacks only 
for(i in c("racmar","racpush", "racopen","racdin")) {
  gss[gss$year %in% 1972:1977,i] <- NA
}

gss <- gss %>% droplevels()

#dichotomize variable

dichotomize <- function(var){
  n_categories = n_distinct(var, na.rm = TRUE)
  # for already binary items recode "yes" (1st level) to 1 and "no" to 0
  if(n_categories == 2){
    return(2 - as.numeric(var))
  }
  # if item has even number of levels recode first half of levels as 1, and second as 0
  if(n_categories %% 2 == 0){
    return(ifelse(as.numeric(var) > n_categories/2, 0, 1))
  } else {
    # if item has odd number of levels recode the middle level to NA, first half of levels as 1, and second as 0
    middle <- ceiling(n_categories/2)
    return(case_when(
      as.numeric(var) == middle ~ NA_real_,
      as.numeric(var) < middle ~ 1,
      as.numeric(var) > middle ~ 0,
      TRUE ~ NA_real_)
    )
  }
}

gss_bin <- gss %>% 
  select(id, year, wgt,
         one_of(issue_list$issue)) %>% 
  mutate(pornlaw = ifelse(pornlaw == "legal", 0 , 1), 
         paidlvdv = ifelse(as.numeric(paidlvdv) < 3, 0, 1),
         # reverse code the items to match the AA default position
         across(c(wrksch, wrkbaby), ~ case_when(
           . == "stay home" ~ 1,
           . == "work full-time" ~ 0, 
           TRUE ~ NA_real_
         ))) %>% 
  mutate_at(issue_list$issue[!issue_list$issue  %in% c("pornlaw", "paidlvdv", "wrksch", "wrkbaby")], 
            dichotomize)

# Save 2018 sample for analysis of the sampling error
gss_2018 <- gss_bin %>% 
  filter(year == 2018) %>% 
  select(id, year, wgt, abany:wrksch)

write_rds(gss_2018, "data/gss-issp-individual-data-2018.rds")

gss_bin <- gss_bin %>% 
  gather(issue, opinion, abany:wrksch) 

gss_aggr <- gss_bin %>% 
  drop_na(opinion) %>%
  group_by(issue, year) %>%
  summarise(mean_opin = weighted.mean(opinion, wgt),
            .groups = "drop")

write_rds(gss_aggr, "data/gss-issp-trends-to-predict.rds")


# 2021 GSS ----------------------------------------------------------------

gss_21 <- read_stata("data/gss2021.dta")

gss_21 <- gss_21 %>% 
  # rename the items without volunteered response to match names from previous GSS
  rename(racopen = racopennv, grass = grassnv) %>% 
  # drop police items for the experimental Y form (with alternative gender neutral wording)
  mutate(polmurdr = ifelse(form == 1, polmurdr, NA_real_),
         polescap = ifelse(form == 1, polescap, NA_real_)) %>% 
  select(id, year, wgt = wtssps, one_of(issue_list$issue))

gss_21_bin <- gss_21 %>% 
  mutate(pornlaw = ifelse(pornlaw == 3, 0, 1), 
         racopen = na_if(racopen, 3),
         across(c(abany:polmurdr, premarsx:xmarsex), dichotomize))

write_rds(gss_21_bin, "data/gss-issp-individual-data-2021.rds")

gss_21_bin <- gss_21_bin %>% 
  gather(issue, opinion, abany:xmarsex) 

gss_21_aggr <- gss_21_bin %>% 
  drop_na(opinion) %>%
  group_by(issue, year) %>%
  summarise(mean_opin = weighted.mean(opinion, wgt),
            mean_opin_no_wgt = mean(opinion),
            .groups = "drop")

write_rds(gss_21_aggr, "data/gss-issp-2021.rds")

