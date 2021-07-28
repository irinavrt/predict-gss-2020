library(tidyverse)
library(haven)
library(labelled)

# The file GSS7218_R3.sav is available at ...

# GSS ---------------------------------------------------------------------

gss_full <- read_sav("../data/GSS7218_R3.sav")

issue_items <- issue_list$issue[issue_list$asked_after_2010 == 1]

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
  select(id, year, wtssall, formwt, sample, oversamp,
         polviews, partyid, wordsum, degree, educ, 
         talkpol, talkpol1, talkpol2, talkpol3,
         discpol, poldisgn, chat12, polinf12,
         poleff13, poleff15, polefy13, polefy15,
         sex, age, race, class, region,  finrela, income, rincome, income16, income, rincome, income16,
         relig, attend, god, reliten,
         one_of(issue_items))

gss <- gss %>% 
  mutate(othshelp_dub = othshelp,
         careself_dub = careself) %>% 
  mutate_at(vars(id, wtssall, oversamp, age, educ, year, polviews,
                 othshelp_dub, careself_dub, 
                 peoptrbl:helpjob), zap_labels) %>% 
  mutate_if(is.labelled, ~fct_relabel(as_factor(.), tolower)) %>% 
  mutate(birth_year = year - age,
         polviews_cont = polviews,
         polviews = cut(polviews, c(0, 3, 4, 7), 
                        labels = c("liberal", "moderate", "conservative")),
         year_r = (year - 1972)/10,
         wgt = wtssall*oversamp) %>% 
  mutate_at(vars(othshelp_dub, peoptrbl),
            funs((5 - .)/4)) %>% 
  mutate_at(vars(givblood:helpjob),
            funs((6 - .)/5)) %>% 
  mutate(selffrst = (selffrst - 1)/4,
         careself_dub = (careself_dub - 1)/4,
         selfless = (6 - selfless)/5)

gss <- gss %>% 
  mutate(oth_att = rowMeans(select(., othshelp_dub, peoptrbl, selffrst)),
         oth_att_full = rowMeans(select(., othshelp_dub, careself_dub, peoptrbl, selffrst)),
         oth_behave = rowMeans(select(., givhmlss, retchnge, cutahead, 
                                  volchrty, givseat,  
                                  carried, directns, loanitem, helphwrk, 
                                  lentto)),
         oth_behave_full = rowMeans(select(., givblood:helpjob)),
         oth_netw_b = rowMeans(select(., helphwrk, lentto, talkedto, helpjob)),
         oth_strng_b = rowMeans(select(., givhmlss, retchnge, cutahead, 
                                       volchrty, givseat,  helpaway,
                                       carried, directns, loanitem))) %>% 
  # none did all behaviors, devide by max to standardize
  mutate_at(vars(oth_behave, oth_behave_full, oth_netw_b, oth_strng_b),
            ~./max(., na.rm = TRUE))


# Recoding levels that are not used
gss <- gss %>% 
  mutate(class = fct_recode(class, NULL = "no class"),
         homosex = fct_recode(homosex, NULL = "other"),
         racchng = fct_recode(racchng, NULL = "wdnt belong"),
         racopen = fct_recode(racopen, NULL = "neither"),
         sexeduc = fct_recode(sexeduc, NULL = "depends"),
         relig = fct_lump(relig, prop = .05),
         reliten = fct_relevel(reliten, "no religion", "not very strong", 
                               "somewhat strong", "strong")) 

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
  select(id, year, year_r, wgt,
         polviews, polviews_cont, partyid, 
         wordsum, degree, educ, 
         oth_att, oth_behave, oth_att_full, oth_behave_full, 
         oth_netw_b, oth_strng_b,
         matches("talkpol"), discpol, poldisgn,
         poleff13, poleff15, polefy13, polefy15,
         relig, attend, god, reliten,
         sex, age, birth_year, race, class, region, finrela, income, rincome, income16, 
         one_of(issue_items)) %>% 
  mutate(pornlaw = ifelse(pornlaw == "legal", 0 , 1), 
         paidlvdv = ifelse(as.numeric(paidlvdv) < 3, 0, 1)) %>% 
  mutate_at(issue_items[!issue_items  %in% c("pornlaw", "paidlvdv")], dichotomize)

# add proportion of liberals
gss_bin <- gss_bin %>% 
  group_by(region) %>% 
  mutate(liberals_prop = weighted.mean(polviews == "liberal",
                                       na.rm = TRUE,
                                       wgt)) %>% 
  ungroup() 


gss_bin <- gss_bin %>% 
  gather(issue, opinion, abany:givinffor) 

gss_aggr <- gss_bin %>% 
  drop_na(opinion) %>%
  group_by(issue, year) %>%
  summarise(mean_opin = weighted.mean(opinion, wgt),
            n_agree = sum(opinion*wgt),
            n_sample = sum(wgt),
            n_disagree = n_sample - n_agree,
            .groups = "drop")

write_rds(gss_aggr, "data/gss-issp-trends-to-predict.rds")


