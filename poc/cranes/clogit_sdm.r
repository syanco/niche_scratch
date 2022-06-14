# scratch script to play with using clogit to produce a dynamic sdm like thingy

library(tidyverse)
library(rstoat)
library(lubridate)
library(survival)
library(glue)

source("src/funs/big_stoat.r")

dat <- read.csv("data/542803025.csv")

dat1 <- dat %>%
  mutate(
    date = as.Date(timestamp, format = "%m/%d/%y %H:%M"),
    wk = format(date, format="%U"),
    yr = year(date)
  ) %>%  
  filter (yr == 2019)


# create conditional dataset
wks <- unique(dat1$wk)

out <- list()

for(i in 1:length(wks)){
  out[[i]] <- dat1 %>% 
    arrange(wk) %>% 
    mutate(strat = wks[i],
           use = case_when(wk == wks[i] ~ 1,
                           T ~ 0)) 
  
}

dat_strat <- do.call("rbind", out) %>% 
  mutate(d = glue("{yr}-{strat}-4"),
         date_int = as_date(d, format = "%Y-%U-%u"),
         lat_r = round(lat, 3), # round coords to reasonable spatial precision
         lon_r = round(lon, 3))

# reduce dataset to unique date X loc combinations to speed up annotation
dat_stoat_in <- dat_strat %>% 
  distinct(lon_r, lat_r, date_int) %>% 
  rename(date = date_int, lng = lon_r, lat = lat_r ) %>% 
  filter(!is.na(date))

# Declare annotation variables
vars <- list(list(id = "JRC/GSW1_3/MonthlyHistory",
                  static = FALSE,
                  reducers = list("mean"),
                  s_buff = 30.1, # i dunno
                  
                  t_buff = 30, 
                  bands = list("water")))

# run stoat annotation
dat_stoat_out <- big_stoat(data = dat_stoat_in, vars = vars)

# join annotation back up to the full dataset

dat_clog <- dat_strat %>% 
  full_join(dat_stoat_out, 
            by = c("lon_r" = "lng", "lat_r" = "lat", "date_int" = "date")) %>% 
  mutate(use = as.integer(use)) %>% 
  filter(!is.na(value)) %>% 
  mutate(
    val_rnd = round(value, 0),
    water = case_when(val_rnd == 0 ~ NA_real_,
                      val_rnd == 1 ~ 0,
                      val_rnd == 2 ~ 1)) %>% 
  filter(water != 0)

ggplot(dat_clog) +
  geom_density(aes(x = water, color = as.factor(use)))

fm <- clogit(use ~ water + strata(strat), data = dat_clog, method = "approximate")
summary(fm)
summary(resid(fm))

ACF(fm)
exp(fm$coefficients)
exp(confint(fm))
plogis(fm$coefficients)
plogis(confint(fm))

##--  SAVE/LOAD  --##
# save.image(file = "analysis/cranes/clog_scratch_06082022.rdata")
load(file = "analysis/cranes/clog_scratch_06082022.rdata")
