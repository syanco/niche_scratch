library(gtools)

ind_ts_dat_rev <- ind_ts_dat %>% 
  group_by(ind) %>% 
  mutate(ts_imp = 1:n())


ggplot(ind_ts_dat %>% ungroup(), aes(x = ts, y = c_range))+
  geom_line(aes(group=ind))+
  facet_wrap(~ind) +
  theme_minimal()

ggplot(ind_ts_dat %>% ungroup(), aes(x = c_dist, y = c_range))+
  geom_line(aes(group=ind))+
  theme_minimal() +
  facet_wrap(~ind)

ggplot(ind_ts_dat %>% ungroup(), aes(x = c_n, y = c_range))+
  geom_line(aes(group=ind))+
  theme_minimal() +
  facet_wrap(~ind)


ggplot(ind_ts_dat_rev , aes(x = ts_imp, y = c_range))+
  geom_point(aes(color=as.factor(ind)))+
  theme_minimal()

ggplot(ind_ts_dat_rev , aes(x = c_dist, y = c_range))+
  geom_point(aes(color=as.factor(ind)))+
  theme_minimal()

ggplot(ind_ts_dat_rev , aes(x = c_n, y = c_range))+
  geom_point(aes(color=as.factor(ind)))+
  theme_minimal()

# Logistic model
#get inits
coef(lm(logit(c_range/2)~c_dist,data=ind_ts_dat_rev))





ggplot(ind_ts_dat %>% ungroup(), aes(x = ts, y = c_dist))+
  geom_point()+
  theme_minimal()
