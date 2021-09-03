dissim_mean <- dissim_ts %>%
  mutate(week = week(ts)) %>% 
  filter(BD_tot != -Inf,
         BD_tot != Inf) %>% 
  group_by(week) %>% 
  summarize(BD_mu = mean(na.omit(BD_tot)))

size_mean <- n_size %>% 
  mutate(week = week(ts)) %>% 
  filter(total != -Inf,
         total != Inf) %>% 
  group_by(week) %>% 
  summarize(size_mu = mean(na.omit(total)))

nsd_wk_mu <- nsd %>% 
  group_by(week) %>% 
  summarize(nsd = sum(netSQ),
            vel = sum(abs(vel)),
            ind = individual_id[1],
            ts = timestamp[1]) %>%
  mutate(ts = ymd_hms(ts),
         ind = as.factor(ind))

dmu <- ggplot(dissim_mean) +
  geom_line(aes(x=week, y= BD_mu))

smu <- ggplot(size_mean) +
  geom_line(aes(x=week, y= size_mu))

nsdmu <- ggplot(nsd_wk_mu) +
  geom_line(aes(x=week, y= nsd))

ggarrange(dmu, smu, nsdmu, ncol = 1)

comb_mean <- dissim_mean %>% 
  left_join(size_mean) %>% 
  left_join(nsd_wk_mu)

ggplot(comb_mean, aes(x=nsd, y = BD_mu)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_mean, aes(x=nsd, y = size_mu)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_mean, aes(x=vel, y = BD_mu)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_mean, aes(x=vel, y = size_mu)) +
  geom_point() +
  geom_smooth(method = "lm")
#############################################

dissim_ts <- ind_ts_dat %>% 
  # filter(ind == idx[i]) %>%    
  mutate(ts = ymd_hms(ts),
         ind = as.factor(ind)) %>% 
  group_by(ind, WKxYR) %>% 
  summarize(ind = ind,
            ts = ts,
            BD_tot = Bhattacharyya_distance[[1]][1],
            day_temp= Bhattacharyya_distance[[1]][2],
            evi = Bhattacharyya_distance[[1]][3],
            ndvi = Bhattacharyya_distance[[1]][4]) %>% 
  filter(BD_tot != -Inf,
         BD_tot != Inf)

n_size <- n_size %>% 
  mutate(ind = as.factor(ind)) %>% 
  filter(total != -Inf,
         total != Inf)

nsd_wk <- nsd %>% 
  group_by(individual_id, WKxYR) %>% 
  mutate(n=n()) %>% #calculate # of observations per interval
  filter(n >= 2) %>% 
  summarize(nsd = sum(netSQ),
            vel = sum(abs(vel)),
            ind = as.factor(individual_id[1]),
            ts = timestamp[1]) %>%
  mutate(ts = ymd_hms(ts),
         ind = as.factor(ind),
         WKxYR = as.character(WKxYR))

comb_wk <- dissim_ts %>% 
  left_join(n_size, by = c("WKxYR", "ind")) %>% 
  left_join(nsd_wk, by = c("WKxYR", "ind"))

ggplot(comb_wk, aes(x=vel, y = total))+
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_wk, aes(x=nsd, y = total))+
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_wk, aes(x=vel, y = BD_tot))+
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_wk, aes(x=nsd, y = BD_tot))+
  geom_point() +
  geom_smooth(method = "lm")



ggplot(comb_wk, aes(x=WKxYR, y = BD_tot, group = ind))+
  geom_line()

lmer(total ~ vel + (1|ind), data = comb_wk)

# +
#   geom_smooth(method = "lm")

ggplot(comb_wk, aes(x=WKxYR, y = BD_tot))+
  geom_point() +
  geom_smooth(method = "lm")

ggplot(comb_wk, aes(x=WKxYR, y = BD_tot))+
  geom_point() +
  geom_smooth(method = "lm")
