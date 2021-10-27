# Use the `niche_mix` conda environment

library(forcats)
library(lubridate)
library(dplyr)
library(RSQLite)
library(DBI)
library(gridExtra)
library(ggplot2)


#---- Functions ----#

# Anticipates vector of sample means, vector of sample variances
estPopVar <- function(sd, means){
  pop_mean <- mean(means)
  n <- length(sd)
  
  vec <- c()
  
  for(i in 1:n){
    vec[i] <- (1/n)*(sd[i]^2 + means[i]^2 - pop_mean^2)
  }
  
  out <- sum(vec)
  return(out)
}

# calculate an individual's total contribution to pop variance 
indContrib <- function(mu_i, sd_i, mu_pop, n){
  cont <- (1/n)*(sd_i^2 + mu_i^2 - mu_pop^2)
  return(cont)
}

# calculate the mean component of ind contribution to pop variance
#(variance contrib is just the ind variance itself)
muContrib <- function(mu_i, mu_pop, n){
  cont <- (1/n)*(mu_i^2 - mu_pop^2)
  return(cont)
}

#---- Load data ----#

# Pull in Pop scale niche data
pop <- read.csv("analysis/cranes/Anthropoides virgo_landsat8_evi_500_16.csv")

# Pull in indivdiual niches

.anno <- "anno_join_2021-09-08"
.wd <- '~/projects/niche_scratch'
.dbPF <- file.path(.wd,'data/anno_move.db')

db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)

message("Loading annotations...")
anno0 <- tbl(db, .anno)
message("Loading species data...")
sp <- tbl(db, "individual")

#read in NSD data created by `calc_nsd.r`
dat_track_nsd <- read.csv("analysis/cranes/nsd.csv")

#---- Perform analysis  ----#

#get species smaple sizes
sp %>% 
  group_by(taxon_canonical_name) %>% 
  summarise(n= n())

dat <- anno0 %>% 
  left_join(sp, by = "individual_id") %>% 
  filter(taxon_canonical_name == "Anthropoides virgo") %>%
  collect() %>% 
  mutate(ind_f = factor(individual_id))


# Calc individual means and vars
ind_sum <- dat %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(ind_f) %>% 
  summarise(mu = mean(na.omit(`value_derived:evi`)),
            sd = sd(na.omit(`value_derived:evi`)))

# make population estimates
(mix_mean <-mean(na.omit(ind_sum$mu))) 
(mix_var <- estPopVar(sd = na.omit(ind_sum$sd), means = na.omit(ind_sum$mu)))

# calc pop estimates from MOL data

#pop mean
(pop_mean <- mean(na.omit(pop$value)))
(pop_var <- var(na.omit(pop$value)))


# do it with lmer
library(lme4)

fm <- lmer(`value_derived:evi` ~ 1 + (1|ind_f), data = dat)
summary(fm)
(vc <- as.data.frame(VarCorr(fm)))
vc_var <- sum(vc$vcov)
vc_mean <- fixef(fm)

# gather into a df
df <- data.frame(mu = c(mix_mean, vc_mean, pop_mean),
                 var = c(mix_var, vc_var, pop_var))
row.names(df) <- c("mixture est", "variance components est", "pop est")
df

# get individual contributions to pop variance
contrib <- c()
mean_con <- c()
for(i in 1:nrow(ind_sum)){
  contrib[i] <- indContrib(ind_sum$mu[i], ind_sum$sd[i], mu_pop = pop_mean, n = nrow(ind_sum))
  mean_con[i] <- muContrib(ind_sum$mu[i], mu_pop = pop_mean, n = nrow(ind_sum))
}
ind_sum$tot_contrib <- contrib
ind_sum$mean_contrib <- mean_con

#check that the variance contribution components
round(((1/nrow(ind_sum))*ind_sum$sd^2)+ind_sum$mean_contrib, digits = 17) == round(ind_sum$tot_contrib, digits = 17)
# they match to within 17 decimals...

# Pull some movement stats for each ind
move_sum <- dat_track_nsd %>% 
  mutate(ind_f = as.factor(individual_id),
         ts = ymd_hms(timestamp)) %>% 
  group_by(ind_f) %>% 
  summarise(max_nsd = max(na.omit(netSQ)),
            max_v = max(na.omit(vel)),
            mean_v = mean(na.omit(vel)),
            mean_lat = mean(y_),
            mean_lon = mean(x_), 
            dur = difftime(max(ts), min(ts), units = "days")) %>% 
  right_join(ind_sum) %>% 
  filter(dur >= 200) #remove individuals with "short" tracks

dat <- dat %>% 
  mutate(ts = ymd_hms(timestamp),
         doy = yday(ts),
         yr = year(ts),
         indyr = factor(paste0(ind_f, yr)))


#---- Plots ----#
library(ggplot2)
library(viridisLite)

ggplot() +
  geom_line(data=dat, aes(x = doy, y = `value_derived:evi`, color = indyr),
            alpha= 0.3) +
  theme(legend.position = "none")


# plot of all ind and pop distributions
ggplot() +
  geom_line(data = dat , aes(x=`value_derived:evi`, group = individual_id, color = individual_id), alpha = 0.2,
            stat = "density") +
  geom_density(dat = pop, aes(x=value), size = 2)



#plot individual means (like a coefficient plot)
move_sum %>% mutate(ind_f = fct_reorder(ind_f, mu, min)) %>%   # reset factors
  ggplot()+
  geom_point(aes(x=mu, y = ind_f, color = tot_contrib))+
  scale_color_viridis_c() +
  geom_vline(data=pop, aes(xintercept=mean(na.omit(value))))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#plot individual vars (like a coefficient plot)
move_sum %>% mutate(ind_f = fct_reorder(ind_f, sd^2, min)) %>%   # reset factors
  ggplot()+
  geom_point(aes(x=sd^2, y = ind_f, color = tot_contrib))+
  scale_color_viridis_c() +
  geom_vline(data=pop, aes(xintercept=var(na.omit(value))))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# plot variance contribution components
ggplot(move_sum)+
  geom_point(aes(x=mean_contrib, y = sd^2, color = mean_lat)) +
  geom_abline(intercept = 0)
  # xlim(c(0,0.05)) +
  # ylim(c(0,0.05))


# ggplot(move_sum, aes(x=mean_lat, y = mean_contrib))+
#   geom_point() +
#   geom_smooth(method= "lm")
# 
# lat0 <- lm(mean_contrib ~ 1, data = move_sum)
# lat1 <- lm(mean_contrib ~ mean_lat, data = move_sum)
# AIC(lat0, lat1)
# summary(lat1)
# 
# ggplot(move_sum, aes(x=mean_lat, y = sd^2))+
#   geom_point() +
#   geom_smooth(method= "lm")
# 
# lat0 <- lm(sd^2 ~ 1, data = move_sum)
# lat1 <- lm(sd^2 ~ mean_lat, data = move_sum)
# AIC(lat0, lat1)
# 


# Movement Plots
ggplot(move_sum)+
  geom_point(aes(x=max_nsd, y = tot_contrib))
ggplot(move_sum)+
  geom_point(aes(x=max_nsd, y = mean_contrib))
ggplot(move_sum)+
  geom_point(aes(x=max_nsd, y = sd^2))
ggplot(move_sum)+
  geom_point(aes(x=max_nsd, y = mean_contrib/(sd^2)))
ggplot(move_sum)+
  geom_point(aes(x=mean_lat, y = mean_contrib/(sd^2)))+
  geom_smooth(aes(x=mean_lat, y = mean_contrib/(sd^2)), method = "lm")
  
summary(lm(mean_contrib/sd^2 ~ mean_lat, data = move_sum))
ggplot(move_sum)+
  geom_point(aes(x=mean_lat, y = tot_contrib))
ggplot(move_sum)+
  geom_point(aes(x=mean_lon, y = tot_contrib))
ggplot(move_sum)+ 
  geom_point(aes(x=mean_lon, y = mean_lat))
ggplot(move_sum)+
  geom_point(aes(x=dur, y = max_nsd))


ggplot(move_sum, aes(x=max_nsd)) +
  geom_density() +
  geom_dotplot()


#---- Intra-Individual Analysis ----#
inds <- unique(dat$ind_f)

out <- list()

for (j in 1:length(inds)){
  targ <- inds[j]
  
  ind_dat <-dat %>% 
    filter(ind_f == targ) %>% 
    mutate(week = week(ts))
  
  # Calc weekly means and vars
  targ_sum <- ind_dat %>%
    mutate(ind_f = as.factor(individual_id),
           week_f = as.character(week)) %>%
    group_by(week_f) %>% 
    summarise(mu = mean(na.omit(`value_derived:evi`)),
              sd = sd(na.omit(`value_derived:evi`)),
              ind_f = ind_f[1],
              week = week[1])
  
  # make ind estimates
  (mix_mean <-mean(na.omit(targ_sum$mu))) 
  (mix_var <- estPopVar(sd = na.omit(targ_sum$sd), means = na.omit(targ_sum$mu)))
  
  # individual direct calc
  (mu_ind <- mean(na.omit(ind_dat$`value_derived:evi`)))
  (sd_ind <- var(na.omit(ind_dat$`value_derived:evi`)))
  
  # get weekly contributions to ind variance
  contrib <- c()
  mean_con <- c()
  for(i in 1:nrow(targ_sum)){
    contrib[i] <- indContrib(targ_sum$mu[i], targ_sum$sd[i], mu_pop = mu_ind, n = nrow(targ_sum))
    mean_con[i] <- muContrib(targ_sum$mu[i], mu_pop = mu_ind, n = nrow(ind_sum))
  }
  targ_sum$tot_contrib <- contrib
  targ_sum$mean_contrib <- mean_con
  targ_sum$week_f <- as.factor(targ_sum$week)
  
  # Pull some movement stats for targ ind
  move_sum_ind <- dat_track_nsd %>% 
    mutate(ind_f = as.character(individual_id),
           week_f = as.character(week),
           ts = ymd_hms(timestamp)) %>% 
    filter(ind_f == targ) %>% 
    group_by(week_f) %>% 
    summarise(max_nsd = max(na.omit(netSQ)),
              max_v = max(na.omit(vel)),
              mean_v = mean(na.omit(vel)),
              mean_lat = mean(y_),
              mean_lon = mean(x_), 
              dur = difftime(max(ts), min(ts), units = "days")) %>% 
    right_join(targ_sum)
  
  out[[j]] <- move_sum_ind
  print(glue("{j} out of {length(inds)}"))
}


out_tot <- do.call("rbind", out)

p1 <- ggplot(data = out_tot) +
  geom_point(aes(x=week, y= mean_contrib))

p2 <- ggplot(data = out_tot) +
  geom_point(aes(x=week, y= sd^2))

p3 <- ggplot(data = out_tot) +
  geom_line(aes(x=week, y= max_nsd, color = ind_f)) +
  theme(legend.position = "none")

# combine plots
comb_plot <- grid.arrange(p1, p2, p3,
                       ncol = 1)


ggplot(data = out_tot)+
  geom_point(aes(x = mean_v, y = sd^2))

  # filter(dur >= 200) #remove individuals with "short" tracks

pop_df <- data.frame(mu = pop_mean,
                     var = pop_var)

#plot individual means (like a coefficient plot)
out_tot %>% mutate(week_f = fct_reorder(week_f, mean_contrib, min)) %>%   # reset factors
  ggplot()+
  geom_line(aes(x=mu, y = week_f))+
  scale_color_viridis_c() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


#plot individual vars (like a coefficient plot)
targ_sum %>% mutate(week_f = fct_reorder(week_f, sd^2, min)) %>%   # reset factors
  ggplot()+
  geom_point(aes(x=sd^2, y = week_f, color = tot_contrib))+
  scale_color_viridis_c() +
  geom_vline(data=pop_df, aes(xintercept=var(na.omit(`value_derived:evi`))))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


ggplot(data = out_tot)+
  geom_line(aes(x=week, y = tot_contrib, group = ind_f))+
  geom_line(aes(x=week, y = sd^2, group = ind_f), color = "red")+
  geom_line(aes(x=week, y = mean_contrib, group = ind_f), color = "blue")

out_mean <- out_tot %>% 
  group_by(week) %>% 
  summarize(tot_contrib = mean(na.omit(tot_contrib)),
            var = mean(na.omit(sd^2)),
            mean_contrib = mean(na.omit(mean_contrib)))
ggplot(data = out_mean)+
  geom_point(aes(x=week, y = tot_contrib))+
  geom_smooth(aes(x=week, y = tot_contrib), color = "black")+
  geom_point(aes(x=week, y = var), color = "red")+
  geom_smooth(aes(x=week, y = var), color = "red")+
  geom_point(aes(x=week, y = mean_contrib), color = "blue") +
  geom_smooth(aes(x=week, y = mean_contrib), color = "blue")
  


  
ggplot(data = move_sum_ind)+
  geom_point(aes(x=mean_lat, y = tot_contrib))
