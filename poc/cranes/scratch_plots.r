library(gtools)

ind_ts_dat <- read_csv("analysis/cranes/niche_accum.csv") %>% 
  group_by(ind) %>% 
  arrange(ts) %>% 
  mutate(ts_imp = 1:n()) %>% 
  ungroup()

ind_ts_dat_rev <- ind_ts_dat %>% 
  group_by(ind) %>% 
  mutate(ts_imp = 1:n()) %>% 
  ungroup() %>% 
  filter(!is.na(c_range)) %>% 
  filter(!is.infinite(c_range),
         !is.infinite(min),
         !is.infinite(max))


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


nsd <- read_csv("analysis/cranes/nsd.csv")
ind_ID <- c(55754621,234543899,146932155,10836626,234544351, 195375963,963160836)
nsd_mod <- nsd %>% 
  filter(individual_id %in% ind_ID)


for(i in 1:length(ind_ID)){

    ggplot(ind_ts_dat, aes(x=ts, color = as.factor(ind))) +
    geom_linerange(aes(ymin=min, ymax = max)) +
    theme_bw() +
    ylim(c(-1,1)) +
    facet_wrap(~as.factor(ind), scale = "free")
  
  ggplot(nsd_mod)+
    geom_line(aes(x=timestamp, y = netSQ))+
    ggtitle("Net Sq. Displacement") +
    theme_bw() +
    facet_wrap(~as.factor(individual_id), scale = "free")
  
}

##############################################

ggplot(ind_ts_dat_rev, aes(x=ts_imp, y = c_dist))+
  geom_point(aes(color = as.factor(ind))) +
  geom_smooth(method = "lm")

ind_ts_dat_rev <- ind_ts_dat_rev %>% 
  mutate(burnin = case_when(ts_imp < 30 ~ "1",
                            ts_imp > 30 ~ "0"))

ggplot(ind_ts_dat_rev , aes(x = ts_imp, y = c_range,
                            color=as.factor(burnin)))+
  geom_point(
    # aes(color=as.factor(ind))
  )+
  geom_smooth(method = "lm")


# Logistic model
#get inits
phi1 <- 3 # guess for asymptote
init_mod <- coef(lm(logit(c_range/phi1)~ts_imp,data=ind_ts_dat_rev))

mod_time <- nls(c_range ~ phi1/(1+exp(-(phi2+phi3*ts_imp))),
                start = list(phi1 = phi1, phi2 = init_mod[1], phi3=init_mod[2]),
                data=ind_ts_dat_rev,
                trace=TRUE)

summary(mod_time)

newx <- c(min(ind_ts_dat_rev$ts_imp):max(ind_ts_dat_rev$ts_imp)) #construct a range of x values bounded by the data
newy <- coef(mod_time)[1]/(1+exp(-(coef(mod_time)[2]+coef(mod_time)[3]*newx))) #predicted mass
predict <- data.frame(newx,newy) #create the prediction data frame#And add a nice plot (I cheated and added the awesome inset jpg in another program)

ggplot(data=ind_ts_dat_rev, aes(x=ts_imp, y=c_range))+
  geom_point(aes(color=as.factor(ind)))+
  theme_bw()+
  labs(x='Days', y='Cummulative EVI Range')+
  # scale_x_continuous(breaks=c(0,250,500,750, 1000,1250))+
  # scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
  # theme(axis.text=element_text(size=18),axis.title=element_text(size=24))+
  geom_line(data = predict, aes(x = newx, y = newy), size=1)

#logarithmic
init_ln <- coef(lm(log(c_range)~ts_imp,data=ind_ts_dat_rev))
mod_ln <- nls(c_range ~ A*log(ts_imp)+B,
              start = list(A = init_ln[2], B = 0),
              data=ind_ts_dat_rev,
              trace=TRUE)
summary(mod_ln)

newy <- logarithmicGrowth(t=newx, A = coef(mod_ln)["A"], B = coef(mod_ln)["B"])
predict <- data.frame(newx,newy) #create the prediction data frame#And add a nice plot (I cheated and added the awesome inset jpg in another program)

ggplot(data=ind_ts_dat_rev, aes(x=ts_imp, y=c_range))+
  geom_point(aes(color=as.factor(ind)))+
  theme_bw()+
  labs(x='Days', y='Cummulative EVI Range')+
  # scale_x_continuous(breaks=c(0,250,500,750, 1000,1250))+
  # scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
  # theme(axis.text=element_text(size=18),axis.title=element_text(size=24))+
  geom_line(data = predict, aes(x = newx, y = newy), size=1)


######################################3
logarithmicGrowth <- function(t, A, B) {
  n <- A*log(t)+B
  return(n)
}
x <- 1:300
A <- 1
B <- 0
y<- logarithmicGrowth(x, A=A, B=B)
plot(x, y, ylim= c(0,10))

#####################################
# Fit logarithmic growth to an individual, extract and plot residuals
newx <- c(min(ind_ts_dat$ts_imp):max(ind_ts_dat$ts_imp)) #construct a range of x values bounded by the data


#fit model
ind_mods <- list()

dat_ls <- ind_ts_dat %>% 
  group_by(ind) %>% 
  group_split()

for(i in 1:length(unique(ind_ts_dat_rev$ind))){
  # get inits from a simple lm
  init_ln <- coef(lm(log(c_range)~ts_imp, data = dat_ls[[i]]))
  
  
  ind_mods[[i]] <- nls(c_range ~ A*log(ts_imp)+B,
                       start = list(A = init_ln[2], B = 0),
                       data=dat_ls[[i]],
                       trace=TRUE)  
  
  readline(prompt="Press [enter] to continue")
  newy <- logarithmicGrowth(t=newx, A = coef(ind_mods[[i]])["A"], 
                            B = coef(ind_mods[[i]])["B"])
  predict <- data.frame(newx,newy) #create the prediction data frame#And add a nice plot (I cheated and added the awesome inset jpg in another program)
  
  print(ggplot(data=dat_ls[[i]], aes(x=ts_imp, y=c_range))+
          geom_point(aes(color=as.factor(ind)))+
          theme_bw()+
          labs(x='Days', y='Cummulative EVI Range')+
          ggtitle(glue("Crane {dat_ls[[i]]$ind}")) +
          # scale_x_continuous(breaks=c(0,250,500,750, 1000,1250))+
          # scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80))+
          # theme(axis.text=element_text(size=18),axis.title=element_text(size=24))+
          geom_line(data = predict, aes(x = newx, y = newy), size=1))
  
  readline(prompt="Press [enter] to continue")
  
  resid_df <- data.frame(x = dat_ls[[i]]$ts_imp,
                         y = resid(ind_mods[[i]]))
  
  print(ggplot(resid_df) +
          geom_point(aes(x = x, y = y)) +
          geom_hline(yintercept = 0) +
          ggtitle(glue("Crane {dat_ls[[i]]$ind}")))
}



