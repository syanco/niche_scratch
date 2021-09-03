#Script to simulate a species niche to test if inclusion of irrelevant variables creates bias in niche size estimation

library(MVNH)
library(tidyverse)
options(scipen = 999)
#simulate environmental variables

#parameters for the variables
#Var 1
mu_v1 <- 0
sd_v1 <- 100

#Var 2
mu_v2 <- 0
sd_v2 <- 100


grid_n <- 100000 # size fo the #world

#create the environment - assume independence between variables
env1 <- rnorm(grid_n, mu_v1, sd_v1)
env2 <- rnorm(grid_n, mu_v2, sd_v2)
# env3 <- rnorm(grid_n, mu_v3, sd_v3)


#set up RSF model
#define coefficients
b0 <- -1  #intercept, doesn't really matter
b1 <- 3
b2 <- 2

#get untransformed model
z <- b0 + b1*env1 + b2*env2 

#transform to logit
pr <- 1/(1+exp(-z))         # pass through an inv-logit function
summary(pr)
hist(pr)


nicheIDX <- sample(1:grid_n, 100, prob = pr, replace = T)

dat1 <- data.frame(var1=env1[nicheIDX],
                   var2=env2[nicheIDX]) %>% 
  as.matrix()


#Do it again but now change beta for one var to 0..
#set up RSF model
#define coefficients
b0.2 <- -1  #interceot, doesn't really matter
b1.2 <- 0
b2.2 <- 2

#get untransformed model
z.2 <- b0.2 + b1.2*env1 + b2.2*env2 

#transform to logit
pr.2 <- 1/(1+exp(-z.2))         # pass through an inv-logit function
summary(pr.2)
hist(pr.2)

nicheIDX.2 <- sample(1:grid_n, 100, prob = pr.2, replace = T)

dat2 <- data.frame(var1=env1[nicheIDX.2],
                   var2=env2[nicheIDX.2]) %>% 
  as.matrix()

dat3 <- data.frame(var2=env2[nicheIDX.2]) %>% 
  as.matrix()

MVNH_det(dat1, log = T)
MVNH_det(dat2, log = T)
MVNH_det(dat3, log = T)

MVNH_det(dat1, log = F)
MVNH_det(dat2, log = F)
MVNH_det(dat3, log = F)


par(mfrow=c(3,1))
hist(env1, xlim= c(-450,450), main = "Env Total")
hist(env1[nicheIDX], xlim= c(-450,450), main = "Selected, beta = 3")
hist(env1[nicheIDX.2], xlim= c(-450,450), main = "Selected, beta = 0")

par(mfrow=c(3,1))
hist(env2, xlim= c(-450,450), main = "Env Total")
hist(env2[nicheIDX], xlim= c(-450,450), main = "Selected, beta = 2")
hist(env2[nicheIDX.2], xlim= c(-450,450), main = "Selected, beta = 2")
