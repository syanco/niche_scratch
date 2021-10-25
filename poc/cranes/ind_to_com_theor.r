
# Calculate the poulation niche mean.
# `ind_mus`: vector of individual means
popMu <- function(ind_mus){
  n <- length(ind_mus)
  sum(ind_mus)/n
}

# Calculate the community niche mean.
# `ws`: vector of population weights - typically relative abundance
# `pop_mus`: vector of population means
comMu <- function(ws, pop_mus){
  stopifnot(length(ws) == length(pop_mus))
  vals <- ws*pop_mus
}

# Calculate the population niche variance
# `sd`: vector of individual standard deviations
# `means`: vector of individual means
popVar <- function(sd, ind_mus){
  pop_mean <- mean(ind_mus)
  n <- length(sd)
  
  vec <- c()
  
  for(i in 1:n){
    vec[i] <- (1/n)*(sd[i]^2 + ind_mus[i]^2 - pop_mean^2)
  }
  
  out <- sum(vec)
  return(out)
}

# Calculate the community niche variance
# `ws`: vector of population weights - typically relative abundance
# `pop_var`: vector of population niche variances
# `means`: vector of population niche means
comVar <- function(ws, pop_var, pop_mus){
  stopifnot(length(ws) == length(pop_mus))
  stopifnot(length(pop_var) == length(pop_mus))
  stopifnot(length(ws) == length(pop_var))
  
  com_mean <- comMu(ws=ws, pop_mus=pop_mus)
  n <- length(pop)
  
  vec <- c()
  
  for(i in 1:n){
    vec[i] <- w[i]*(pop_var[i] + pop_mus[i]^2 - com_mean^2)
  }
  
  out <- sum(vec)
  return(out)
}

normalDensity <- function(mu, sd){
  x <- seq(-100, 100, length=10000)
  
  #create a vector of values that shows the height of the probability distribution
  #for each value in x
  y <- dnorm(x=x, mean = mu, sd = sd)
  
  return(data.frame(x=x, y=y))
}


#---- SIMS ----#

mu_is <- rnorm(n=100, mean=0, sd=100)
var_is <- rgamma(shape= 9, 100)

pop_mu <- popMu(mu_is)
pop_var <- popVar(sd=sqrt(var_is), mu_is)

den <- normalDensity(pop_mu, sqrt(pop_var))
plot(den$x, den$y, type = "l")


#==============================================================================#

#---- Using ind means ----#


#-- Functions --#

comVar2 <- function(mean_pop_var, mean_pop_diff){
  mean_pop_var + mean_pop_diff
}

popVar2 <- fucntion(mean_ind_var, mean_ind_diff){
  mean_ind_var + mean_ind_diff
}

meanIndDiff <- function(mu_vec){
  #get population mean
  pop_mu <- mean(mu_vec)
  
  #get num individuals
  n <- length(mu_vec)
  
  #init empty vector to hold ind diffs
  v <- c()
  #calc ind diffs
  for(i in 1:n){
    (1/n)*(mu_vec[i]^2-pop-mu^2)
  }
  
  #take the sum of the weighted vector
  out <- sum(v)
  
  return(out)
}

meanPopDiff <- function(mu_vec, w_vec){
  #get population mean
  com_mu <- mean(mu_vec)
  
  #get num individuals
  n <- length(mu_vec)
  
  #init empty vector to hold ind diffs
  v <- c()
  #calc ind diffs
  for(i in 1:n){
    w_vec*(mu_vec[i]^2-com-mu^2)
  }
  
  #take the sum of the weighted vector
  out <- sum(v)
  
  return(out)
}


#-- Demo --#

# Individual trade-off within population
pop_var <- 10

mean_ind_var <- seq(0,10, by = 0.1)

mean_mu_diff <- c()
for(i in 1:length(mean_ind_var)){
  mean_mu_diff[i] <- pop_var-mean_ind_var[i]
}

plot(x=mean_ind_var, y=mean_mu_diff)
