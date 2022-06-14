# Test using a jags-based hmm to parse dbbmm var into states
library(jagsUI)
library(move)

load("data/542803025.rdata")

sink(file = "src/jags/jags_dbbmmvar_hmm.txt")

cat("
    model{
    
    #PRIORS
    for(k in 1:n.states){ #gamma priors for each of k states
      
      #Make priors for gamma parameters as function of moments
      a[k] <- pow(m[k],2) / pow(sd[k],2) #prior for variance shape param
      b[k] <- m[k] / pow(sd[k],2) #pior for variance rate param
      
      m[k] ~ dunif(.01,100) #pror for mean variance
      sd[k] ~ dunif(0.01,100) #prior for SD of variance
    }#k
    
    #treat time-varying transition probabilities as homogenous
    for(t in 1:steps){
      g11[t] <- mean.g11
      g12[t] <- mean.g12
      g21[t] <- mean.g21
      g22[t] <- mean.g22
    } #t
    
    #priors for transition probabilities
    mean.g11 ~ dbeta(1,1)
    mean.g12 ~ dbeta(1,1)
    mean.g21 ~ dbeta(1,1)
    mean.g22 ~ dbeta(1,1)
    
    #State Transition Matrix
    for(t in 1:steps){
      ps[1,1,t] <- g11[t]
      ps[1,2,t] <- g12[t]
      ps[2,1,t] <- g21[t]
      ps[2,2,t] <- g22[t]
    } #t
    
    #LIKELIHOOD
    
    for(i in 1:n.ind) {
      idx[i,1] <- 1
      
      for(t in 2:steps){ #loop over time steps
        
        #Process Model
        idx[i,t] ~ dcat(ps[idx[i,t-1], ,t-1])
        
        #Observation Model
        v[i,t] ~ dgamma(a[idx[i,t]], b[idx[i,t]]) # likelihood for steps
        
      }#t
    } #i

    }", fill = T)
sink()


# ==== Run Models ==== #

# #- Declare model input
# 
# mean.phiW <- if(ctf$wPar==0){rep(0.99, 2)}else{NA}
# mean.phiM <- if(ctf$mPar==0){rep(0.99, 2)}else{NA}
# mean.phiB <- if(ctf$bPar==0){matrix(0.99, nrow = 3, ncol = 2)}else{NA}
# muW <- if(ctf$wPar==0){NA}else{rep(0, 2)}
# betaW <- if(ctf$wPar==0){NA}else{rep(0, 2)}
# muM <- if(ctf$mPar==0){NA}else{rep(0, 2)}
# betaM <- if(ctf$mPar==0){NA}else{rep(0, 2)}
# muB <- if(ctf$bPar==0){NA}else{matrix(0, nrow = 3, ncol = 2)}
# betaB <- if(ctf$bPar==0){NA}else{matrix(0, nrow = 3, ncol = 2)}
# 
# # Initial values
# inits <- function(){
#   l <- list(
#     mean.phiW = mean.phiW,
#     mean.phiM = mean.phiM,
#     mean.phiB = mean.phiB,
#     mean.pW = 0.25,
#     mean.pM = 0.25,
#     mean.pB = 0.5,
#     z = z.init,
#     muW = muW,
#     betaW = betaW,
#     muM = muM,
#     betaM = betaM,
#     muB = muB,
#     betaB = betaB
#   )
#   out <- l[!is.na(l)]
#   return(out)
# }
# 
# var <- if(is.na(ctf$var)){NA}else{get(ctf$var)}

# y <- as.numeric(na.omit(tmp_out$`dBBMM Variance`@means))
# 
# v <- matrix(y, nrow = 6)

v <- dat1[,2:ncol(dat1)]

## Bundle data
jags.data <- list(v = v,
                  n.ind = nrow(v),
                  steps = ncol(v),
                  n.states = 2)


### Parameters
params <- c("a", "b", "m", "sd", "idx")

### MCMC settings
nc <- 4
ni <- 5000
nb <- 2000
# nt <- 500
# na <- 1000

# Fit model and view output
mod <- jagsUI::jags(data = jags.data, 
                    # inits = inits, 
                    parameters.to.save = params, 
                    model.file = "src/jags/jags_dbbmmvar_hmm.txt", 
                    n.chains = nc, 
                    n.iter = ni, 
                    n.burnin = nb,
                    # n.thin = nt, 
                    # n.adapt = na, 
                    parallel = TRUE)

# save(mod, file = glue("{.outPF}/model_{ctf$model_no}.rdata"))

mod
plot(mod)
