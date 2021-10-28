
#create simple not in function (negation of `%in%`)
`%notin%` <- Negate(`%in%`)

# function to calculate individual-specific time-varying niche size
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNTS <- function(ind_ID, data, interval, min_obs = 2, vars, log=F){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    group_by_(interval) %>% #group by the supplied time interval
    mutate(n=n()) %>% #calculate # of observations per interval
    filter(n >= min_obs) %>% #filter out those below the supplied minimum
    ungroup() #ungroup for next operation
  
  #calculate niches
  indMVNH <- dat_ind %>%
    group_by_(interval) %>% #group by time interval
    group_split() %>% #split into lists of DFs by interval
    lapply(., FUN = function(x){ #apply across elements of ths list
      x %>% select(all_of(vars)) %>% #select only the supplied vars
        MVNH_det(data = ., log = log) #calculate niche hypevolume
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  #create matching identifier information to link back up to the niche estimates
  indID <- dat_ind %>% 
    group_by_(interval) %>% #group by time interval
    group_split() %>% # split into lists of DFs by interval
    lapply(., FUN = function(x){#apply across elements of the list
      #summarize the identifier info (just use first record)
      x %>% summarize(WKxYR = WKxYR[1], # week X yr factor
                      ind = individual_id[1], #ind ID
                      ts = timestamp[1]) #interval starting time stamp
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  out <- cbind(indMVNH, indID) #bind the niche estimates back to the identifiers
  return(out) #return the combined DF for the target individual
}


# function to calculate individual-specific time-varying niche dissimilarity
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNDissim <- function(ind_ID, data, interval, min_obs = 2, vars){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  # .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    group_by_(interval) %>% #group by the supplied time interval
    mutate(n=n()) %>% #calculate # of observations per interval
    filter(n >= min_obs) %>% #filter out those below the supplied minimum
    ungroup() #ungroup for next operation
  
  #get individual total niche to use in calculating dissim below.
  niche0 <- dat_ind %>% 
    select(all_of(vars)) %>% 
    na.omit()
  
  #calculate niches
  indMVNH <- dat_ind %>%
    group_by_(interval) %>% #group by time interval
    group_split() %>% #split into lists of DFs by interval
    lapply(., FUN = function(x){ #apply across elements of this list
      tryCatch(x %>% select(all_of(vars)) %>% #select only the supplied vars
                 na.omit() %>% # omit NAs
                 MVNH_dissimilarity(db1 = ., db2 = niche0), #calculate niche dissim
               error=function(e) {list(Bhattacharyya_distance = NA,
                                       Mahalanobis_distance = NA,
                                       Determinant_ratio = NA)}
      ) # tryCatch
    } #f(x)
    ) %>% # lapply
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  #create matching identifier information to link back up to the niche estimates
  indID <- dat_ind %>% 
    group_by_(interval) %>% #group by time interval
    group_split() %>% # split into lists of DFs by interval
    lapply(., FUN = function(x){#apply across elements of the list
      #summarize the identifier info (just use first record)
      x %>%  
        summarize(WKxYR = WKxYR[1], # week X yr factor
                  ind = individual_id[1], #ind ID
                  ts = timestamp[1] #interval starting time stamp
        ) 
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  tryCatch(
    out <- cbind(indMVNH, indID) #bind the niche estimates back to the identifiers
  ) # tryCatch
  
  return(out) #return the combined DF for the target individual
}


# function to calculate individual-specific time-varying niche size
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNicheAccum <- function(ind_ID, data, interval, min_obs = 2, vars){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    arrange({{interval}}) %>% 
    group_by({{interval}}) %>% #group by the supplied time interval
    # mutate(n=n()) %>% #calculate # of observations per interval
    # filter(n >= min_obs) %>% #filter out those below the supplied minimum
    summarise(n = n(),
              max = max({{vars}}, na.rm = T),
              min = min({{vars}}, na.rm = T),
              range = max-min) %>%
    ungroup()
  # 
  # nsd_dat <- dat_ind %>% 
  #   make_track(.x = lon, .y = lat, .t = ts, order_by_ts = F, 
  #            crs = CRS("+init=epsg:4326"), all_cols = T) %>% 
  #   mutate(netSQ=nsd(.),
  #          nsd_accum = netSQ+dplyr::lag(netSQ)) #add NSD variable to each individual DF
  # 
  
  return(dat_ind) #return the combined DF for the target individual
}


#-- Individual to Population Functions  --##

# Anticipates vector of sample means, vector of sample variances
estPopVar <- function(var, means, pop_mean){
  n <- length(var)
  
  vec <- c()
  
  for(i in 1:n){
    vec[i] <- (1/n)*(var[i] + means[i]^2 - pop_mean^2)
  }
  
  out <- sum(vec)
  return(out)
}

# calculate an individual's total contribution to pop variance 
indContrib <- function(mu_i, var_i, mu_pop, n){
  # n <- length(mu_i)
  cont <- (1/n)*(var_i + mu_i^2 - mu_pop^2)
  return(cont)
}

# calculate the mean component of ind contribution to pop variance
#(variance contrib is just the ind variance itself)
muContrib <- function(mu_i, mu_pop, n){
  cont <- (1/n)*(mu_i^2 - mu_pop^2)
  return(cont)
}
