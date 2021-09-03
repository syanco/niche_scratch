#!/usr/bin/env Rscript  
# chmod 744 script_template.r #Use to make executable

# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

'
Template

Usage:
script_template <dat> <out> [--seed=<seed>] [-b] [-t]
script_template (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-b --rollback   Rollback transaction if set to true.
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  .seed <- NULL
  .rollback <- TRUE
  .test <- TRUE
  rd <- here::here
  
  .datPF <- file.path(.wd,'analysis/dat_anno.csv')
  #TODO: figure out what the output will be
  # .outPF <- file.path(.wd,'analysis/track_annotated.csv')
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .rollback <- as.logical(ag$rollback)
  .test <- as.logical(ag$test)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/funs/input_parse.r'))
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
}

#---- Initialize Environment ----#
.seed <- ifelse(is.null(.seed),5326,as.numeric(.seed))

set.seed(.seed)
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(MVNH)
    # library(plyr)
    library(lubridate)
    library(dplyr)
    library(amt)
    library(grid)
    library(gridExtra)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Local parameters ----#

#---- Load control files ----#

#---- Load data ----#
message('Loading data...')
options(scipen=999)
dat0 <- read_csv("data/Migration timing in barnacle geese (Barents Sea) (data from Kölzsch et al. and Shariatinajafabadi et al. 2014)-5326125163732046073.csv",
                 col_types = "???????????????nnnnnnn"
)

#load data
dat <- dat0 %>% 
  #rename columns to easier names to use...
  rename(ind=`individual-local-identifier`,
         temp = `ECMWF Interim Full Daily SFC Temperature (2 m above Ground)`,
         gpp = `MODIS Land GPP 500m 8d Terra GPP`,
         lai = `MODIS Land LAI & FPAR 500m 4d Combined LAI`,
         fpar = `MODIS Land LAI & FPAR 500m 4d Combined FPAR`,
         ndvi= `MODIS Land Vegetation Indices 1km 16d Terra NDVI`,
         # `GPM Precipitation V05B Precipitation (Calibrated)`,
         elev = `height-above-msl`)%>% 
  #remove observations with high groundspeed
  filter(`ground-speed`<10) %>% 
  #create a few new vars:
  dplyr::mutate(week = week(timestamp), #week var
                year = year(timestamp), #year var
                WKxYR = glue("{year}{week}")) #year X week combo var 

#extract just the niche variables and a few key identifiers
dat_niche <- dat %>% 
  select(ind, timestamp, WKxYR, temp, gpp, lai, fpar, ndvi, elev) %>%
  drop_na()

#create `amt` formatted movement track for calculating NSD
dat_track <- dat %>% 
  # filter(individual_id == i) %>% 
  mutate(ts = ymd_hms(timestamp)) %>% #reformate timestamp as lubridate
  arrange(ts) %>% #order by ts
  # convert to `amt` track
  make_track(.x = `location-long`, .y = `location-lat`, .t = ts, order_by_ts = F, 
             crs = CRS("+init=epsg:4326"), all_cols = T)

#calculate NSD for each individual
dat_track_nsd <- dat_track %>% 
  group_by(ind) %>% #group by individual
  group_split() %>% #convert to a list of DFs by individual
  lapply(., FUN = function(x){ #apply across the leist
    x %>%   mutate(netSQ=nsd(.),
                   vel = netSQ-dplyr::lag(netSQ)) #add NSD variable to each individual DF
  }) %>% 
  do.call("rbind", .) #bind the DFs back together into a single DF

# tot_det <- dat_niche %>% 
#   select(temp,
#          gpp
#          # lai,
#          # fpar,
#          # ndvi,
#          # `GPM Precipitation V05B Precipitation (Calibrated)`,
#          # elev
#   )  %>% 
#   MVNH_det(data = ., log = T)
# 
# sum(tot_det[2:4]) - as.numeric(tot_det[1])
# # tot_det
# 
# ind_niche <- split(dat_niche, dat_niche$ind)
# 
# ind_det_l <- lapply(ind_niche, FUN = function(x){
#   x %>% select(temp, gpp) %>% 
#     MVNH_det(data = ., log = T)
# })
# 
# ind_det <- do.call("rbind", ind_det_l) %>% 
#   as.data.frame() %>% 
#   rownames_to_column

# expand_grid(ind_det$rowname, ind_det$rowname, .name_repair = "universal")

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
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(ind==ind_ID) %>% #select target individual
    group_by_(interval) %>% #group by the supplied time interval
    mutate(n=n()) %>% #calculate # of observations per interval
    filter(n >= min_obs) %>% #filter out those below the supplied minimum
    ungroup() #ungroup for next operation
  
  #calculate niches
  indMVNH <- dat_ind %>%
    group_by_(interval) %>% #group by time interval
    group_split() %>% #split into lists of DFs by interval
    lapply(., FUN = function(x){ #apply across elements of ths list
      x %>% select(!!!.vars) %>% #select only the supplied vars
        MVNH_det(data = ., log = log) #calculate niche hypevolume
    }) %>%
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  #create matching identifier information to link back up to the niche estimates
  indID <- dat_ind %>% 
    group_by_(interval) %>% #group by time interval
    group_split() %>% # split into lists of DFs by interval
    lapply(., FUN = function(x){#apply across elements of ths list
      #summarize the identifier info (just use first record)
      x %>% summarize(WKxYR = WKxYR[1], # week X yr factor
                      ind = ind[1], #ind ID
                      ts = timestamp[1]) #interval starting timestamp
    }) %>%
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  out <- cbind(indMVNH, indID) #bind the niche estimates back to the identifiers
  return(out) #return the combined DF for the target indivudal
}

#apply the above fxn across all individuals and bind into a df for all 
ind_ts_dat <- lapply(as.list(unique(dat_niche$ind)), FUN = indNTS, data = dat_niche, interval = "WKxYR",
                     vars=c(temp, gpp), log = T) %>%
  do.call("rbind", .)

# TODO:  I think this is unnecesary, remove?
# # join
# dat_niche_out <- dat_niche %>%
#   group_by(ind, WKxYR) %>% 
#   mutate(n=n()) %>%
#   filter(n >= 2) %>%
#   summarise() %>% 
#   left_join(ind_ts_dat, by = c("ind"="ind", "WKxYR"="WKxYR")) %>% 
#   ungroup()
# 
# #summarize nsd by week
# dat_comb <- dat_track_nsd %>% 
#   group_by(ind, WKxYR) %>% 
#   summarize(m_netSQ = mean(netSQ),
#             ts = t_[1],
#             ind = ind[1]) %>% 
#   right_join(dat_niche_out)


# conbine niche DF with NSD DF
dat_comb <- dat_track_nsd %>% #start with NSD
  group_by(ind, WKxYR) %>% #group by individual and weeek
  summarize(m_netSQ = mean(netSQ), #calc NSD as the weekly mean
            m_vel = mean(vel),
            ts = t_[1], #grab the first timestamp of the interval
            ind = ind[1]) %>% #grab individual ID
  right_join(ind_ts_dat) #join to the niche size dataframe



#rolling dissimilarity

#fx breaks data into a list of lists first order is individual, 2nd order is interval
# extracts specified variables
breakDat <- function(ind_ID, data, interval, min_obs = 2, vars){
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  
  data %>% 
    filter(ind==ind_ID) %>% #extract target individual
    group_by_(interval) %>% #group by the time interval
    mutate(n=n()) %>% # count number of obs in each interval
    filter(n >= min_obs) %>% #filter those below minimum threshold
    group_split() %>% #split into a list of DFs by grouping
    # select(temp, gpp) %>%
    lapply(., FUN = function(x){ #apply across list
      x %>% select(!!!.vars) # extract the eselected variables
    })
  #function returns the list of intervals for the target individual
}

#apply breakDat across individuls gathering only the niche variables
dat_break <- lapply(as.list(unique(dat_niche$ind)), FUN = breakDat, data = dat_niche, interval = "WKxYR",
                    vars = c(temp, gpp, elev))
#same as above but gather time and identifier information to match niches
break_key <- lapply(as.list(unique(dat_niche$ind)), FUN = breakDat, data = dat_niche, interval = "WKxYR",
                    vars = c(ind, timestamp, WKxYR))

# fxn to get dissimilarity between niche at each timestep andniche at first
#TODO: I think this doesn't work b/c of how `MVNH_dissimilarity` formats output
rollingDissim <- function(data){
  tryCatch(lapply(data, FUN = function(x){
    # .vars <- enquos(vars)
    diss_ts <- list()
    # x <- lapply(data, FUN = function(x, .vars){
    #   x %>% 
    #     select(!!!.vars) 
    # },.vars=.vars)
    niche0 <- x[[1]]
    for(i in 2:length(x)){
      diss <- MVNH_dissimilarity(niche0, x[[i]])
      BD <- data.frame(as.list(diss$Bhattacharyya_distance), metric = "BD")
      MD <- data.frame(as.list(diss$Mahalanobis_distance), metric = "MD")
      DR <- data.frame(as.list(diss$Determinant_ratio), metric = "DR")
      diss_ts[[i-1]] <- rbind(BD, MD, DR)
    }
    return(diss_ts)
  }
  # , vars = c(temp, gpp, elev)
  ))
  return(diss_ts)
}  

ind_dissim <- rollingDissim(data=dat_break[[1]])


diss_ts <- list()
# x <- lapply(data, FUN = function(x, .vars){
#   x %>% 
#     select(!!!.vars) 
# },.vars=.vars)
niche0 <- dat_break[[1]][[1]]
for(i in 2:length(dat_break[[1]])){
  diss <- MVNH_dissimilarity(niche0, dat_break[[1]][[i]])
  BD <- data.frame(as.list(diss$Bhattacharyya_distance), metric = "BD")
  MD <- data.frame(as.list(diss$Mahalanobis_distance), metric = "MD")
  DR <- data.frame(as.list(diss$Determinant_ratio), metric = "DR")
  diss_ts[[i-1]] <- rbind(BD, MD, DR)
}

key_ts <- list()
for(i in 2:length(break_key[[1]])){
  key_ts[[i-1]] <- break_key[[1]][[i]] %>% 
    summarize(ind = rep(ind[1], 3),
              timestamp = rep(timestamp[1], 3),
              WKxYR = rep(WKxYR[1], 3))
}


ts_comb_one <- list()
for(j in 1:length(key_ts)){
  ts_comb_one[[j]] <- cbind(diss_ts[[j]], key_ts[[j]])
}
ts_comb_one_df <- bind_rows(ts_comb_one)


dis_70833 <- ggplot(ts_comb_one_df)+
  geom_line(aes(x=timestamp, y = total))+
  ylab("Niche Dissimilarity")+
  theme(legend.position = "none")
(nsd_70833 <- ggplot(dat_comb[dat_comb$ind == "78033" & !is.na(dat_comb$m_netSQ),])+
    geom_line(aes(x=ts, y=m_netSQ, color = ind))+
    ylab("Net Sq. Displacement")+
    theme(legend.position = "none"))
# (nsd_70833 <- ggplot()+
#     geom_point(data=dat_track_nsd[dat_track_nsd$ind == "78033",], aes(x=t_, y=netSQ, color = ind))+
#     ylab("Net Sq. Displacement")+
#     theme(legend.position = "none"))

size_70833 <- ggplot()+
  geom_line(data=dat_niche_out[dat_niche_out$ind == "78033",], aes(x=ts, y=total, color = ind))+
  ylab("Niche Size")+
  theme(legend.position = "none")
(size_70833 <- ggplot(dat_comb[dat_comb$ind == "78033",])+
    geom_line(aes(x=ts, y=total))+
    geom_line(aes(x=ts, y=gpp), color = "darkgreen")+
    geom_line(aes(x=ts, y=temp), color = "red")+
    geom_line(aes(x=ts, y=cor), color = "blue")+
    ylab("Niche Size")+
    theme(
      # legend.position = "none"
      ))

gl <- list(nsd_70833, dis_70833, size_70833)

dissim_fig <- arrangeGrob(
  grobs = gl,
  heights = unit(c(4,4,4), c("in", "in", "in")),
  #  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1, 1),
                        c(2, 2),
                        c(3, 3))
)

ggsave("analysis/scenario1/figs/comb_70833.pdf", dissim_fig, height = 12)


ggplot(data=dat_comb) +
  geom_point(aes(x=abs(m_vel), y = cor))

# +
#   geom_smooth()

dissim_key <- lapply(break_key, FUN = function(x){
  key_ts <- list()
  for(i in 2:length(x)){
    key_ts[[i-1]] <- x[[i]] %>% 
      summarize(ind = rep(ind[1], 3),
                timestamp = rep(timestamp[1], 3),
                WKxYR = rep(WKxYR[1], 3))
  }
  return(key_ts)
})

dissim_tot <- list()
for(i in 1:length(dissim_key)){
  tmp <- list()
  for(j in 1:length(dissim_key[[i]])){
    tmp[[j]] <- cbind(ind_dissim[[i]][[j]], dissim_key[[i]][[j]])
  }
  dissim_tot[[i]] <- tmp
}
dissim_tot_df <- bind_rows(dissim_tot)



# PLOTS
ggplot()+
  geom_line(data=dat_niche_out, aes(x=ts, y=total, color = ind))+
  facet_wrap(~ind)


ggplot()+
  geom_line(data=dat_track_nsd, aes(x=t_, y=netSQ, color = ind))+
  facet_wrap(~ind)

ggplot(dat_comb)+
  geom_point(aes(x=m_netSQ, y = total, color = ind))

dissim_tot_df %>% 
  filter(metric == "BD") %>% 
  ggplot()+
  geom_line(aes(x=timestamp, y=total, color = ind))
# testFilter <- function(ind_ID, data){
#   out <- data %>%
#     filter(ind==ind_ID)
#   print(ind_ID)
#   return(out)
# }

# lapply(as.list(unique(dat_niche$ind)), FUN = testFilter, data = dat_niche)
########################################################

dat_ind_ts <- dat_niche %>% 
  filter(ind=="78033") %>% 
  group_by(interval) %>% 
  mutate(n=n()) %>% 
  filter(n >= 2) %>% 
  ungroup()


det_78033 <- dat_ind_ts %>% 
  group_split("week") %>% 
  # select(temp, gpp) %>% 
  lapply(., FUN = function(x){
    x %>% select(all_of(vars)) %>% 
      MVNH_det(data = ., log = F)
  })

det_78033_df <- do.call("rbind", det_78033) %>% 
  as.data.frame() %>% 
  rownames_to_column



plot(x=det_78033_df$rowname, y = det_78033_df$total, type = "l")

MVNH_dissimilarity(db1=ind_niche[[4]][,2:3], db2=ind_niche[[2]][,2:3])



chulls <- dat_niche %>% 
  group_by(ind) %>% 
  do(.[chull(.$temp, .$gpp), ])

ggplot()+
  # geom_point(data = dat_niche[,1:3], aes(x=temp, y = gpp, color = ind))+
  geom_polygon(data=chulls, aes(x=temp, y = gpp, color = ind), fill = NA)


#====

#---- Perform analysis ----#


#---- Save output ---#
message(glue('Saving to {.outPF}'))

dir.create(dirname(.outPF),recursive=TRUE,showWarnings=FALSE)

datout %>% write_csv(.outPF)

#---- Finalize script ----#

if(!.test) {
  library(uuid)
  
  .parPF <- file.path(.wd,"run_params.csv")
  
  #Update repo and pull out commit sha
  saveParams(.parPF)
}

message(glue('Script complete in {diffmin(t0)} minutes'))