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

dat_niche <- dat0 %>% 
  rename(ind=`individual-local-identifier`,
         temp = `ECMWF Interim Full Daily SFC Temperature (2 m above Ground)`,
         gpp = `MODIS Land GPP 500m 8d Terra GPP`,
         lai = `MODIS Land LAI & FPAR 500m 4d Combined LAI`,
         fpar = `MODIS Land LAI & FPAR 500m 4d Combined FPAR`,
         ndvi= `MODIS Land Vegetation Indices 1km 16d Terra NDVI`,
         # `GPM Precipitation V05B Precipitation (Calibrated)`,
         elev = `height-above-msl`)%>% 
  filter(`ground-speed`>0) %>% 
  dplyr::mutate(week = week(timestamp)) %>% 
  select(ind, week, temp, gpp, lai, fpar, ndvi, elev) %>% 
  drop_na() 
# ind_mean <- dat_niche %>% 
#   group_by(individual_id) %>% 
#   summarise(mu_night = mean(na.omit(value_lst_night)),
#             mu_day = mean(na.omit(value_lst_day))) %>% 
#   select(mu_night, mu_day)
tot_det <- dat_niche %>% 
  select(temp,
         gpp
         # lai,
         # fpar,
         # ndvi,
         # `GPM Precipitation V05B Precipitation (Calibrated)`,
         # elev
  )  %>% 
  MVNH_det(data = ., log = T)

sum(tot_det[2:4]) - as.numeric(tot_det[1])
tot_det

ind_niche <- split(dat_niche, dat_niche$ind)

ind_det_l <- lapply(ind_niche, FUN = function(x){
  x %>% select(temp, gpp) %>% 
    MVNH_det(data = ., log = T)
})

ind_det <- do.call("rbind", ind_det_l) %>% 
  as.data.frame() %>% 
  rownames_to_column

# expand_grid(ind_det$rowname, ind_det$rowname, .name_repair = "universal")


indNTS <- function(ind_ID, data, interval, min_obs = 2, vars, log=F){
  dat_ind <- data %>% 
    filter(ind==ind_ID) %>%
    group_by_(interval) %>%
    mutate(n=n()) %>%
    filter(n >= min_obs) %>%
    ungroup()

  indMVNH <- dat_ind %>%
    group_by_(interval) %>%
    group_split() %>%
    # select(temp, gpp) %>%
    lapply(., FUN = function(x){
      x %>% select_(vars) %>%
        MVNH_det(data = ., log = log)
    }) %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    rownames_to_column %>%
    mutate(ind_ID)
  return(indMVNH)
}

ind_ts_dat <- lapply(as.list(unique(dat_niche$ind)), FUN = indNTS, data = dat_niche, interval = "week",
       vars=c("temp","gpp", "elev"), log = T) %>%
  do.call("rbind", .)



ggplot()+
  geom_line(data=ind_ts_dat, aes(x=as.numeric(rowname), y=total, color = ind_ID))+
  facet_wrap(~ind_ID)

testFilter <- function(ind_ID, data){
  out <- data %>%
    filter(ind==ind_ID)
  print(ind_ID)
  return(out)
}

lapply(as.list(unique(dat_niche$ind)), FUN = testFilter, data = dat_niche)
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