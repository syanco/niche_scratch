#!/usr/bin/env Rscript --vanilla
# chmod 744 niche_contrib_by_seg.r #Use to make executable

# TODO: The whole top is wrong - description, docopts, etc

# ==== Setup ====

'
Outputs a dataframe of individual- season-specific niche means and variances.
Rows represent a single season for a single individual.

Usage:
ind_season_moments.r <db> <out> <ctf> 
ind_season_moments.r (-h | --help)

Parameters:
  db: path to movement databse. 
  out: path to output directory.
  ctf: path to segementation control file.
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  .rollback <- TRUE
  rd <- here::here
  
  .outPF <- file.path(.wd,'analysis/cranes/ind_to_pop/seasonal_seg')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfs <- file.path(.wd, "ctfs/crane_segmentation")
  .anno <- "anno_join_2021-08-16"
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  # .seed <- ag$seed
  # .rollback <- as.logical(ag$rollback)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .ctfs <- makePath(ag$ctf)
  .anno <- "anno_join_2021-08-16"
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(file.path(.wd,'src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    # library(raster)
    library(bayesmove)
  }))

#Source all files in the auto load funs directory
list.files(rd('funs/auto'),full.names=TRUE) %>%
  walk(source)

segs <- list.files(.ctfs, pattern = "*.csv")


#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))


#---- Perform analysis ----#
message("Gathering movement data...")

inds <- tbl(db, "individual") %>% 
  filter(taxon_canonical_name == "Anthropoides virgo") %>% 
  collect() %>% 
  pull(individual_id) 

# Grab one ind, sort by time within ind
evt0 <- tbl(db, .anno) %>% 
  filter(individual_id %in% inds) %>% 
  collect() %>% 
  mutate(date = strptime(timestamp, format = "%Y-%m-%d %H:%M:%OS"),
         #some functions in package don't allow you to tell it the individual identifier
         id = individual_id) %>%   
  group_by(individual_id) %>% 
  arrange(date)
  
####==== bayesmove ====####
#calc intervals, step length, turn angle
tracks<- prep_data(dat = evt0, coord.names = c("lon","lat"), id = "id")

#round the times to interval, use 5 minute tolerance
##TODO: need better way to find track interval...
tracks<- round_track_time(dat = tracks, id = "id", int = 1200, tol = 300, 
                          time.zone = "UTC", units = "secs") 

indsum <- tracks %>% 
  filter(dt == 1200) %>% 
  group_by(id) %>% 
  summarize(min = min(date),
            max = max(date),
            dur = as.period(max-min),
            n = n()) %>% 
  filter(dur > 100) %>% 
  filter(n > 200)

#filter tracks by duration (set in chunck above)
tracks <- tracks %>% 
  filter(id %in% indsum$id)

tracks.list<- df_to_list(dat = tracks, ind = "id")

# #remove inds with fewer than 100 fixes
# # TODO: revisit this filter threshold
# tracks.list <- tracks.list[lapply(tracks.list, FUN = function(x){
#   (min(x$date, na.rm = T)-max(x$date, na.rm = T))/86400
# # print(x$date)}
# # > 100
# })
# ]

# Filter observations
tracks_filt.list<- filter_time(dat.list = tracks.list, int = 1200)

#TODO:  these are the defaults from the vingette...
# Define bin number and limits for turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Define bin number and limits for step lengths
# dist.bin.lims=quantile(tracks[tracks$dt == 1200,]$step,
#                        c(0,0.05,0.1,0.15,0.2,0.25,0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
#                          0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1), na.rm=T)
dist.bin.lims=quantile(tracks[tracks$dt == 1200,]$step,
                       c(0,0.99, 1), na.rm=T)

# dist.bin.lims=seq(0, max(tracks[tracks$dt == 1200,]$step, na.rm = T), length.out = 6)

# Assign bins to observations
tracks_disc.list<- map(tracks_filt.list,
                       discrete_move_var,
                       lims = list(dist.bin.lims, angle.bin.lims),
                       varIn = c("step", "angle"),
                       varOut = c("SL", "TA"))

#select only vars needed (I guess?)
tracks.list2<- map(tracks_disc.list,
                   subset,
                   select = c(id, SL, TA))

# Segment...
# Define hyperparameter for prior distribution
alpha<- 1

# Set number of iterations for the Gibbs sampler
ngibbs<- 10000

# Set the number of bins used to discretize each data stream
nbins<- c(2,8)

progressr::handlers(progressr::handler_progress(clear = FALSE))
future::plan("multicore", workers = 2)  #run all MCMC chains in parallel
#refer to future::plan() for more details
dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)

future::plan(future::sequential)  #return to single core


traceplot(data = dat.res, type = "nbrks")
traceplot(data = dat.res, type = "LML")


MAP.est<- get_MAP(dat = dat.res$LML, nburn = 5000)
MAP.est
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))

plot_breakpoints(data = tracks.list, as_date = T, var_names = c("step","angle"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)"), brkpts = brkpts)

tracks.seg<- assign_tseg(dat = tracks.list, brkpts = brkpts)

####==== segclust2d ====####

library(segclust2d)
library(bayesmove)

#calc intervals, step length, turn angle
tracks<- prep_data(dat = evt0, coord.names = c("lon","lat"), id = "id")

evt0 %>% 
  group_by(individual_id) %>% 
  summarise(n=n())

evt1 <- evt0 %>% 
  filter(individual_id == inds[5])

x <- segclust(as.data.frame(evt1),
              # Kmax = 10, 
              lmin = 10,
              ncluster = c(2,3,4),
              # order.var = "date",
              seg.var = c("lon", "lat"))
