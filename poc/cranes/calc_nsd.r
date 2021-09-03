#!/usr/bin/env Rscript --vanilla
# chmod 744 calc_nsd.r #Use to make executable

# This script calculates time-dynamic net squared displacement (NSD) across 
# individual animals.

# ==== Setup ====

'
Calculate time dynamic net squared displacement across multiple individuals

Usage:
join_annos.r <db> <out> 
join_annos.r (-h | --help)

Parameters:
  dat: path to input csv file. 
  out: path to output directory.
  
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
  
  .outPF <- file.path(.wd,'analysis/cranes/nsd.csv')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .rollback <- as.logical(ag$rollback)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(rd('src/startup.r'))
#TODO: fix this file path after moving the function script
source(rd("src/poc/cranes/niche_funs.r"))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(amt)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#
message("Loading annotations...")
evt0 <- tbl(db, "event") %>% 
  collect()
message("Disconnecting from databse...")
dbDisconnect(db)

message("Gathering movement data...")

#load data
evt_mod <- evt0 %>% 
  #remove observations with high groundspeed
  #TODO: eventually remove this and make a separate script that cleans data
  # earlier in the workflow
  filter(`ground_speed`<10) %>% 
  #create a few new vars:
  dplyr::mutate(week = week(timestamp), #week var
                year = year(timestamp), #year var
                WKxYR = glue("{year}{week}")) #year X week combo var 


# NOTE/TODO: keeping this intermediary object in code for now in case we want to spit 
# the `amt` track object out the side fo this script for other analyses
message("Making tracks...")

dat_track_nsd <- evt_mod %>% 
  # filter(individual_id == i) %>% 
  mutate(ts = ymd_hms(timestamp)) %>% #reformat timestamp as lubridate
  arrange(ts) %>% #order by ts
  group_by(individual_id) %>% #group by individual
  group_split() %>% #convert to a list of DFs by individual
  lapply(., FUN = function(x){ #apply across the leist
    x %>%   
      # convert to `amt` track
      make_track(.x = lon, .y = lat, .t = ts, order_by_ts = F, 
                 crs = CRS("+init=epsg:4326"), all_cols = T) %>% 
      mutate(netSQ=nsd(.),
             vel = netSQ-dplyr::lag(netSQ)) #add NSD variable to each individual DF
  }) %>% 
  do.call("rbind", .) #bind the DFs back together


#---- Finalize script ----#

message(glue("Saving output to: {.outPF}"))
write_csv(dat_track_nsd, path = .outPF)

message(glue('Script complete in {diffmin(t0)} minutes'))
