#!/usr/bin/env Rscript --vanilla
# chmod 744 calc_niches.r #Use to make executable

# This script is used to calculate time-dynamic niche accumulation curves.

# ==== Setup ====

'
Calculate time-dynamic niches accumulation curves.

Usage:
niche_accumulation.r <db> <out> <ctf> <anno> 
niche_accumulation.r (-h | --help)

Parameters:
  dat: path to input csv file. 
  out: path to output directory.
  ctf: path to control file
  anno: name of annotation table in db

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
  
  .outPF <- file.path(.wd,'analysis/cranes/niche_sizes.csv')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfPF <- file.path(.wd, "analysis/cranes/ctfs/anno_vars.csv")
  .anno <- "anno_join_2021-09-08"
  
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
  .ctfPF <- makePath(ag$ctf)
  .anno <- as.character(ag$anno)
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
    library(MVNH)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Control Files ----#
ctf <- read_csv(.ctfPF)

#get run rows from ctf
runs <- ctf %>% 
  filter(run == 1)

.variables <- runs %>% 
  select(variable) 

.values <- .variables %>% 
  sapply(., FUN = function(x){glue("value_{x}")}) %>% 
  c()

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#
message("Loading annotations...")
anno0 <- tbl(db, .anno) %>% 
  collect()
message("Disconnecting from databse...")
dbDisconnect(db)

message("Gathering niche data...")

#load data
anno_mod <- anno0 %>% 
  #remove observations with high groundspeed
  #TODO: eventually remove this and make a separate script that cleans data
  # earlier in the workflow
  filter(`ground_speed`<10) %>% 
  #create a few new vars:
  dplyr::mutate(week = week(timestamp), #week var
                year = year(timestamp), #year var
                WKxYR = glue("{year}{week}")) #year X week combo var 

message("Calculating niche sizes...")
ind_ts_dat <- lapply(as.list(unique(anno_mod$individual_id)), FUN = indNTS, 
                     data = anno_mod, interval = "WKxYR",
                     vars=.values, log = T) %>%
  bind_rows()


#---- Finalize script ----#

message(glue("Saving output to: {.outPF}"))
write_csv(ind_ts_dat, path = .outPF)

message(glue('Script complete in {diffmin(t0)} minutes'))