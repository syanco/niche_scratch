# Scratch script to test out covariate dBBMM methods


#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  rd <- here::here
  
  # .outPF <- file.path(.wd,'analysis/cranes/nsd.csv')
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