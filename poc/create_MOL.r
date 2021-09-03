#!/usr/bin/env Rscript
# chmod 744 script_template.r #Use to make executable

# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

'
Template

Usage:
script_template <out> [--seed=<seed>] [-t]
script_template (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-s --seed=<seed>  Random seed. Defaults to 5280 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '/home/syanco/projects/niche_scratch'
  .seed <- NULL
  .rollback <- TRUE
  .test <- TRUE
  rd <- here::here
  
  .outPF <- file.path(.wd, "analysis/dat_MOL.csv")
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
  
  .outPF <- makePath(ag$out)
}

#---- Initialize Environment ----#
.seed <- ifelse(is.null(.seed),5280,as.numeric(.seed))

set.seed(.seed)
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Local parameters ----#
.dbPF <- file.path(.wd,"analysis/movedb/data/move.db")

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))


#---- Load data ----#
message("Connected to database, processing data...")
indtb <- tbl(db,'individual')
evttb <- tbl(db,'event')

datout <- evttb %>% 
  select(individual_id, event_id, timestamp, lat, lon) %>% 
  inner_join(indtb, by = "individual_id") %>% 
  as.data.frame() %>% 
  mutate(date = date(ymd_hms(timestamp)))

#---- Save output ---#
message(glue('Saving to {.outPF}'))

datout %>% write_csv(.outPF)

#---- Finalize script ----#

if(!.test) {
  library(git2r)
  library(uuid)
  
  .runid <- UUIDgenerate()
  .parPF <- file.path(.wd,"run_params.csv")
  
  #Update repo and pull out commit sha
  repo <- repository(rd('src'))
  
  rstat <- status(repo)
  if(length(rstat$staged) + 
     length(rstat$unstaged) + 
     length(rstat$untracked) > 0) {
    add(repo,'.')
    commit(repo, glue('script auto update. runid: {.runid}'))
  }
  
  
  .git_sha <- sha(repository_head(repo))
  
  #Save all parameters to csv for reproducibility
  #TODO: write this to a workflow database instead
  saveParams(.parPF)
}

dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))