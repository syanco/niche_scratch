#!/usr/bin/env Rscript --vanilla
# chmod 744 join_annos.r #Use to make executable

# This script can be used to join annotations back to a mosey_db-style database
# and add then joined table back to the db.  it anticipates a control file that 
# selects which elements of the database should be included in the join actions.

# ==== Setup ====

'
Join annotations to movement data.

Usage:
join_annos.r <db> <out> <ctf> [-b]
join_annos.r (-h | --help)

Parameters:
  dat: path to input csv file. 
  out: path to output directory.
  ctf: path to control file

Options:
-h --help     Show this screen.
-v --version     Show version.
-b --rollback   Rollback transaction if set to true.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch/'
  .rollback <- TRUE
  rd <- here::here
  
  .outPF <- file.path(.wd,'analysis/cranes_anno.csv')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfPF <- file.path(.wd, "analysis/cranes/ctfs/anno_vars.csv")
  
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
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Control Files ----#
ctf <- read_csv(.ctfPF)

#get run rows from ctf
runs_dyn <- ctf %>% 
  filter(run == 1 & dyn == 1)

runs_stat <- ctf %>% 
  filter(run == 1 & dyn == 0)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

#get event table
evt <- tbl(db, "event") %>% 
  collect()
evt$event_id <- as.numeric(evt$event_id)

message("Collecting annotations...")

#-- Dynamic Annotations

if(nrow(runs_dyn) > 0){
  #pull the required tables into memory and store in list
  annos_dyn <- list()
  for(i in 1:nrow(runs_dyn)){
    annos_dyn[[i]] <- tbl(db, runs_dyn$table[i]) %>% 
      filter(variable == !!runs_dyn$variable[i],
             s_buff == !!runs_dyn$s_buff[i],
             t_buff == !!runs_dyn$t_buff[i]) %>% 
      rename_with(.fn = ~c(glue("variable_{runs_dyn$variable[i]}"),
                           glue("s_buff_{runs_dyn$variable[i]}"),
                           glue("product_{runs_dyn$variable[i]}"),
                           glue("total_scenes_{runs_dyn$variable[i]}"),
                           glue("total_value_scenes_{runs_dyn$variable[i]}"),
                           glue("pixel_count_{runs_dyn$variable[i]}"),
                           glue("value_{runs_dyn$variable[i]}"),
                           glue("stdev_{runs_dyn$variable[i]}"),
                           glue("t_buff_{runs_dyn$variable[i]}"),
                           glue("notes_{runs_dyn$variable[i]}")),
                  .cols = c(variable, s_buff, product, total_scenes, 
                            total_value_scenes, pixel_count, value, stdev, t_buff, 
                            notes)) %>% 
      collect()
  } # i
} # if 

#-- Static Annotations

if(nrow(runs_stat) > 0){
  #pull the required tables into memory and store in list
  annos_stat <- list()
  for(i in 1:nrow(runs_stat)){
    annos_stat[[i]] <- tbl(db, runs_stat$table[i]) %>% 
      filter(variable == !!runs_stat$variable[i],
             s_buff == !!runs_stat$s_buff[i],
             t_buff == !!runs_stat$t_buff[i]) %>% 
      rename_with(.fn = ~c(glue("variable_{runs_stat$variable[i]}"),
                           glue("s_buff_{runs_stat$variable[i]}"),
                           glue("product_{runs_stat$variable[i]}"),
                           glue("total_scenes_{runs_stat$variable[i]}"),
                           glue("total_value_scenes_{runs_stat$variable[i]}"),
                           glue("pixel_count_{runs_stat$variable[i]}"),
                           glue("value_{runs_stat$variable[i]}"),
                           glue("stdev_{runs_stat$variable[i]}"),
                           glue("t_buff_{runs_stat$variable[i]}"),
                           glue("notes_{runs_stat$variable[i]}")),
                  .cols = c(variable, s_buff, product, total_scenes, 
                            total_value_scenes, pixel_count, value, stdev, t_buff, 
                            notes)) %>% 
      collect()
  } # i
} # if
message("Disconnecting from databse...")
dbDisconnect(db)


#-- Dynamic Joins
if(exists("annos_dyn")){
  message("Starting dynamic joins...")
  
  # do the first join outside loop
  dyn_join_1 <- annos_dyn[[1]]
  
  #starting with second join, loop through joins
  if(length(annos_dyn) > 1){
    for(i in 2:length(annos_dyn)){
      assign(glue("dyn_join_{i}"),
             get(glue("dyn_join_{i-1}")) %>% 
               full_join(annos_dyn[[i]], by = "event_id", 
                         suffix = c(glue("_{runs_dyn$variable[i-1]}"), glue("_{runs_dyn$variable[1]}")))
      ) #assign
      annos_dyn[[i-1]] <- NA #eliminate previous df from list to save memory
      rm(list=(glue("dyn_join_{i-1}")))
    } #for(i)
  } #if
} else {
  message("No dynamic annotations to join, moving to static annotations...")
}


#-- Static Joins
if(exists("annos_stat")){
  message("Starting static joins...")
  
  # do the first join outside loop
  stat_join_1 <- annos_stat[[1]]
  
  #starting with second join, loop through joins
  if(length(annos_stat) > 1){
    for(i in 2:length(annos_stat)){
      assign(glue("stat_join_{i}"),
             get(glue("stat_join_{i-1}")) %>% 
               full_join(annos_stat[[i]], by = "event_id", 
                         suffix = c(glue("_{runs_stat$variable[i-1]}"), glue("_{runs_stat$variable[1]}")))
      ) #assign
      annos_stat[[i-1]] <- NA #eliminate previous df from list to save memory
      rm(list=(glue("stat_join_{i-1}")))
    } #for(i)
  } #if
  
  #recover temporal duplicates for static layers
  stat_int <- evt %>% 
    select(event_id, lat, lon) %>% #xtract event and lat/lon infor from events 
    right_join(get(glue("stat_join_{length(annos_stat)}"))) %>% # join to annos
    select(-event_id) %>% #drop the events
    full_join(evt, by = c("lat", "lon")) #join to full events with lat/lon
  
  evt_join_final <- stat_int %>% 
    full_join(get(glue("dyn_join_{length(annos_dyn)}")))
} else {
  message("No static annotations to join...")
  evt_join_final <- evt %>% 
    full_join(get(glue("dyn_join_{length(annos_dyn)}")))
}

message("Joins complete...")

#-- Write the joined data back to database
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

dbBegin(db)

dbWriteTable(conn= db, value= evt_join_final, 
             name = glue("anno_join_{Sys.Date()}"), overwrite = T)

#---- Finalize script ----#
if(.rollback) {
  message('Rolling back transaction because this is a test run.')
  dbRollback(db)
} else {
  message(glue("Committing table 'join_{Sys.Date()}' to database."))
  dbCommit(db)
}

message("Disconnecting from databse...")
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
