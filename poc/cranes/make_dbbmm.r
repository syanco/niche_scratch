#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script generates individual dynamic brownian bridge models and associated 
# UDs .

# TODO:  The dBBMM paramaters (e.g., window size, margin, error, etc.) are 
# currently hardcoded.  Could be passed in as options to the script.
# 
# ==== Setup ====

'
Estimate Dynamic Brownian Bridge Movement Models (dBBMs) for pre-segemented 
migratory periods across multiple individuals

Usage:
make_dbbmm.r <db> <out> <nc>
make_dbbmm.r (-h | --help)

Parameters:
  db: path to movement databse. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  rd <- here::here
  
  .outPF <- file.path(.wd,'analysis/cranes/')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .nc <- ag$nc
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'/src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(raster)
    library(move)
    library(doMC)
    library(foreach)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual') %>% 
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)


# # Read in output log
# log <- read_csv(glue("{.outPF}/dbbmm_log.csv"))


registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
# foreach(j = 1:length(inds)) %do% {
foreach(j = 1:length(ind), .errorhandling = "pass", .inorder = F) %dopar% {
  # # check whether individual has been previously considered
  if(file.exists(glue("{.outPF}/dbbmms/dbbmm_{ind[j]}.rdata"))){
    message(glue("Individual {ind[j]} already modelled, moving on..."))
  } else {
    # if(log %>% 
    #    filter(ind_id == ind[j] & year == yearvec[i]) %>% 
    #    nrow() == 0){
    message(glue("Starting individual {ind[j]}..."))  
    
    #---- Initialize database ----#
    message(glue("Initializing database connection for individual {ind[j]}..."))
    
    invisible(assert_that(file.exists(.dbPF)))
    db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
    invisible(assert_that(length(dbListTables(db))>0))
    
    
    #---- Perform analysis ----#
    message(glue("Gathering movement data for individual {ind[j]}..."))
    
    evt0 <- tbl(db, "event")
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    
    message("Filtering data and manipulating dates...")
    
    # prep data
    evt_mod <- evt0 %>% 
      # extract ind
      filter(individual_id == !!ind[j]) %>% 
      # extract year
      # filter(yr == !!yearvec[i]) %>%
      collect() %>% 
      mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T")) %>% 
      # sort by timestamp
      arrange(timestamp)
    
    dbDisconnect(db) 
    
    # TODO: this is an arbitrary minimum... check
    if(nrow(evt_mod) < 25){
      message(glue("Insufficient records for individual {ind[j]}, movin on..."))
      
      # # Make entry in log file
      # outlog <- matrix(c(scientificname, ind[j], studyid, yearvec[i], "dbbm", 
      #                    glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"),
      #                    0, as.character(Sys.Date())), 
      #                  nrow = 1)
      # write.table(outlog, glue("{.outPF}/dbbmm_log.csv"), append = T, row.names = F, 
      #             col.names = F, sep = ",")
      
    } else {
      
      #-- Fit dBBMMs
      message("Estimating dBBMMs...")
      tryCatch({
        # make minimal df for `move`    
        evt_tmp <- evt_mod %>% 
          select(lon, lat, timestamp)
        
        evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, trt = evt_tmp$trt,
                       proj=CRS("+proj=longlat"))
        # burstid <- factor(evt_tmp$wk[1:(n.locs(evt_mv)-1)])
        #id "intended fix rate"
        fixmed <- median(timeLag(x=evt_mv, units="mins"))
        # evt_burst <- burst(evt_mv, burstid)
        evt_mv_t <- spTransform(evt_mv, center = T)
        dbb_var <- brownian.motion.variance.dyn(object = evt_mv_t, 
                                                # TODO check on error modeling
                                                location.error = if(any(is.na(evt_mod$horizontal_accuracy))){
                                                  #TODO: fixed at 5m - maybe reset to the global mean horiz accuracy in the data?
                                                  rep(5, n.locs(evt_mv_t))}else{
                                                    evt_mod$horizontal_accuracy},
                                                
                                                #TODO: think about these...
                                                margin = 11, 
                                                window.size = 31)
        # remove any segments with gaps > 3x the intended fix rate
        dbb_var@interest[timeLag(evt_mv_t,"mins")>(fixmed*3)] <- FALSE
        
        dbbm <- brownian.bridge.dyn(dbb_var, 
                                    #TODO: need to come up with a way to select this more dynamically based on the scale of the data
                                    # raster = 100,
                                    # If we have horiz accuracy use that, otherwise use fixed accuracy
                                    location.error = if(any(is.na(evt_mod$horizontal_accuracy))){
                                      #TODO: fixed at 5m - maybe reset to the global mean horiz accuracy in the data?
                                      rep(5, n.locs(evt_mv_t))}else{
                                        evt_mod$horizontal_accuracy} ,
                                    # location.error = rep(10, n.locs(evt_mv_t)),
                                    time.step = (fixmed/15),
                                    # burstType = levels(burstid),
                                    #TODO: below is somewhat arbitrary
                                    dimSize = 1000,
                                    ext = 10,
                                    margin = 11, window.size = 31)
      }, error = function(e){cat(glue("ERROR: unspecified error in fitting dBBMM for ind {ind[j]}}", 
                                      "\n"))})
      tryCatch({
        if(exists("dbbm")){
          # # Get UD volume
          # vol <- getVolumeUD(dbbm)
          # 
          # # Make 95% mask
          # mask95 <- vol
          # mask95[mask95>0.95] <- NA
          # mask95[mask95<0.95] <- 1
          # 
          # # Clip out 95% UD
          # ud95 <- dbbm*mask95
          
          # write the dbbmm objects to list
          tmp_out <- list("dBBMM Variance" = dbb_var,
                          "dBBMM Object" = dbbm,
                          # "Contours" =  raster2contour(dbbm),
                          # "UD Volume" = vol,
                          # "95% Mask" = mask95,
                          # "95% UD" = ud95,
                          "events" = evt_mod
          )
        } else {
          message(glue("No dbbmm object in memory, nothing written to tmp for individual {ind[j]}"))
        }
      }, error = function(e){cat("ERROR: couldnt write dBBMM objects to tmp", 
                                 "\n")})
      
      #-- Save individual output
      
      message(glue("Writing output for individual {ind[j]} to file..."))
      tryCatch({
        if(exists("tmp_out")){
          save(tmp_out,
               file = glue("{.outPF}/dbbmms/dbbmm_{ind[j]}.rdata")
          )
          
          # # Make entry in log file
          # outlog <- matrix(c(scientificname, ind[j], studyid, yearvec[i], "dbbm", 
          #                    glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"),
          #                    1, as.character(Sys.Date())), 
          #                  nrow = 1)
          # write.table(outlog, glue("{.outPF}/dbbmm_log.csv"), append = T, row.names = F, 
          #             col.names = F, sep = ",")
        }else {
          message(glue("No tmp list in memory, nothing written to file for {ind[j]}"))
        }
      }, error = function(e){cat("ERROR: couldnt save tmp_out to file", 
                                 "\n")})
    } # fi
  } # ifelse file exists
}#i (end loop through years) : #j (end loop through individuals)

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
