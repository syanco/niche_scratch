#!/usr/bin/env Rscript --vanilla
# chmod 744 niche_contrib_by_seg.r #Use to make executable

# This script estimates individual contribution to population niche breadth by 
# season (previously estimated).

# ==== Setup ====

'
Estimates individual niche breadth contributions by seasonal segment

Usage:
niche_contrib_by_seg.r <db> <out> <ctf> 
niche_contrib_by_seg.r (-h | --help)

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
  
  .outPF <- file.path(.wd,'analysis/cranes/niche_contrib')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfs <- file.path(.wd, "ctfs/crane_segmentation")
  
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
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(file.path(.wd,'src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(ctmm)
    library(raster)
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
evt0 <- tbl(db, "event") %>% 
  collect()

species <- tbl(db, "individual") %>% 
  group_by(taxon_canonical_name) %>% 
  summarise(n=n()) %>% 
  collect() %>% 
  filter(n>10,
         !is.na(taxon_canonical_name)) #TODO: this removes 81 NA spp - figure out the species...

tot <- list() # make top level list for storign results
for(s in 1:length(species$taxon_canonical_name)){
  message(glue("Starting species {species$taxon_canonical_name[s]}..."))
  
  ind <- tbl(db, "individual") %>% 
    filter(taxon_canonical_name == !!species$taxon_canonical_name[s]) %>% 
    collect() %>% 
    pull(individual_id)  
  
  # for(j in 1:length(ind)){
  #   print(read_csv(file.path( 
  #     .ctfs,segs[which(str_extract(segs, "[0-9]+") == ind)[1]]
  #   )))}
  
  for(j in 1:length(ind)){
    message(glue("Starting individual {ind[j]}..."))
    
    message("Filtering data and manipulating dates...")
    # prep data
    evt_mod <- evt0 %>% 
      #filter to a single individual
      #TODO: remove   this     once scriiipt moves into production
      filter(individual_id == ind[j]) %>% 
      #remove observations with high groundspeed
      #TODO: eventually remove this and make a separate script that cleans data
      # earlier in the workflow
      filter(is.na(`ground_speed`) | `ground_speed`<10,
             is.na(gps_hdop) | gps_hdop < 5) %>% 
      #create a few new vars:
      dplyr::mutate(week = week(timestamp), #week var
                    year = year(timestamp), #year var
                    WKxYR = glue("{year}{week}"),
                    # timestamp = ymd_hms(timestamp)
                    timestamp = date(timestamp)
      ) #year X week combo var 
    
    #-- Load Segmentation Info
    message("Getting seasonal segmentations...")
    # TODO: I don't know, I suspect there's a more elegant way to do this...
    
    # grab the segmentation file by matching filename to ind
    if(ind[j] %in% str_extract(segs, "[0-9]+")){
      seg_temp <- read_csv(file.path( 
        .ctfs,segs[which(str_extract(segs, "[0-9]+") == ind[j])[1]]
      )) %>%
        arrange(Date) %>% # sort by date
        #remove stopovers
        filter(Status %in% c("Start Fall", "End Fall", "Start Spring", "End Spring"))
      
    } else {
      next
    }
    
    #check for segmentation file data
    if(nrow(seg_temp) < 1){
      message("Insufficient segmentation data, moving to next individual!")
      next
    }
    
    # #check if the first transition starts a stationary period, remove if not
    # if(seg_temp[1,"Status"] == "Start Fall" | seg_temp[1,"Status"] == "Start Spring"){
    #   seg_temp <- seg_temp[-1,]
    # }
    # 
    # #check if the last row finishes a stationary period, remove if not
    # if(nrow(seg_temp) %% 2 != 0){seg_temp <- seg_temp[-nrow(seg_temp),]}
    
    
    # add a variable for the stationary code
    seg_full <- seg_temp %>% 
      mutate(res_period = case_when(Status == "End Fall" ~ "Winter",
                                    Status == "Start Spring" ~ "Winter",
                                    Status == "End Spring" ~ "Summer",
                                    Status == "Start Fall" ~ "Summer"),
             mig_period = case_when(Status == "Start Fall" ~ "Fall",
                                    Status == "End Fall" ~ "Fall",
                                    Status == "Start Spring" ~ "Spring",
                                    Status == "End Spring" ~ "Spring"),
             res_first = res_period == lead(res_period),
             res_last = res_period == lag(res_period),
             mig_first = mig_period == lead(mig_period),
             mig_last = mig_period == lag(mig_period)) 
    
    #- Resident Periods
    
    seg_res <- seg_full %>%
      #remove first or last row if there's an odd number of cuts
      {if(!.$res_first[1]){slice(., -1)}else{.}} %>% 
      {if(!.$res_last[nrow(.)]){slice(., -nrow(.))}else{.}} %>% 
      #drop mig_period col
      select(-mig_period) %>%
      #create id col to keep the pivot working
      mutate(id = ntile(n =nrow(.)/2)) %>% 
      # pivot wider so we can loop through
      pivot_wider(values_from = Date, names_from = Status)
    
    # Winter
    seg_wint <- seg_res %>% 
      filter(res_period == "Winter") 
    
    #Summer
    seg_summ <- seg_res %>% 
      filter(res_period == "Summer") 
      
    #- Migratory Periods
    seg_mig <- seg_full %>%
      #remove first or last row if there's an odd number of cuts
      {if(!.$mig_first[1]){slice(., -1)}else{.}} %>% 
      {if(!.$mig_last[nrow(.)]){slice(., -nrow(.))}else{.}} %>% 
      #drop mig_period col
      select(-res_period) %>%
      #create id col to keep the pivot working
      mutate(id = ntile(n =nrow(.)/2)) %>% 
      # pivot wider so we can loop through
      pivot_wider(values_from = Date, names_from = Status)
    
    seg_fall <- seg_mig %>% 
      filter(mig_period == "Fall") 
    
    seg_spring <- seg_mig %>% 
      filter(mig_period == "Spring") 
    
    #-
    
#TODO:  Pick up here...
    #-- Winter AKDEs
    message("Estimating AKDEs for winter period(s)...")
    
    # init empty list for results
    wint_hr <- list()
    
    #TODO: replace 1's with i's below
    for(i in 1:nrow(seg_wint)){
      #format data as `telemetry` object for ctmm
      evt_tmp <- evt_mod %>% 
        filter(timestamp > seg_wint$`End Fall`[i] & timestamp < seg_wint$`Start Spring`[i]) %>% 
        arrange(timestamp) %>% 
        rename(location.long = lon,
               location.lat = lat,
               individual.local.identifier = individual_id) 
      
      if(nrow(evt_tmp) == 0){
        message("No suitable records during season, moving to next...")
        next 
      } else{
        message(glue("{nrow(evt_tmp)} records found..."))
        evt_telem <- as.telemetry(evt_tmp)
      } 
      
      # get initial acf guess
      guess <- ctmm.guess(evt_telem, interactive = F)
      
      # acf model selection using guess as init
      fit <- ctmm.select(evt_telem, guess)
      
      # fit 
      akde <- akde(evt_telem, fit)
      
      
      # write the akde object to list
      tmp_out <- list()
      tmp_out[[1]] <- akde
      tmp_out[[2]] <- raster(akde, DF = "PDF")
      tmp_out[[3]] <- as.sf(akde, level.UD=.95)
      tmp_out[[4]] <- mask(tmp_out[[2]], tmp_out[[3]])
      wint_hr[[i]] <- tmp_out
    }
    
    #-- Summer AKDEs
    message("Estimating AKDEs for summer period(s)...")
    
    # init empty list for results
    summ_hr <- list()
    
    #TODO: replace 1's with i's below
    for(i in 1:nrow(seg_summ)){
      #format data as `telemetry` object for ctmm
      evt_tmp <- evt_mod %>% 
        filter(timestamp > seg_summ$`End Spring`[i] & timestamp < seg_summ$`Start Fall`[i]) %>% 
        arrange(timestamp) %>% 
        rename(location.long = lon,
               location.lat = lat,
               individual.local.identifier = individual_id) 
      
      if(nrow(evt_tmp) == 0){
        message("No suitable records during season, moving to next...")
        next 
      } else{
        message(glue("{nrow(evt_tmp)} records found..."))
        evt_telem <- as.telemetry(evt_tmp)
      } 
      # get initial acf guess
      guess <- ctmm.guess(evt_telem, interactive = F)
      
      # acf model selction using guess aas init
      fit <- ctmm.select(evt_telem, guess)
      
      # fit 
      akde <- akde(evt_telem, fit)
      
      # write the akde object to list
      tmp_out <- list()
      tmp_out[[1]] <- akde
      tmp_out[[2]] <- raster(akde, DF = "PDF")
      tmp_out[[3]] <- as.sf(akde, level.UD=.95)
      tmp_out[[4]] <- mask(tmp_out[[2]], tmp_out[[3]])
      summ_hr[[i]] <- tmp_out
    }
    
    #-- Save individual output
    
    # declare output destination
    .outTMP <- file.path(.outPF, scientificname, ind[j])
    
    if(length(wint_hr) > 0){
      message(glue("Writing winter UDs for individual {ind[j]} to file..."))
      for(i in 1:length(wint_hr)){
        tryCatch({
          writeRaster(x=wint_hr[[i]][[4]], 
                      format = "GTiff",
                      filename = glue("{.outTMP}/akde_{ind[j]}_wint_{year(seg_wint$`End Fall`[i])}-{year(seg_wint$`Start Spring`[i])}.tif"),
                      overwrite = T)
          # Make entry in log file
          if(file.exists(glue("{.outTMP}/akde_{ind[j]}_wint_{year(seg_wint$`End Fall`[i])}-{year(seg_wint$`Start Spring`[i])}.tif"))){
            outlog <- matrix(c(scientificname, ind[j], "akde", 
                               glue("akde_{ind[j]}_wint_{year(seg_wint$`End Fall`[i])}-{year(seg_wint$`Start Spring`[i])}.tif"),
                               as.character(Sys.Date())), 
                             nrow = 1)
            write.table(outlog, glue("{.outPF}/log.csv"), append = T, row.names = F, col.names = F,
                        sep = ",")
          } #fi
        }, error = function(e){cat("ERROR: unspecified error saving .tiff, probably was not estimated",
                                   "\n")})
      } #i
    } else {#if
      message("No winter UDs to write...") 
    } # else
    
    if(length(summ_hr) > 0){
      message(glue("Writing summer UDs for individual {ind[j]} to file..."))
      for(i in 1:length(summ_hr)){
        tryCatch({
          writeRaster(x=summ_hr[[i]][[4]], 
                      format = "GTiff",
                      filename = glue("{.outTMP}/akde_{ind[j]}_summ_{year(seg_summ$`End Fall`[i])}-{year(seg_summ$`Start Spring`[i])}.tif"),
                      overwrite = T)
          # Make entry in log file
          if(file.exists(glue("{.outTMP}/akde_{ind[j]}_summ_{year(seg_summ$`End Fall`[i])}-{year(seg_summ$`Start Spring`[i])}.tif"))){
            outlog <- matrix(c(scientificname, ind[j], "akde", 
                               glue("akde_{ind[j]}_summ_{year(seg_summ$`End Fall`[i])}-{year(seg_summ$`Start Spring`[i])}.tif"),
                               as.character(Sys.Date())), 
                             nrow = 1)
            write.table(outlog, glue("{.outPF}/log.csv"), append = T, row.names = F, col.names = F,
                        sep = ",")
          } #fi
        }, error = function(e){cat("ERROR: unspecified error saving .tiff, probably was not estimated",
                                   "\n")})
        
      } #i
    } else {#if
      message("No summer UDs to write...")
    } # else
  } #j (end loop through individuals)
} # s en loop through species

#---- Finalize script ----#

message("Disconnecting from databse...")
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
