#!/usr/bin/env Rscript --vanilla
# chmod 744 niche_contrib_by_seg.r #Use to make executable

# This script estimates individual- and season-specific univariate niche 
# parameters (mean and var). Output will be a dataframe where rowsare a single 
# season for a single individual


# For testing use s = 3; j = 7; i = 1
# TODO:
#   * Variable is currently hardcoded, use  strategy from `calc_niche_sizes.r`


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
evt0 <- tbl(db, .anno) %>% 
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
  
  # create list to store btw ind, w/i species results
  sp_ls <- list()
  
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
    if(nrow(seg_temp) < 2){
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
                                    Status == "End Spring" ~ "Spring")) %>% 
      arrange(Date)
    
    # TODO:  This process will still allow segmentation files with the correct
    #   sequence of stop-start, but across years, so encompassing 5 seasons
    #   need to build catch for that Grus grus (sp=3) ind=7 has this situation
    
    # Winter
    seg_wint <- seg_full %>%
      # getb only winter obs
      filter(res_period == "Winter") %>%
      {if(nrow(.) < 2){.}else{
        # must start with season "openning" transition
        {if(.$Status[1] != "End Fall"){slice(., -1)}else{.}} %>% 
          # and end with season closing transition
          {if(.$Status[nrow(.)] != "Start Spring"){slice(., -nrow(.))}else{.}} %>% 
          # create numeric codes for each season and a season complete "checksum"
          mutate(status_code = case_when(Status == "End Fall" ~ 1,
                                         Status == "Start Spring" ~ 2),
                 lagsum = status_code + lag(status_code),
                 leadsum = status_code + lead(status_code)) %>% 
          # remove any ont of place transitions
          filter(!((Status == "End Fall" & (leadsum == 2)) | 
                     (Status == "Start Spring" & (lagsum == 4)))) %>% 
          # create repeating id columnto assist the pivot
          {if(nrow(.) > 1){
            mutate(., id = rep(seq(1, nrow(.)/2, by = 1), each = 2))}else{
              .
            }} %>% 
          # remove unneeded columns
          select(-c(X1, status_code, leadsum, lagsum, mig_period)) %>% 
          # pivot to wide format (1 row per season)
          pivot_wider(values_from = Date, names_from = Status) %>% 
          mutate(run=1)
      }}
    
    #Summer
    seg_summ <- seg_full %>%
      # getb only summer obs
      filter(res_period == "Summer") %>% 
      {if(nrow(.) < 2){.}else{
        #must start with season "openning" transition
        {if(.$Status[1] != "End Spring"){slice(., -1)}else{.}} %>%
          # and end with season closing transition
          {if(.$Status[nrow(.)] != "Start Fall"){slice(., -nrow(.))}else{.}} %>% 
          # create numeric codes for each season and a season complete "checksum"
          mutate(status_code = case_when(Status == "End Spring" ~ 1,
                                         Status == "Start Fall" ~ 2),
                 lagsum = status_code + lag(status_code),
                 leadsum = status_code + lead(status_code)) %>% 
          # remove any ont of place transitions
          filter(!((Status == "End Spring" & (leadsum == 2)) | 
                     (Status == "Start Fall" & (lagsum == 4)))) %>% 
          #create repeating id columnto assist the pivot
          {if(nrow(.) > 1){
            mutate(., id = rep(seq(1, nrow(.)/2, by = 1), each = 2))}else{
              .
            }} %>% 
          # remove unneeded columns
          select(-c(X1, status_code, leadsum, lagsum, mig_period)) %>% 
          #pivot to wide format (1 row per season)
          pivot_wider(values_from = Date, names_from = Status) %>% 
          mutate(run = 1)
      }}
    
    #- Migratory Periods
    seg_spring <- seg_full %>%
      # getb only summer obs
      filter(mig_period == "Spring") %>% 
      {if(nrow(.) < 2){.}else{
        #must start with season "openning" transition
        {if(.$Status[1] != "Start Spring"){slice(., -1)}else{.}} %>% 
          # and end with season closing transition
          {if(.$Status[nrow(.)] != "End Spring"){slice(., -nrow(.))}else{.}} %>% 
          # create numeric codes for each season and a season complete "checksum"
          mutate(status_code = case_when(Status == "Start Spring" ~ 1,
                                         Status == "End Spring" ~ 2),
                 lagsum = status_code + lag(status_code),
                 leadsum = status_code + lead(status_code)) %>% 
          # remove any ont of place transitions
          filter(!((Status == "Start Spring" & (leadsum == 2)) | 
                     (Status == "End Spring" & (lagsum == 4)))) %>% 
          #create repeating id columnto assist the pivot
          {if(nrow(.) > 1){
            mutate(., id = rep(seq(1, nrow(.)/2, by = 1), each = 2))}else{
              .
            }} %>% 
          # remove unneeded columns
          select(-c(X1, status_code, leadsum, lagsum, res_period)) %>% 
          #pivot to wide format (1 row per season)
          pivot_wider(values_from = Date, names_from = Status) %>% 
          mutate(run=1)
      }}
    
    seg_fall <- seg_full %>%
      # getb only summer obs
      filter(mig_period == "Fall") %>% 
      {if(nrow(.) < 2){.}else{
        #must start with season "openning" transition
        {if(.$Status[1] != "Start Fall"){slice(., -1)}else{.}} %>% 
          # and end with season closing transition
          {if(.$Status[nrow(.)] != "End Fall"){slice(., -nrow(.))}else{.}} %>% 
          # create numeric codes for each season and a season complete "checksum"
          mutate(status_code = case_when(Status == "Start Fall" ~ 1,
                                         Status == "End Fall" ~ 2),
                 lagsum = status_code + lag(status_code),
                 leadsum = status_code + lead(status_code)) %>% 
          # remove any ont of place transitions
          filter(!((Status == "Start Fall" & (leadsum == 2)) | 
                     (Status == "End Fall" & (lagsum == 4)))) %>% 
          #create repeating id columnto assist the pivot
          {if(nrow(.) > 1){
            mutate(., id = rep(seq(1, nrow(.)/2, by = 1), each = 2))}else{
              .
            }} %>% 
          # remove unneeded columns
          select(-c(X1, status_code, leadsum, lagsum, res_period)) %>% 
          #pivot to wide format (1 row per season)
          pivot_wider(values_from = Date, names_from = Status) %>% 
          mutate(run = 1)
      }}
    
    
    #-- Winter Niches
    message("Gathering individual stats for winter period(s)...")
    
    # init empty list for results
    wint_out <- list()
    
    if(nrow(seg_wint) > 0 & "run" %in% colnames(seg_wint)){
      for(i in 1:nrow(seg_wint)){
        
        # extract season of interest
        evt_tmp <- evt_mod %>% 
          filter(timestamp > seg_wint$`End Fall`[i] & timestamp < seg_wint$`Start Spring`[i])
        
        # Move to next season in the loop
        if(nrow(evt_tmp) == 0){
          message("No suitable records during season, moving to next...")
          next 
        }  
        
        # write mean and var to tmp
        tmp_out <- data.frame(
          season = "Winter",
          ind = as.character(ind[j]),
          mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
          var = var(na.omit(evt_tmp$`value_derived:evi`)),
          n = nrow(evt_tmp),
          ts = evt_tmp$timestamp
        )
        wint_out[[i]] <- tmp_out
      } #i
    } # fi      
    
    #-- Summer Niches
    message("Gathering individual stats for summer period(s)...")
    
    # init empty list for results
    summ_out <- list()
    
    if(nrow(seg_summ) > 0 & "run" %in% colnames(seg_summ)){
      for(i in 1:nrow(seg_summ)){
        
        #format data as `telemetry` object for ctmm
        evt_tmp <- evt_mod %>% 
          filter(timestamp > seg_summ$`End Spring`[i] & timestamp < seg_summ$`Start Fall`[i])
        
        # Move to next season in the loop
        if(nrow(evt_tmp) == 0){
          message("No suitable records during season, moving to next...")
          next 
        }  
        
        # write mean and var to tmp
        tmp_out <- data.frame(
          season = "Summer",
          ind = as.character(ind[j]),
          mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
          var =var(na.omit(evt_tmp$`value_derived:evi`)),
          n = nrow(evt_tmp),
          ts = evt_tmp$timestamp
        )
        summ_out[[i]] <- tmp_out
      } #i
    } # fi   
    
    #-- Spring Niches
    message("Gathering individual stats for spring period(s)...")
    
    # init empty list for results
    spring_out <- list()
    
    if(nrow(seg_spring) > 0 & "run" %in% colnames(seg_spring)){
      for(i in 1:nrow(seg_spring)){
        
        #format data as `telemetry` object for ctmm
        evt_tmp <- evt_mod %>% 
          filter(timestamp > seg_spring$`Start Spring`[i] & timestamp < seg_spring$`End Spring`[i])
        
        # Move to next season in the loop
        if(nrow(evt_tmp) == 0){
          message("No suitable records during season, moving to next...")
          next 
        }  
        
        # write mean and var to tmp
        tmp_out <- data.frame(
          season = "Spring",
          ind = as.character(ind[j]),
          mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
          var =var(na.omit(evt_tmp$`value_derived:evi`)),
          n = nrow(evt_tmp),
          ts = evt_tmp$timestamp
        )
        spring_out[[i]] <- tmp_out
      } #i
    } # fi   
    
    #-- Fall Niches
    message("Gathering individual stats for fall period(s)...")
    
    # init empty list for results
    fall_out <- list()
    
    if(nrow(seg_fall) > 0 & "run" %in% colnames(seg_fall)){
      for(i in 1:nrow(seg_fall)){
        
        #format data as `telemetry` object for ctmm
        evt_tmp <- evt_mod %>% 
          filter(timestamp > seg_fall$`Start Fall`[i] & timestamp < seg_fall$`End Fall`[i])
        
        # Move to next season in the loop
        if(nrow(evt_tmp) == 0){
          message("No suitable records during season, moving to next...")
          next 
        }  
        
        # write mean and var to tmp
        tmp_out <- data.frame(
          season = "Fall",
          ind = as.character(ind[j]),
          mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
          var =var(na.omit(evt_tmp$`value_derived:evi`)),
          n = nrow(evt_tmp),
          ts = evt_tmp$timestamp
        )
        fall_out[[i]] <- tmp_out
      } #i
    } # fi   
    
    #-- Collate outputs
    message("Collating output...")
    
    ind_ls <- list(
      wint = if(length(wint_out) > 0){do.call("rbind", wint_out)}else{NULL},
      spring = if(length(spring_out) > 0){do.call("rbind", spring_out)}else{NULL},
      summ = if(length(summ_out) > 0){do.call("rbind", summ_out)}else{NULL},
      fall = if(length(fall_out) > 0){do.call("rbind", fall_out)}else{NULL}
    )
    
    #add sp name to the df, catch any null dfs
    sp_ls[[j]] <- do.call("rbind", ind_ls) %>%
      {if(!is.null(.)){
        mutate(., species = species$taxon_canonical_name[s])
      }else{.}}
    
  } #j (end loop through individuals)
  
  #-- Save output per species
  
  message(glue("Writing csv to {.outPF}/{species$taxon_canonical_name[s]}_seas_i2p_{Sys.Date()}.csv"))
  
  out <- do.call("rbind", sp_ls)
  
  write.csv(out, file = glue("{.outPF}/{species$taxon_canonical_name[s]}_seas_i2p_{Sys.Date()}.csv"))
  
} # s en loop through species

#---- Finalize script ----#

message("Disconnecting from databse...")
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
