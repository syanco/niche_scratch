# Function to convert segmentation cutpoints into seasonal annotations
# Aticipates segmentation files as output by Fabiola Ianarilli's shiny app
# 
# df = event dataframe - anticipates mosey_db schema
# segs = vector of filepaths to the segmentation files
# inds = vector of individual ids to process, if NULL (default) considers all inds in df


cuts2anno <- function(df, segs, inds = NULL){
  
  require(glue)
  
  # get list of unique individuals
  if(is.null(inds)){
    inds <- unique(df$individual_id)
  }
  # init list to store individual results
  out_ind <- list()
  
  for(j in 1:length(inds)){
    message(glue("Starting individual {inds[j]}..."))
    
    message("Filtering data and manipulating dates...")
    
    # create working df
    evt_mod <- df %>% 
      # extract individual
      filter(individual_id == inds[j]) %>% 
      # create year and date fields
      mutate(date = date(timestamp),
             yr = year(date))
    
    #-- Load Segmentation Info
    message("Getting seasonal segmentations...")
    
    # grab the segmentation file by matching filename to ind
    if(inds[j] %in% str_extract(segs, "[0-9]+")){ #if the file exists...
      seg_temp <- read_csv(file.path( # ...read it in...
        .ctfs,segs[which(str_extract(segs, "[0-9]+") == inds[j])[1]] # ...matching individual_id
      )) %>%
        arrange(Date) %>% # sort by date
        #remove stopovers
        filter(Status %in% c("Start Fall", "End Fall", "Start Spring", "End Spring")) %>% 
        mutate(yr = year(Date)) # add year col
      
    } else {
      # otherwise skip to next individual
      message(glue("No segmentation file found for individual {inds[j]}, moving on!"))
      next
    }
    
    #init list to store ea. years results for an ind
    out_yr <- list()
    
    # Loop through the calendar years in the data
    for(y in 1:length(unique(evt_mod$yr))){
      message(glue("Starting year {y}..."))
      
      # get yr-specific events and segs
      evt_yr <- evt_mod %>%
        filter(yr == unique(evt_mod$yr)[y])
      seg_yr <- seg_temp %>% 
        filter(yr == unique(evt_mod$yr)[y])
      
      #check for segmentation file data - set season to NA if segmentation is empty
      if(nrow(seg_yr) < 1){
        
        message("Insufficient segmentation data, moving to next year!")
        
        evt_out <- evt_yr %>% 
          mutate(winter = NA,
                 spring = NA,
                 summer = NA,
                 fall = NA) %>% 
          mutate(checksum = winter+spring+summer+fall)
        
        out_yr[[y]] <- evt_out  
        
      } else {
        
        # Unpack date thresholds
        ssp <- seg_yr$Date[seg_yr$Status == "Start Spring"]
        esp <- seg_yr$Date[seg_yr$Status == "End Spring"]
        sfa <- seg_yr$Date[seg_yr$Status == "Start Fall"]
        efa <- seg_yr$Date[seg_yr$Status == "End Fall"]
        
        evt_out <- evt_yr
        
        # Use >=/<= for res periods, >/< for mig periods
        # Annotate Winter
        if(length(ssp) > 0 & length(efa) > 0){
          evt_out <- evt_out %>% 
            mutate(winter = case_when(date <= ssp | date >= efa ~ 1,
                                      TRUE ~ 0))
        } else {
          if(length(ssp) > 0){
            evt_out <- evt_out %>% 
              mutate(winter = case_when(date <= ssp ~ 1,
                                        TRUE ~ 0))
          } else{
            if(length(efa) > 0){
              
              evt_out <- evt_out %>% 
                mutate(winter = case_when(date >= efa ~ 1,
                                          TRUE ~ 0))
            } else {
              evt_out <- evt_out %>% 
                mutate(winter = 0)
              
            }
          }
        }
        
        # Annotate Summer 
        if(length(esp) > 0 & length(sfa) > 0){
          evt_out <- evt_out %>% 
            mutate(summer = case_when(date >= esp & date <= sfa ~ 1,
                                      TRUE ~ 0))
        } else {
          if(length(sfa) > 0){
            evt_out <- evt_out %>% 
              mutate(summer = case_when(date <= sfa ~ 1,
                                        TRUE ~ 0))
          } else{
            if(length(esp) > 0){
              
              evt_out <- evt_out %>% 
                mutate(summer = case_when(date >= esp ~ 1,
                                          TRUE ~ 0))
            } else {
              evt_out <- evt_out %>% 
                mutate(summer = 0)
            }
          }
        }
        
        # Annotate Spring
        if(length(ssp) > 0 & length(esp) > 0){
          evt_out <- evt_out %>% 
            mutate(spring = case_when(date > ssp & date < esp ~ 1,
                                      TRUE ~ 0))
        } else {
          if(length(ssp) > 0){
            evt_out <- evt_out %>% 
              mutate(spring = case_when(date > ssp ~ 1,
                                        TRUE ~ 0))
          } else{
            if(length(esp) > 0){
              
              evt_out <- evt_out %>% 
                mutate(spring = case_when(date < esp ~ 1,
                                          TRUE ~ 0))
            } else {
              evt_out <- evt_out %>% 
                mutate(spring = 0)
            }
          }
        }
        
        # Annotate Fall
        if(length(sfa) > 0 & length(efa) > 0){
          evt_out <- evt_out %>% 
            mutate(fall = case_when(date > sfa & date < efa ~ 1,
                                    TRUE ~ 0))
        } else {
          if(length(sfa) > 0){
            evt_out <- evt_out %>% 
              mutate(fall = case_when(date > sfa ~ 1,
                                      TRUE ~ 0))
          } else{
            if(length(efa) > 0){
              
              evt_out <- evt_out %>% 
                mutate(fall = case_when(date < efa ~ 1,
                                        TRUE ~ 0))
            } else {
              evt_out <- evt_out %>% 
                mutate(fall = 0)
            }
          }
        }
        
        # add a check sum
        evt_out <- evt_out %>% 
          mutate(checksum = winter+spring+summer+fall)
        
        #test code
        if(any(evt_out$checksum > 1)){
          warning(glue("Some events were assigned to two or more season (checksum > 2)for inividual {inds[j]}
                       Check records for evt_out$date[evt_out$checksum > 1"))
        }
        
        #stash year results inot year list
        out_yr[[y]] <- evt_out 
        
      } # else - check for any rows in seg_yr
      
    } # y - end year loop
    
    # bind multiple years for a single individual
    out_ind[[j]] <- do.call("rbind", out_yr)
    
    message(glue("Individual {inds[j]} complete..."))
  } # j - end ind loop
  
  # bind multiple individuals back into single df
  out <- do.call("rbind", out_ind)
  
  message("Annotation complete, savign results.")
  return(out)
  
} #function