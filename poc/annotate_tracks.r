#!/usr/bin/env Rscript 
# chmod 744 script_template.r #Use to make executable

# TODO: Update header and breezy elements


'
Template

Usage:
script_template <dat> <out> [--seed=<seed>] [-b] [-t]
script_template (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-b --rollback   Rollback transaction if set to true.
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  .seed <- NULL
  .rollback <- TRUE
  .test <- TRUE
  rd <- here::here
  
  .datPF <- file.path(.wd,'analysis/dat_MOL.csv')
  .outPF <- file.path(.wd,'analysis/track_annotated.csv')
  
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
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
}

#---- Initialize Environment ----#
.seed <- ifelse(is.null(.seed),5326,as.numeric(.seed))

set.seed(.seed)
t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(rstoat)
    library(lubridate)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Load Control files ----#
vars <- read.csv(file.path(.wd, 'analysis/ctfs/stoat_vars.csv')) %>% 
  filter(run == 1) %>% 
  select(var) %>% 
  as_vector("character") %>% 
  unname()

#---- Load data ----#
message('Loading data...')
dat0 <- read.csv(.datPF)

#====

#---- Perform analysis ----#
message("Starting annotation...")
dat0$date <- as_date(dat0$date)

if(nrow(dat0) > 1000){ #OTF annotator can only do 1000 records at once
  
  #number of df segments that will be needed
  n_segs <- ceiling(nrow(dat0)/1000)
  
  if((nrow(dat0) %% 1000) == 0){ #if dividing df by 1000 results in 0 remainder
    #generate seqence of segment ending indices
    end_seq <- seq(1000, nrow(dat0), by = 1000)
  } else { #if nrow(df)/1000 results in remainder
    #generate ending indices and add irregular last segment index
    end_seq <- c(seq(1000, nrow(dat0), by = 1000), nrow(dat0))
  }
  #generate sequence of segment beginning indicess
  beg_seq <- seq(1, floor(nrow(dat0)/1000)*1000+1, 1000)
  
  #init empty lists for segmented data and annotations
  d_list <- list()
  anno_list <- list()
  
  for (i in 1:n_segs){ #loop through segements
    d_list[[i]] <- dat0[beg_seq[i]:end_seq[i],] #put df segments into list
    # annotate each df
    anno_list[[i]] <- start_annotation_simple(d_list[[i]], 
                                              layers = vars, 
                                              coords = c("lon", "lat"),
                                              date = "date")
  }
  
  # anno_list<- lapply(d_list, start_annotation_simple, 
  #                    layers = 'landsat8-evi-100-16', coords = c("lon", "lat"))
  anno_df <- do.call("rbind", anno_list)
} else {
  anno_df <- start_annotation_simple(dat0, vars, 
                                     coords = c("lon", "lat"), date = "date")
  
}

dat_out <- anno_df %>% 
  pivot_wider(id_cols = c(individual_id, timestamp), names_from = variable, 
              values_from = c(s_buff, t_buff, value, stdev, 
                              valid_pixel_count)) %>%
  full_join(dat0)

#---- Save output ---#
message(glue('Saving to {.outPF}'))

dat_out %>% write_csv(.outPF)
 
#---- Finalize script ----#

if(!.test) {
  library(uuid)

  .runid <- UUIDgenerate()
  .parPF <- file.path(.wd,"run_params.csv")

  #Save all parameters to csv for reproducibility
  saveParams(.parPF)
}

message(glue('Script complete in {diffmin(t0)} minutes'))
