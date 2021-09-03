#!/usr/bin/env Rscript --vanilla
# chmod 744 plot_niches.r #Use to make executable

# This script calculated time-dynamic niche overlaps across individual animals.
# The overlap is currently measured as the overlap between a time-specific 
# niche and the individual grand mean niche.

# ==== Setup ====

'
Calculate time dynamic niches across multiple individuals

Usage:
plot_niches.r <dat> <out>
plot_niches.r (-h | --help)

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
  
  .outPF <- file.path(.wd,'analysis/cranes')
  .datPF <- file.path(.wd,'analysis/cranes')
  # .ctfPF <- file.path(.wd, "analysis/cranes/ctfs/anno_vars.csv")
  # .anno <- "anno_join_2021-08-17"
  
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
  .datPF <- makePath(ag$dat)
  # .ctfPF <- makePath(ag$ctf)
  # .anno <- as.character(ag$anno)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(rd('src/startup.r'))
#TODO: fix this file path after moving the function script
source(rd("src/poc/cranes/niche_funs.r"))

suppressWarnings(
  suppressPackageStartupMessages({
    library(lubridate)
    library(gridExtra)
    library(ggpubr)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Control Files ----#


#---- Load Data ----#
message("Loading data...")

n_size <- read_csv(file.path(.outPF, "niche_sizes.csv"), col_types = cols(WKxYR = col_character()))
n_dissim <- load(file.path(.outPF, "niche_dissims.Rdata"))
nsd <- read_csv(file.path(.outPF, "nsd.csv"))

#---- Perform analysis ----#
message("Building plots...")

#get index of crane IDs
idx <- unique(n_size$ind)

#-- Dissimilarity Plots
for(i in 1:length(idx)){
  message(glue("Plotting crane {idx[i]} ({i}/{length(idx)})"))
  # Gather data
  dissim_ts <- ind_ts_dat %>% 
    filter(ind == idx[i]) %>%    
    mutate(ts = ymd_hms(ts),
           ind = as.factor(ind)) %>% 
    group_by(WKxYR) %>% 
    summarize(ind = ind,
              ts = ts,
              BD_tot = Bhattacharyya_distance[[1]][1],
              day_temp= Bhattacharyya_distance[[1]][2],
              evi = Bhattacharyya_distance[[1]][3],
              ndvi = Bhattacharyya_distance[[1]][4])
  
  niche_comb <- n_size %>% 
    filter(ind == idx[i]) %>%
    mutate(ts = ymd_hms(ts),
           ind = as.factor(ind)) %>% 
    full_join(dissim_ts, by = c("ind", "WKxYR"))
  
  # produce plot objects
  size_plot <- ggplot(niche_comb) +
    geom_line(aes(x=ts.x, y = total)) +
    ggtitle("Log(niche size)")

  dis_plot <- ggplot(niche_comb) +
    geom_line(aes(x=ts.x, y = BD_tot)) +
    ggtitle("Niche Dissimilarity")

  nsd_plot <- nsd %>%
    filter(individual_id == idx[i]) %>%
    ggplot()+
    geom_line(aes(x=timestamp, y = netSQ))+
    ggtitle("Net Sq. Displacement")

  # combine plots
  comb_plot <- ggarrange(size_plot, dis_plot, nsd_plot,
                          ncol = 1)

  #save plots
  ggsave(filename = file.path(.outPF, glue("niche_plots/np_{idx[i]}.pdf")))
}

message(glue('Script complete in {diffmin(t0)} minutes'))

