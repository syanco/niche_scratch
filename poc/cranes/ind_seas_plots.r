#!/usr/bin/env Rscript --vanilla
# chmod 744 ind_seas_plots.r #Use to make executable

# This script builds plots of individual- and season-specific niche metrics 
# (univariate).  Anticipates data as formatted by ind_seasonal_moments.r (df of
# component niche means and vars with rows as component distributions)

# ==== Setup ====

'
Builds plots of individual- and season-specific niche metrics

Usage:
ind_seas_plots.r <dat> <out> 
ind_seas_plots.r (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc


#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/niche_scratch'
  rd <- here::here
  
  # path to folder with input .csv's
  .datPF <- file.path(.wd,'analysis/cranes/ind_to_pop/seasonal_seg')
  .outPF <- file.path(.wd,'analysis/cranes/ind_to_pop/figs')
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/funs/input_parse.r'))
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
}


#---- Initialize Environment ----#
t0 <- Sys.time()

source(file.path(.wd,'src/startup.r'))
source(file.path(.wd, "src/poc/cranes/niche_funs.r"))

suppressWarnings(
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gridExtra)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#TODO: check whether I actually like this theme...
theme_set(theme_eda)


#---- Load data ----#
message('Loading data...')

#get list of inout files (1 per sp)
fls <- list.files(.datPF, full.names = T)

#init empty list for results
dat0 <- list()

for(i in 1:length(fls)){
  tryCatch({ # catch empty csv files to let loop continue
    dat0[[i]] <- read.csv(fls[i])
  },error = function(e){cat(glue("ERROR: couldn't read {fls[i]} \n
                                 Moving to nest input file. \n"))})
}


#---- Perform analysis ----#

#-- Calculate individual contributions to population variance
message("\n Calculating indivdual variance contributions")

# init empty list to store species results
sp_c <- list()

for(s in 1:length(dat0)){
  #get species data
  message(glue("Loading data from {fls[s]}"))
  dat_sp <- dat0[[s]] 
  
  # check for no dataconditon and move to next if found 
  if(is.null(dat_sp)){
    message("No records found moving to next species...")
    next}
  
  sp_c[[s]] <- dat_sp %>% 
    drop_na() %>% 
    group_by(season) %>%
    nest() %>% 
    mutate(
      # get the grand mean for each season (mean of individual means)
      grand_mean = map_dbl(data, ~{mean(.x$mean)}),
      # estimate population scale variance
      var_est = map2_dbl(.x=data, 
                         .y=grand_mean, 
                         ~{estPopVar(
                           mean=.x$mean, 
                           var=.x$var, 
                           pop_mean = .y)}),
      # estimate individual contribution to population variance
      ind_components = map2(
        .x=data, # iterate across the dataframes
        .y=grand_mean, # and corresponding pop grand means
        ~{mutate(.x, # mutate each df
                 ind_cont = indContrib( # to add ind contribution to var
                   .x$mean, 
                   .x$var, 
                   n = nrow(.x), 
                   mu_pop = .y),
                 ind_mu_contrib = muContrib(
                   .x$mean, 
                   n = nrow(.x), 
                   mu_pop = .y),
                 ind_var_contrib = var/nrow(.x))})
    ) %>% 
    unnest(c(season)) %>% 
    ungroup() %>% 
    # set season factor order and sort tibble by season
    mutate(season = fct_relevel(season, "Winter", "Spring", "Summer")) %>% 
    arrange(season)
  
} # s

#--  "Coefficient Plots" for Each Species
message("\n Making 'coefficient plots' for each species...")

# code to get min and max values ttto set common x axis across plots
minx <- c()
maxx <- c()
for(s in 1:length(sp_c)){
  if(is.null(sp_c[[s]])){
    message(glue("\n No records found for {fls[s]}, moving to next species..."))
    next
  }else{
    d <- do.call("rbind", sp_c[[s]]$ind_components)
    minx[s] <- min(na.omit(d$ind_cont))*.9
    maxx[s] <- max(na.omit(d$ind_cont))*1.1
    
  }
}



# init empty list to store species plots
sp_c_plot <- list()

for(s in 1:length(sp_c)){
  # check for no dataconditon and move to next if found 
  if(is.null(sp_c[[s]])){
    message(glue("\n No records found for {fls[s]}, moving to next species..."))
    next
  }else{
    message(glue("\n Making plots for {fls[s]}")) 
  }
  
  season_l <- list()
  for(j in 1:4) { # for each season
    message(glue("Plotting {sp_c[[s]]$season[j]}"))
    contribs <- list()
    
    message("Total niche contribution plot...")
    contribs[[1]] <- sp_c[[s]] %>% 
      pluck("ind_components", j) %>% 
      ungroup() %>% 
      mutate(ind_f = as.factor(ind),
             ind_fro = fct_reorder(.f=ind_f, .x=ind_cont, .fun=min),
             id=1:nrow(.),
             id_f = fct_reorder(as.factor(id), ind_cont, min)) %>%   # reset factors
      ggplot()+
      geom_point(aes(x=ind_cont, y = id_f, color = ind_cont))+
      scale_color_viridis_c() +
      geom_vline(data=NULL, aes(xintercept=0))+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none") +
      xlim(c(minx[s], maxx[s]))+
      xlab("")
    
    message("Niche contribution via mean plot...")
    contribs[[2]] <- sp_c[[s]] %>% 
      pluck("ind_components", j) %>% 
      ungroup() %>% 
      mutate(ind_f = as.factor(ind),
             ind_fro = fct_reorder(.f=ind_f, .x=ind_mu_contrib, .fun=min),
             id=1:nrow(.),
             id_f = fct_reorder(as.factor(id), ind_mu_contrib, min)) %>%   # reset factors
      ggplot()+
      geom_point(aes(x=ind_mu_contrib, y = id_f, color = ind_mu_contrib))+
      scale_color_viridis_c() +
      geom_vline(data=NULL, aes(xintercept=0))+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none") +
      xlim(c(minx[s], maxx[s]))+
      xlab("")
    
    message("Niche contribution via variance plot...")
    contribs[[3]] <- sp_c[[s]] %>% 
      pluck("ind_components", j) %>% 
      ungroup() %>% 
      mutate(ind_f = as.factor(ind),
             ind_fro = fct_reorder(.f=ind_f, .x=ind_var_contrib, .fun=min),
             id=1:nrow(.),
             id_f = fct_reorder(as.factor(id), ind_var_contrib, min)) %>%   # reset factors
      ggplot()+
      geom_point(aes(x=ind_var_contrib, y = id_f, color = ind_var_contrib))+
      scale_color_viridis_c() +
      geom_vline(data=NULL, aes(xintercept=0))+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none") +
      xlim(c(minx[s], maxx[s]))+
      xlab("")
    
    season_l[[j]] <- contribs
    
  } # j
  sp_c_plot[[s]] <- season_l
} # s


#-- Write plots to file
message("\n Writing pdfs...")

pdf(file.path(.outPF, glue("seasonal_ind_contrib_{Sys.Date()}.pdf")),
    width = 8, height = 10.5, paper = "letter")

for(s in 1:length(sp_c_plot)){
  # check for no dataconditon and move to next if found 
  if(is.null(sp_c_plot[[s]])){
    message("No records found moving to next species...")
    next
  }else{
    message(glue("\n Combining plots for {fls[s]}")) 
  }
  
  season_plot <-list()  
  
  # make the left column first
  season_plot[[1]] <- arrangeGrob(
    arrangeGrob(sp_c_plot[[s]][[1]][[1]], left = "Total"),
    arrangeGrob(sp_c_plot[[s]][[1]][[2]], left = "Via Mean"),
    arrangeGrob(sp_c_plot[[s]][[1]][[3]], left = "Via Variance"),
    ncol = 1,
    top = as.character(sp_c[[s]]$season[1]))
  
  for(j in 2:4){
    season_plot[[j]] <- arrangeGrob(
      arrangeGrob(sp_c_plot[[s]][[j]][[1]]),
      arrangeGrob(sp_c_plot[[s]][[j]][[2]]),
      arrangeGrob(sp_c_plot[[s]][[j]][[3]]),
      ncol = 1, 
      top = as.character(sp_c[[s]]$season[j]))
    
  }
  
  out <- arrangeGrob(season_plot[[1]],
                     season_plot[[2]],
                     season_plot[[3]],
                     season_plot[[4]], 
                     ncol = 4,
                     #get sp for the plot title
                     #TODO: this is a ridiculous way to do this.  It was very fast to
                     # write it into ind_seasonal_moments.r this way, but holy hell...
                     #just look at that.  I'll never go back and change it, this is 
                     # just here so you know that I know that this is shit.
                     top = sp_c[[s]]$data[[1]]$species[1])
  grid.newpage()
  grid.draw(out)
  
}

dev.off()

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
