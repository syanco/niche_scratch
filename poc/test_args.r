#!/usr/bin/env Rscript 
# chmod 744 script_template.r #Use to make executable

# TODO: Update header and breezy elements


'
Template

Usage:
script_template <vars>
script_template (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-v --vars=<vars> Input variables
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
  vars <- ag$vars
  
  source(rd('src/funs/input_parse.r'))
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
}

print(vars)
print(as.character(vars))
print(c(vars))
print(c(as.character(vars)))