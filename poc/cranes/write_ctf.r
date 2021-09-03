      ########################################
      ####                                ####
      ####  Make Annotation Control File  ####
      ####                                ####
      ####        Scott Yanco, PhD        ####
      ####      scott.yanco@yale.edu      ####
      ####                                ####
      ########################################

# Script to generate a control file for joining annotations (used in script: 
      # join_annos.r).  Kind of crummy hard code that will never be fixed b/c
      # really this is meant to just be a bit faster than making the ctf by 
      # hand... Don't @ me.

# define path variables
.wd <- '~/projects/niche_scratch/'
.outPF <- file.path(.wd,"analysis/cranes/ctfs/anno_vars.csv")
.dbPF <- file.path(.wd,'data/anno_move.db')


suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
  }))

#connect to db
invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF)
invisible(assert_that(length(dbListTables(db))>0))

# create vec of tables, remove non-annotation tables 
# !!! NOTE: hardcoded - needs to be updates for different annotation sets and 
#  db structures  !!!
tables <- dbListTables(db)[-c(1,2,8,9,10,11,12,14,15)]

#init empty list
varcomb <- list()

#get unique combos from each table
for(i in 1:length(tables)){
  #as.character(tables[i])
  varcomb[[i]] <- tbl(db, tables[i]) %>%
    #select relevant columns !!! works for STOAT annotations only !!!
    select(product, variable, s_buff, t_buff) %>% 
    distinct() %>% 
    collect() %>% mutate(table = tables[i])
}

#bind them together
ctf <- do.call("rbind", varcomb) %>% 
  mutate(run = 0, # deafult the run status to 0
         dyn = 1) # distinguish dynamic v statis vars - edit this by hand in csv 

#write out the csv that will be the control file
write_csv(ctf, path = .outPF)

#disconnect from db
dbDisconnect(db)
