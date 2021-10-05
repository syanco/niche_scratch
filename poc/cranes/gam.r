#!/usr/bin/env Rscript --vanilla
# chmod 744 calc_niches.r #Use to make executable

# This script is used to calculate time-dynamic niche accumulation curves.

# ==== Setup ====

'
Calculate time-dynamic niches accumulation curves.

Usage:
niche_accumulation.r <db> <out> <ctf> <anno> 
niche_accumulation.r (-h | --help)

Parameters:
  dat: path to input csv file. 
  out: path to output directory.
  ctf: path to control file
  anno: name of annotation table in db

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
  
  .outPF <- file.path(.wd,'analysis/cranes/niche_accum.csv')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfPF <- file.path(.wd, "analysis/cranes/ctfs/anno_vars.csv")
  .anno <- "anno_join_2021-09-08"
  
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
  .anno <- as.character(ag$anno)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source(rd('src/startup.r'))
#TODO: fix this file path after moving the function script
source(rd("src/poc/cranes/niche_funs.r"))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(geosphere)
    library(mgcv)
    # library(MVNH)
  }))
conflict_prefer("lag", "dplyr")
#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Control Files ----#
ctf <- read_csv(.ctfPF)

#get run rows from ctf
runs <- ctf %>% 
  filter(run == 1)

.variables <- runs %>% 
  select(variable) 

.values <- .variables %>% 
  sapply(., FUN = function(x){glue("value_{x}")}) %>% 
  as.character()

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#
message("Loading annotations...")
anno0 <- tbl(db, .anno)
message("Loading species data...")
sp <- tbl(db, "individual")
message("Disconnecting from databse...")

dat <- anno0 %>% 
  left_join(sp, by = "individual_id") %>% 
  filter(taxon_canonical_id == "Grus grus")

message("Gathering niche data...")

#load data
anno_mod <- dat %>% 
  #remove observations with high groundspeed
  #TODO: eventually remove this and make a separate script that cleans data
  # earlier in the workflow
  filter(`ground_speed`<10) %>% 
  #create a few new vars:
  dplyr::mutate(week = week(timestamp), #week var
                year = year(timestamp), #year var
                WKxYR = glue("{year}{week}"),
                doy = yday(timestamp)) %>%  #year X week combo var 
  rename(evi = `value_derived:evi`, #make some easier to use names
         sp = taxon_canonical_name) %>% 
  mutate(sp_f = as.factor(sp)) %>% # create factor version of species
  group_by(individual_id, doy) %>% 
  summarize(evi_m = mean(evi),
            sp_f = sp_f[1],
            individual_id = individual_id[1],
            indYr = glue("{individual_id}{year}"),
            lon_m = mean(lon),
            lat_m = mean(lat),
            g_speed_m = mean(ground_speed),
            ts = timestamp[1])

#read in NSD data created by `calc_nsd.r`
dat_track_nsd <- read.csv("analysis/cranes/nsd.csv")


##---- Fit som gams by individual ----##

#get ids to work with, only individuals w/ 300+ days
ids <- anno_mod %>% 
  group_by(individual_id) %>% 
  summarize(dur = max(as.Date(ts))-min(as.Date(ts)),
            sp = sp_f[1]) %>% 
  filter(dur>=300)

#check how many of each sp in data
ids %>% 
  group_by(sp) %>% 
  summarise(n = n())


# extract join anno data and nsd data and extract target individuals
anno_Av <- anno_mod %>% 
  left_join(dat_track_nsd, 
            by = c("individual_id"="individual_id", "ts" = "timestamp")) %>% 
  filter(individual_id %in% ids$individual_id) %>% 
  mutate(ind = as.factor(individual_id))

# Fit gams
m0 <- gam(evi_m ~ 1, data = anno_Av, method = 'REML')

m <- gam(evi_m ~ 0 + ind + s(doy, bs = "cc", by = ind),
         data = anno_Av, method = 'REML')



##----      ----#

anno_mod %>%  group_by(sp_f, individual_id) %>% summarise(ind = individual_id[1]) %>% 
  ungroup() %>% 
  group_by(sp_f) %>% 
  summarise(n = n())

m0 <- gam(evi_m ~ 1, data = anno_mod, method = 'REML')

m <- gam(evi_m ~ 0 + sp_f + s(doy, bs = "cc", by = sp_f, k = 5),
         data = anno_mod, method = 'REML')
m_lat <- gam(evi_m ~ 0 + sp_f + s(lat_m, bs = "tp", by = sp_f, k = 5),
         data = anno_mod, method = 'REML')

mRE <- gam(evi_m ~ 0 + sp_f + s(doy, bs = "cc", by = sp_f, k = 5) + s(individual_id, bs = "re"),
         data = anno_mod, method = 'REML')

summary(m)
plot(m_lat)

AIC(m0, m, m_lat)

# function to compare factor smooths
smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}

# fits for doy as pred
pdat <- expand.grid(doy = 1:366,
                    sp_f = unique(anno_mod$sp_f),
                    individual_id = unique(anno_mod$individual_id))
p <- predict(m, newdata = pdat, se.fit = T)

pred_df <- cbind(pdat, p) %>% 
  mutate(lwr = fit - 2*se.fit,
         upr = fit + 2*se.fit)
# fits and dat on grid of plots by spp
ggplot() +
  geom_line(data=anno_mod, 
             aes(x=doy, y = evi_m, color = as.factor(indYr)),
             alpha = 0.5)+
  geom_ribbon(data = pred_df, aes(x = doy, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = pred_df, aes(x = doy, y = fit)) +
  facet_wrap(~ sp_f, ncol = 2) +
  # coord_cartesian(ylim = c(-30,30)) +
  labs(x = "Day of year", y = 'EVI')+
  theme(legend.position = "none")

# only fits, all on one
ggplot() +
  geom_ribbon(data = pred_df, aes(x = doy, ymin = lwr, ymax = upr, group = sp_f), alpha = 0.2) +
  geom_line(data = pred_df, aes(x = doy, y = fit, color = sp_f)) +
  labs(x = "Day of year", y = 'EVI')

# fits for lat as pred
pdat <- expand.grid(lat_m = seq(min(anno_mod$lat_m), max(anno_mod$lat_m), by = 0.1),
                    sp_f = unique(anno_mod$sp_f),
                    individual_id = unique(anno_mod$individual_id))
p <- predict(m_lat, newdata = pdat, se.fit = T)

pred_df <- cbind(pdat, p) %>% 
  mutate(lwr = fit - 2*se.fit,
         upr = fit + 2*se.fit)
# fits and dat on grid of plots by spp
ggplot() +
  geom_line(data=anno_mod, 
            aes(x=lat_m, y = evi_m, color = as.factor(indYr)),
            alpha = 0.5)+
  geom_ribbon(data = pred_df, aes(x = lat_m, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(data = pred_df, aes(x = lat_m, y = fit)) +
  facet_wrap(~ sp_f, ncol = 2) +
  # coord_cartesian(ylim = c(-30,30)) +
  labs(x = "Latitude", y = 'EVI')+
  theme(legend.position = "none")


# only fits, all on one
ggplot() +
  geom_ribbon(data = pred_df, aes(x = lat_m, ymin = lwr, ymax = upr, group = sp_f), alpha = 0.2) +
  geom_line(data = pred_df, aes(x = lat_m, y = fit, color = sp_f)) +
  labs(x = "Day of year", y = 'EVI')

# xp <- predict(m, newdata = pdat, type = 'lpmatrix')

comp1 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Grus grus', 'sp_f')
comp2 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Anthropoides virgo', 'sp_f')
# comp3 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Grus vipeo', 'sp_f')
comp4 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Balearica pavonina', 'sp_f')
comp5 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Anthropoides paradiseus', 'sp_f')
# comp6 <- smooth_diff(m, pdat, 'Grus nigricollis', 'Gruidae', 'sp_f')
comp7 <- smooth_diff(m, pdat, 'Grus grus', 'Anthropoides virgo', 'sp_f')
# comp8 <- smooth_diff(m, pdat, 'Grus grus', 'Grus vipeo', 'sp_f')
comp9 <- smooth_diff(m, pdat, 'Grus grus', 'Balearica pavonina', 'sp_f')
comp10 <- smooth_diff(m, pdat, 'Grus grus', 'Anthropoides paradiseus', 'sp_f')
comp11 <- smooth_diff(m, pdat, 'Anthropoides virgo', 'Balearica pavonina', 'sp_f')
comp12 <- smooth_diff(m, pdat, 'Anthropoides virgo', 'Anthropoides paradiseus', 'sp_f')
comp13 <- smooth_diff(m, pdat, 'Balearica pavonina', 'Anthropoides paradiseus', 'sp_f')

comp <- cbind(doy = 1:366,
              rbind(comp1, comp2, comp4, comp5, comp6, comp7, comp9, comp10, 
                    comp11, comp12, comp13))
comp <- comp[complete.cases(comp),]

ggplot(comp, aes(x = doy, y = diff, group = pair)) +
  # geom_point()
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~ pair, ncol = 2) +
  # coord_cartesian(ylim = c(-30,30)) +
  labs(x = NULL, y = 'Difference in EVI')


## ---- 

# Repeat GAM but use coarse mig strategy coding
mig_key = data.frame(sp_f = c("Grus grus", "Grus nigricollis", "Grus vipeo",
                              "Anthopoides paradiseus", "Anthopoides virgo", 
                              "Balearica pavonina"),
                     mig_stat = as.factor(c("M", "M", "M", "R", "M", "R" )))

anno_mod_1 <- anno_mod %>% 
  left_join(mig_key, by = "sp_f") %>% 
  filter(!is.na(mig_stat))

m.1 <- gam(evi_m ~ mig_stat + s(doy, bs = "cc", by = mig_stat),
         data = anno_mod_1, method = 'REML')


pdat <- expand.grid(doy = 1:366,
                    mig_stat = c("M", "R"))

compMR <- cbind(doy = 1:366, smooth_diff(m.1, pdat, 'M', 'R', 'mig_stat'))

#diffplot
ggplot(compMR, aes(x = doy, y = diff, group = pair)) +
  # geom_point()
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~ pair, ncol = 2) +
  # coord_cartesian(ylim = c(-30,30)) +
  labs(x = NULL, y = 'Difference in EVI')

## ---- 
# Just grab one species and use individual as the factor var

anno_Gg <- anno_mod %>% 
  filter(sp_f == "Anthropoides virgo") %>% 
  mutate(ind_f = as.factor(individual_id))

m0_Gg <- gam(evi_m ~ 1, data = anno_Gg, method = 'REML')

m_Gg <- gam(evi_m ~ 0 + ind_f + s(doy, bs = "cc", by = ind_f),
         data = anno_Gg,  method = 'REML')

AIC(m0_Gg, m_Gg)

## ---- 

#TODO: this is just for testing remove eventually
ind_ID <- c(55754621,234543899,146932155,10836626,234544351, 195375963,963160836)
anno_mod_ind <- anno_mod %>% 
  filter(individual_id %in% ind_ID)

#TODO: Right now the variable is hardcoded - should be able to find a way to get
# this from .values... 

dat <- anno_mod_ind %>% 
  filter(individual_id == ind_ID[1]) %>% 
  rename(evi = `value_derived:evi`)

m <- gam(evi ~ s(doy, bs = "cc"), data = dat, method = 'REML')
m0 <- gam(evi ~ 1, data = dat, method = 'REML')

summary(m)
summary(m0)
x_new <- 1:366
y_pred <- predict(m, data.frame(doy = x_new), interval = "confidence")
pred_df <- data.frame(x_new = x_new, y_pre = y_pred)


datmod <- cbind(dat, y_pred)
ggplot() + 
  geom_point(data = dat, aes(x=doy, y=evi, color = as.factor(year))) + 
  geom_line(data=pred_df, aes(x=x_new, y = y_pred))

gam.check(m)

for(i in 1:length(ind_ID)){
  
}

#---- Finalize script ----#
dbDisconnect(db)

message(glue("Saving output to: {.outPF}"))
write_csv(ind_ts_dat, file = .outPF)

message(glue('Script complete in {diffmin(t0)} minutes'))