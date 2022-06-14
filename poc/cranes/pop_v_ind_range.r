###########################
#                         #
# Pop vs Ind Niche Dissim #
#                         #
###########################

# Scratch script to test out ideas of comparing population niche dissim with individual.

# LIBRARIES

library(MVNH)
library(jsonlite)
library(httr)
library(RSQLite)
library(DBI)
library(tidyverse)
library(lubridate)
library(glue)

# FUNCTIONS

get_scenarios <- function () {
  resp <- httr::GET(paste0(BASE_URL, '/list/scenarios'))
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc), simplifyVector = T))
}

get_species_scenario <- function (sciname, product, variable, s_buff, t_buff, limit=50000) {
  resp <- httr::POST(paste0(BASE_URL, '/species/metrics/scatter'), body = list(
    mode = 'temporal',
    scientificname = sciname,
    limit= limit,
    yaxis = list(
      product = product,
      variable = variable,
      temporal = t_buff, 
      spatial = s_buff
    )
  ), encode = "json",)
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc))$rows)
}

`%notin%` <- Negate(`%in%`)

source("src/funs/seg2anno.r")

# INITS

BASE_URL <- 'https://stoat-otf.azurewebsites.net'
enc = "UTF-8"
.dbPF <- file.path('data/anno_move.db')
.anno <- "anno_join_2021-09-08"
.ctfs <- file.path("ctfs/crane_segmentation")


# Pull STOAT Pre-Annotations
#
# View possible annotations
get_scenarios()

# get_species_scenario('Grus grus', 'landsat8', 'evi', 250, 16)

#(while STOAT is doen...)
# Pull in Pop scale niche data
pop <- read.csv("analysis/cranes/Anthropoides virgo_landsat8_evi_500_16.csv") %>% 
  mutate(date = as.Date(eventdate),
         doy = yday(date))


# Read in cranes

#get segmentation filenames
segs <- list.files(.ctfs, pattern = "*.csv")

db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)

# Get target inds from one sp
inds <- tbl(db, "individual") %>%
  filter(taxon_canonical_name == "Anthropoides virgo") %>% 
  pull(individual_id)

# create df of annotations
anno0 <- tbl(db, .anno) %>% 
  filter(individual_id %in% inds) %>%     
  filter(is.na(`ground_speed`) | `ground_speed`<10) %>% 
  filter(is.na(gps_hdop) | gps_hdop < 5) %>% 
  collect()

# add season annotations based on the manual segemntation cutpoints
anno_seasons <- cuts2anno(df = anno0, segs = segs, inds = inds) %>% 
  #only retain good segmentations
  filter(checksum == 1)

anno_res <- anno_seasons %>% 
  filter(winter == 1 | summer == 1)

# get individual dissimilarities

# init empty df
ind_out <- data.frame()

# get reduced list of inds (only those who made it through above)
ind2 <- unique(anno_res$individual_id)

for (i in 1:length(ind2)){
  summ <- anno_res %>% 
    filter(individual_id==ind2[i] & summer == 1) %>% 
    select(`value_derived:evi`) %>% 
    na.omit()
  wint <- anno_res %>% 
    filter(individual_id==ind2[i] & winter == 1)%>% 
    select(`value_derived:evi`) %>% 
    na.omit()
  
  if(nrow(summ) < 1 | nrow(wint) < 1){
    next
  } else {
    
    out <- data.frame(ind = ind2[i],
                      ind_dis = as.numeric(MVNH_dissimilarity(summ, wint)$Bhattacharyya_distance['total']),
                      md = as.numeric(MVNH_dissimilarity(summ, wint)$Mahalanobis_distance['total']),
                      dr = as.numeric(MVNH_dissimilarity(summ, wint)$Determinant_ratio['total']),
                      mean_summ = mean(summ$`value_derived:evi`),
                      mean_wint= mean(wint$`value_derived:evi`)
    )
    
    ind_out <- rbind(ind_out, out)
    
  }
}

# Get phenology
phen <- data.frame()
for(j in 1:length(ind2)){
  # grab the segmentation file by matching filename to ind
  if(ind2[j] %in% str_extract(segs, "[0-9]+")){ #if the file exists...
    seg_temp <- read_csv(file.path( # ...read it in...
      .ctfs,segs[which(str_extract(segs, "[0-9]+") == ind2[j])[1]] # ...matching individual_id
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
  phen <- rbind(phen, seg_temp)
}

# Add day of year
phen <- phen %>% 
  mutate(doy = yday(Date))

# calc median transition days
mssp <- median(phen %>% filter(Status == "Start Spring") %>% pull(doy))
mesp <- median(phen %>% filter(Status == "End Spring") %>% pull(doy))
msfa <- median(phen %>% filter(Status == "Start Fall") %>% pull(doy))
mefa <- median(phen %>% filter(Status == "End Fall") %>% pull(doy))

# get pop-scale seasonal data
gbif_summ <- pop %>% 
  filter(doy > mesp & doy < msfa)
gbif_wint <- pop %>% 
  filter(doy < mssp | doy > mefa)

# calc pop means
pop_stats <- data.frame(season = c("summer", "winter"),
                        mean = c(mean(gbif_summ$value), mean(gbif_wint$value)),
                        variance = c(var(gbif_summ$value), var(gbif_wint$value)),
                        stdv = c(sd(gbif_summ$value), sd(gbif_wint$value)))
pop_mean_summ <- mean(gbif_summ$value)
pop_mean_wint <- mean(gbif_wint$value)

# Calc pop dissimilarity
pop_diss <- MVNH_dissimilarity(gbif_summ %>% select(value),
                               gbif_wint %>% select(value))$Bhattacharyya_distance['total']


# Plots
ggplot()+
  geom_density(data = ind_out, aes(x = ind_dis))+
  geom_vline(aes(xintercept = pop_diss), color = "darkred")+
  theme_minimal()

ind_out_long <- ind_out %>% 
  pivot_longer(c(mean_summ, mean_wint), names_to = "case", values_to = "mean") %>% 
  mutate(season = case_when(case == "mean_summ" ~ "summer",
                            case == "mean_wint" ~ "winter"),
         ind = as.factor(ind))

ggplot(pop_stats)+
  geom_point(aes(x=season, y=mean)) +
  geom_errorbar(aes(x= season, ymin = mean-stdv, ymax = mean+stdv), width = 0.1)+
  theme_minimal()+
  ylim(0,0.5)+
  theme(legend.position = "none")

ggplot(ind_out_long)+
  # geom_point(aes(x=season, y=mean, color = ind)) +
  geom_path(aes(x= season, y = mean, group = ind, color = ind))+
  ylim(0,0.5)+
  theme_minimal()+
  theme(legend.position = "none")

