################################################
#                                              #
#    Download Map of Life (MOL) Annotations    #
#                                              #
################################################

# 

library(jsonlite)
library(httr)
BASE_URL <- 'https://stoat-otf.azurewebsites.net'
enc = "UTF-8"
get_scenarios <- function () {
  resp <- httr::GET(paste0(BASE_URL, '/list/scenarios'))
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc), simplifyVector = T))
}

get_scenarios()

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

get_species_scenario('Grus japonensis', 'landsat8', 'evi', 100, 16)