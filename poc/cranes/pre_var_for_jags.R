library(tidyverse)
library(lubridate)
library(glue)

fp <- list.files("analysis/cranes/dbbmms/", full.names = T)

out <- list()

for(i in 1:length(fp)){
  load(fp[i])
  d <- data.frame(tmp_out$events) %>% 
    mutate(doy = yday(timestamp))
  
  v <- tmp_out$`dBBMM Variance`@means
  
  if(length(v) == nrow(d)){
    d$dbbvar <- v  
  } else{
    print(glue("mismatch on {i}"))
    d$dbbvar <- NA
  }
  
  
  out[[i]] <- d
}

tot <- do.call("rbind", out)


# why is this producing only 2 rows...
dat1 <- tot %>% 
  group_by(doy, individual_id) %>% 
  summarize(v = mean(dbbvar, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = doy, values_from = v) %>% 
  as.matrix()

