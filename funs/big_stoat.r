# Function 


big_stoat <- function(data, vars){
  
  if(nrow(data) > 1000){ #OTF annotator can only do 1000 records at once
    
    #number of df segments that will be needed
    n_segs <- ceiling(nrow(data)/1000)
    
    if((nrow(data) %% 1000) == 0){ #if dividing df by 1000 results in 0 remainder
      #generate seqence of segment ending indices
      end_seq <- seq(1000, nrow(data), by = 1000)
    } else { #if nrow(df)/1000 results in remainder
      #generate ending indices and add irregular last segment index
      end_seq <- c(seq(1000, nrow(data), by = 1000), nrow(data))
    }
    #generate sequence of segment beginning indicess
    beg_seq <- seq(1, floor(nrow(data)/1000)*1000+1, 1000)
    
    #init empty lists for segmented data and annotations
    d_list <- list()
    anno_list <- list()
    
    for (i in 1:n_segs){ #loop through segements
      d_list[[i]] <- data[beg_seq[i]:end_seq[i],] #put df segments into list
      # annotate each df
      anno_list[[i]] <- start_annotation_simple(events = d_list[[i]], 
                                                layers = vars)
      Sys.sleep(1)
    }
    
    # anno_list<- lapply(d_list, start_annotation_simple, 
    #                    layers = 'landsat8-evi-100-16', coords = c("lon", "lat"))
    anno_df <- do.call("rbind", anno_list)
  } else {
    anno_df <- start_annotation_simple(events = data, layers = vars)
    
    
  }
  
  return(anno_df)
}


