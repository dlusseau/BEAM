
#' Creates a table of metiers of risk per species from ICES WGBYC dataset (over all years available)
#' 
#' @description
#' `risktable_making` takes bycatch data obtained in the ICES WGBYC annual data call, 
#' and uses those data to generate a table of the metiers of risk for each species.
#' A metier is considered of risk for a species if in the ICES dataset there has 
#' been any bycatch record in any ecoregion.
#' Need to work on the dataset containing data for all years - here extracted in 2026
#' Using the list of species of interest, updated to 2026

risktable_making <- function(bycatch1, ecoreg_species){
  library(dplyr)
  library(readxl)
  
  # a. Metiers of risk over ecoregions and years
  metier_sps <- bycatch1 %>% 
    distinct(species, metierL4) %>% 
    arrange(species, metierL4)
  
  # d. Filter with the list of species of interest
  risktable_assessmentperiod <- metier_sps %>%
    filter(species %in% ecoreg_species$species)
  
  # save results
  return(risktable_assessmentperiod)
  write.csv(risk_17_25_filt, "data/risktable_17_25.csv", row.names = FALSE)
}

