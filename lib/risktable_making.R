
#' Creates a table of metiers of risk per species from ICES WGBYC dataset (over all years available)
#' 
#' @description
#' `risktable_making` takes bycatch data obtained in the ICES WGBYC annual data call, 
#' and uses those data to generate a table of the metiers of risk for each species.
#' A metier is considered of risk for a species if in the ICES dataset there has 
#' been any bycatch record in any ecoregion.
#' Need to work on the dataset containing data for all years - here extracted in 2026
#' Using the list of species of interest, updated to 2026

library(dplyr)
library(readxl)

# Load the data
fishing <- read.csv("data/D1_wgbyc26.csv", sep = ";")
bycatch_allyears <- read.csv("data/D3_wgbyc26.csv", sep = ";")

# Load the list of species of interest for WGBYC (updated to 2026)
thelist <- read_excel("data/ICES_ETP_bycatch_species_2026.xlsx")


# a. Metiers of risk over ecoregions and years
metier_sps <- bycatch_allyears %>% 
  distinct(species, metierL4) %>% 
  arrange(species, metierL4)

# b. Metiers used in fishing between 2017:2025
fishing_17_25 <- fishing %>% 
  filter(year %in% c(2017:2025)) %>% 
  distinct(metierL4)

# c. From the metiers used for fishing in our period of interest 2017:2025 
#     which are of risk for each sps
risk_17_25 <- metier_sps %>% 
  inner_join(fishing_17_25, by = "metierL4") %>% 
  select(metierL4, species)

# d. Filter with the list of species of interest
risk_17_25_filt <- risk_17_25 %>%
  filter(species %in% thelist$Scientific_name)

write.csv(risk_17_25_filt, "data/risktable_17_25.csv", row.names = FALSE)