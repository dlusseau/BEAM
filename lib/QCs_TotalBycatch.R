## ---------------------------
##
## Script name: Final QCs in total bycatch estimates and colour coding
##
## Purpose of script: Final quality checks and colour coding
##
## Author: Paula Gutiérrez-Muñoz
## Mail: paula.gutierrez@ieo.csic.es
##
## Date Created: 2024-09-19
##
## R version 4.3.1 (2023-06-16 ucrt) Beagle Scouts
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# Libraries
library(dplyr)

# Data
#total_bycatch <- read.csv("total_bycatch_version3_2021_2023.csv", sep = ";")
#fishing <- read.csv("data/D1_to_D3_post_QC_16Sept_1900/fish_effort_wgbyc_2017_2023_clean_16Sept.csv",
#                   sep = ";")
#monitoring <- read.csv("data/D1_to_D3_post_QC_16Sept_1900/monitor_wgbyc_2017_2023_clean_16Sept.csv",
#                       sep = ";")

# AM: to integrate with the unified beam workflow, I think better use these vars:
# of course, this might require some refactoring of the code below, since beam
# now converts all column names to lowercase.
total_bycatch <- tot2
fishing <- all2
monitoring <- obs2

# 0. Create table with fishing-monitoring DaS -----------------------------
# Table with DaSFishing and DaSMonitored, quantitative and qualitative, per 
#   ecoregion, metier level 4 and vessel length
fishing_eco_metier <- fishing %>% 
  group_by(ecoregion, metierL4, vesselLengthRange) %>% 
  summarize(DaSFish = sum(daysAtSeaF))
fishing_eco_metier$DaSFish_ql <- TRUE

monitoring_eco_metier <- monitoring %>% 
  group_by(ecoregion, metierL4, vesselLengthRange) %>% 
  summarize(DaSMon = sum(daysAtSeaOb))
monitoring_eco_metier$DaSMon_ql <- TRUE

fishing_monitoring_eco_metier <- merge(fishing_eco_metier, monitoring_eco_metier,
                                      by = c("ecoregion", "metierL4", "vesselLengthRange"),
                                      all = TRUE)

fishing_monitoring_eco_metier$DaSFish_ql[is.na(fishing_monitoring_eco_metier$DaSFish_ql)] <- FALSE
fishing_monitoring_eco_metier$DaSMon_ql[is.na(fishing_monitoring_eco_metier$DaSMon_ql)] <- FALSE


# QC1 - MONITORING EFFORT VS FISHING EFFORT -------------------------------
# More monitoring effort than fishing effort (*if dbsg was done correctly this should not be a problem)
# Could happen that when there are simultaneous monitoring (e.g., EM and SO)
# Green = no, fishing effort >= monitoring effort
# Red = yes, monitoring effort > fishing effort
total_bycatch$QC1 <- "green"
for(a in 1:nrow(total_bycatch)){
  subset <- fishing_monitoring_eco_metier[fishing_monitoring_eco_metier$ecoregion %in% 
                                            total_bycatch$ecoregion[a] & 
                                            fishing_monitoring_eco_metier$metierL4 %in% 
                                            total_bycatch$metierL4[a], ]
  if(sum(subset$DaSMon, na.rm = T) >= sum(subset$DaSFish, na.rm = T)){total_bycatch$QC1[a] <- "red"}
}


# QC2 - DATA AVAILABILITY FOR TOTAL BYCATCH ESTIMATION --------------------
# Data availability to calculate total BPUE by ecoregion, metier lvl 4 and species
# Green = if there is both monitoring effort and total fishing effort
# Yellow = if there is NO fishing effort but YES monitoring effort
# Red = if there is NO monitoring effort and NO total fishing effort
total_bycatch$QC2 <- "green"
for(b in 1:nrow(total_bycatch)){
  subset <- fishing_monitoring_eco_metier[fishing_monitoring_eco_metier$ecoregion %in% 
                                           total_bycatch$ecoregion[b] & 
                                            fishing_monitoring_eco_metier$metierL4 %in% 
                                           total_bycatch$metierL4[b], ]
  if(any(subset$DaSMon_ql == TRUE & subset$DaSFish_ql == FALSE)){total_bycatch$QC2[b] <- "yellow"}else if(
    all(subset$DaSMon_ql == FALSE & subset$DaSFish_ql == FALSE)){total_bycatch$QC2[b] <- "red"}
}


# QC3 - FACTORS INFLUENCING BPUE (RE RETAINED IN THE MODEL) ---------------
# Factors influencing BPUE (random effects retained in the model)
# If there is any factor influencing BPUE, do 
#   we have all the levels of the effect (e.g., no BPUE for one country)?
# Green = yes, all levels of the factor retained are available
# Yellow = not all levels available, but random effect is year
# Red = no, no all levels of the factor are available (except when the RE is year) OR
#   sampling protocol does not match monitoring method
# *tb.message reflects this
total_bycatch$QC3 <- "green"
for(c in 1:nrow(total_bycatch)){
  if(total_bycatch$tb.message[c] %in% "random levels for at least one random effect not ok" &
     total_bycatch$model[c] %in% "n_ind ~ 1 + (1 | year)"){total_bycatch$QC3[c] <- "yellow"}
  if((total_bycatch$tb.message[c] %in% "random levels for at least one random effect not ok" ||
         total_bycatch$tb.message[c] %in% "samplingProtocol or monitoringMethod not ok") &  
         total_bycatch$model[c] != "n_ind ~ 1 + (1 | year)"){total_bycatch$QC3[c] <- "red"}
}


# Save results ------------------------------------------------------------
write.csv(total_bycatch, "results/total_bycatch_version3_2021_2023_QCs.csv")


