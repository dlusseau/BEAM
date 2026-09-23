#' Final quality checks and colour coding
#' 
#' @description
#' `add_qc123` Adds columns QC1, QC2 and QC3 to the provided data.table.
#' 
#' @details
#' QC1: Observed vs. total effort makes sense
#'      - green = fishing effort <= monitoring effort for all vessel categories
#'      - red = monitoring effort > fishing effort for one or more vessel categories
#'      
#' QC2: Data availability for total bycatch estimation
#'      - green = there is both monitoring effort and total effort for all vessel categories
#'      - yellow = there is monitoring effort, but no total effort for one or more vessel categories
#'      - red = there is neither monitoring effort, nor total effort for one or more vessel categories
#'  
#' QC3: Availability of levels of random effects (REs) that are part of the model
#'      - green = there is data for all levels of all REs, if the model uses any REs
#'      - yellow = some levels lack data, but the only RE is year
#'      - red = there is no data for one or more levels of at least one RE
#'       
#' @param tot The data.table to be annotated; i.e. the result of a call to `calc_total()`.
#' @param obs A data.table containing all observed fishing effort, usually `obs2`
#' @param all A data.table containing all fishing effort, usually `all2`.
#' @param cols The columns specified when running `calc_bpue`. 
#' @returns A data.table, i.e. a copy of `tot`, with new columns for QC1, QC2 and QC3.
#' @seealso [calc_total()]
#' @export#' Final quality checks and colour coding
#' 
#' @description
#' `add_qc123` Adds columns QC1, QC2 and QC3 to the provided data.table.
#' 
#' @details
#' QC1: Observed vs. total effort makes sense
#'      - green = fishing effort <= monitoring effort for all vessel categories
#'      - red = monitoring effort > fishing effort for one or more vessel categories
#'      
#' QC2: Data availability for total bycatch estimation
#'      - green = there is both monitoring effort and total effort for all vessel categories
#'      - yellow = there is monitoring effort, but no total effort for one or more vessel categories
#'      - red = there is neither monitoring effort, nor total effort for one or more vessel categories
#'  
#' QC3: Availability of levels of random effects (REs) that are part of the model
#'      - green = there is data for all levels of all REs, if the model uses any REs
#'      - yellow = some levels lack data, but the only RE is year
#'      - red = there is no data for one or more levels of at least one RE
#'       
#' @param tot The data.table to be annotated; i.e. the result of a call to `calc_total()`.
#' @param obs A data.table containing all observed fishing effort, usually `obs2`
#' @param all A data.table containing all fishing effort, usually `all2`.
#' @param cols The columns specified when running `calc_bpue`. 
#' @param year_ass The year to be assessed with the QCs. Numeric.
#' @returns A data.table, i.e. a copy of `tot`, with new columns for QC1, QC2 and QC3.
#' @seealso [calc_total()]
#' @export
#' @example
#' Read data:
#' tot <- fread("data/tot1.csv")
#' obs3 <- fread("data/obs3.csv)
#' all2 <- fread("data/all2.csv")
#' Run QC checks:
#' tot2 <- add_qc(tot = tot, obs = obs3, all = all2, cols = c("ecoregion", "year",
#' "country", "metierl4", "metierl5", year_ass = 2023)

add_qc123 <- function(tot, obs, all, cols, year_ass) {
  
  # aggregate data on the columns specified, AND the vessellength column.
  # why do we also split by vessellength?
  # this was discussed during a plenary session at WGBYC 2023. The consensus then 
  # was that including vessellength was a more precautionary approach, since 
  # pooling data from different vessel length groups may not be appropriate.
  
  # calculate das for both dataset per groupping
  qc_obs <- obs[, .(das = sum(daysatsea)), c(cols, "vessellength_group")]
  qc_all <- all[, .(das = sum(daysatseaf)), c(cols, "vessellength_group")]
  # merge
  qc <- merge(qc_obs, qc_all, by = c(cols, "vessellength_group"), all = TRUE, suffixes = c("_obs", "_all"))
  qc <- qc[!is.na(vessellength_group)] # ignore NAs for now.
  qc[is.na(das_obs), das_obs := 0]
  qc[is.na(das_all), das_all := 0]
  
  tot <- copy(tot)
  # these default values come out this way (red, red, green) because of the 
  # logic in the following code
  tot[, c("QC1", "QC2", "QC3") := list("red", "red", "green")]
  
  # QC1 - MONITORING EFFORT VS FISHING EFFORT -------------------------------
  # More monitoring effort than fishing effort (*if dbsg was done correctly this should not be a problem)
  # Could happen that when there are simultaneous monitoring (e.g., EM and SO)
  # Green = fishing effort >= monitoring effort
  # Red = fishing effort < monitoring effort 
  qc[, QC1 := ifelse(all(das_obs <= das_all), "green", "red"), cols]
  
  # keep only data for the year of assessment
  if("year" %in% cols){qc1_summary <- qc[year %in% year_ass][, year := NULL]}
  # tot has less columns than qc, so we need to summarize qc to match tot
  qc1_summary <- qc1_summary[, .(QC1 = ifelse(any(QC1 == "red"), "red", "green")), by = .(ecoregion, metierl4)]
  
  # summarize QC1 to match tot cols -> summarize results per grouping columns (ecoregion, metierl4)
  tot[qc1_summary, on = c("ecoregion", "metierl4"), QC1 := i.QC1]
  #tot[qc, on = cols, QC1 := i.QC1]
  
  # QC2 - DATA AVAILABILITY FOR TOTAL BYCATCH ESTIMATION --------------------
  # Data availability to calculate total BPUE by ecoregion, metier lvl 4 and species
  # Green = if there is both monitoring effort and total fishing effort
  # Yellow = if there is NO fishing effort but YES monitoring effort
  # Red = if there is NO monitoring effort and NO total fishing effort
  qc[, QC2 := ifelse(all(das_obs > 0 & das_all > 0), "green", "yellow"), cols]
  qc[, QC2 := ifelse(all(das_obs == 0 & das_all == 0), "red", QC2), cols]
  
  # keep only data for the year of assessment
  if("year" %in% cols){qc2_summary <- qc[year %in% year_ass][, year := NULL]}
  # summarize QC2 to match tot cols -> summarize results per grouping columns (ecoregion, metierl4)
  qc2_summary <- qc[, .(QC2 = ifelse(any(QC2 == "red"), "red",
                                     ifelse(any(QC2 == "yellow"), "yellow", "green"))), by = .(ecoregion, metierl4)]
  tot[qc2_summary, on = c("ecoregion", "metierl4"), QC2 := i.QC2]
  
  # QC3 - FACTORS INFLUENCING BPUE (RE RETAINED IN THE MODEL) ---------------
  # Factors influencing BPUE (random effects retained in the model)
  # If there is any factor influencing BPUE, do 
  #   we have all the levels of the effect (e.g., no BPUE for one country)?
  # Green = yes, all levels of the factor retained are available
  # Yellow = not all levels available, but random effect is year
  # Red = no, no all levels of the factor are available (except when the RE is year) OR
  #   sampling protocol does not match monitoring method
  # *tb.message reflects this
  tot[message == "random levels for at least one random effect not ok" & model == "n_ind ~ 1 + (1 | year)", QC3 := "yellow"]
  tot[message == "random levels for at least one random effect not ok" | 
        (message == "samplingProtocol or monitoringMethod not ok" & model != "n_ind ~ 1 + (1 | year)"), QC3 := "red"]
  
  tot
}

#' @example
#' Read data:
#' tot <- fread("data/tot.csv")
#' obs2 <- fread("data/obs2.csv)
#' all2 <- fread("data/all2.csv")
#' Run QC checks:
#' tot2 <- add_qc(tot = tot, obs = obs2, all = all2, cols = c("ecoregion", "metierl4"))

# TODO: Why are we passing it obs3? Should we change the function to take obs3 rather than obs2?
#add_qc123 <- function(tot, obs, all, cols) {
    
    # aggregate data on the columns specified, AND the vessellength column.
    # # why do we also split by vessellength?
    # # this was discussed during a plenary session at WGBYC 2023. The consensus then 
    # # was that including vessellength was a more precautionary approach, since 
    # # pooling data from different vessel length groups may not be appropriate.
    # qc_obs <- obs[, .(das = sum(daysatseaob)), c(cols, "vessellength_group")]
    # qc_all <- all[, .(das = sum(daysatseaf)), c(cols, "vessellength_group")]
    
    # qc <- merge(qc_obs, qc_all, by = cols2, all = TRUE, suffixes = c("_obs", "_all"))
    # qc <- qc[!is.na(vessellength_group)] # ignore NAs for now.
    # qc[is.na(das_obs), das_obs := 0]
    # qc[is.na(das_all), das_all := 0]
    
    # tot <- copy(tot)
    # # these default values come out this way (red, red, green) because of the 
    # # logic in the following code
    # tot[, c("QC1", "QC2", "QC3") := list("red", "red", "green")]
    
    # # QC1 - MONITORING EFFORT VS FISHING EFFORT -------------------------------
    # # More monitoring effort than fishing effort (*if dbsg was done correctly this should not be a problem)
    # # Could happen that when there are simultaneous monitoring (e.g., EM and SO)
    # # Green = no, fishing effort >= monitoring effort
    # # Red = yes, monitoring effort > fishing effort
    # qc[, QC1 := ifelse(all(das_obs <= das_all), "green", "red"), cols]
    # tot[qc, on = cols, QC1 := i.QC1]
    
    # # QC2 - DATA AVAILABILITY FOR TOTAL BYCATCH ESTIMATION --------------------
    # # Data availability to calculate total BPUE by ecoregion, metier lvl 4 and species
    # # Green = if there is both monitoring effort and total fishing effort
    # # Yellow = if there is NO fishing effort but YES monitoring effort
    # # Red = if there is NO monitoring effort and NO total fishing effort
    # qc[, QC2 := ifelse(all(das_obs > 0 & das_all > 0), "green", "yellow"), cols]
    # qc[, QC2 := ifelse(all(das_obs == 0 & das_all == 0), "red", QC2), cols]
    # tot[qc, on = cols, QC2 := i.QC2]
    
    # # QC3 - FACTORS INFLUENCING BPUE (RE RETAINED IN THE MODEL) ---------------
    # # Factors influencing BPUE (random effects retained in the model)
    # # If there is any factor influencing BPUE, do 
    # #   we have all the levels of the effect (e.g., no BPUE for one country)?
    # # Green = yes, all levels of the factor retained are available
    # # Yellow = not all levels available, but random effect is year
    # # Red = no, no all levels of the factor are available (except when the RE is year) OR
    # #   sampling protocol does not match monitoring method
    # # *tb.message reflects this
    # tot[message == "random levels for at least one random effect not ok" & model == "n_ind ~ 1 + (1 | year)", QC3 := "yellow"]
    # tot[message == "random levels for at least one random effect not ok" | 
            # (message == "samplingProtocol or monitoringMethod not ok" & model != "n_ind ~ 1 + (1 | year)"), QC3 := "red"]
    
    # tot
# }
# ## ---------------------------
# ##
# ## Script name: Final QCs in total bycatch estimates and colour coding
# ##
# ## Purpose of script: Final quality checks and colour coding
# ##
# ## Author: Paula GutiC)rrez-MuC1oz
# ## Mail: paula.gutierrez@ieo.csic.es
# ##
# ## Date Created: 2024-09-19
# ##
# ## R version 4.3.1 (2023-06-16 ucrt) Beagle Scouts
# ##
# ## ---------------------------
# ##
# ## Notes:
# ##   
# ##
# ## ---------------------------

# # Libraries
# #library(dplyr)
# # 

# # Data
# #total_bycatch <- read.csv("total_bycatch_version3_2021_2023.csv", sep = ";")
# #fishing <- read.csv("data/D1_to_D3_post_QC_16Sept_1900/fish_effort_wgbyc_2017_2023_clean_16Sept.csv",
# #                   sep = ";")
# #monitoring <- read.csv("data/D1_to_D3_post_QC_16Sept_1900/monitor_wgbyc_2017_2023_clean_16Sept.csv",
# #                       sep = ";")

# # AM: to integrate with the unified beam workflow, I think better use these vars:
# # of course, this might require some refactoring of the code below, since beam
# # now converts all column names to lowercase.
# # total_bycatch <- bycatch
# # fishing <- all2
# # monitoring <- obs2

# # 0. Create table with fishing-monitoring DaS -----------------------------
# # Table with DaSFishing and DaSMonitored, quantitative and qualitative, per 
# #   ecoregion, metier level 4 and vessel length
# # 
# # # 0. Create table with fishing-monitoring DaS -----------------------------
# # # Table with DaSFishing and DaSMonitored, quantitative and qualitative, per 
# # #   ecoregion, metier level 4 and vessel length
# # fishing_eco_metier <- fishing %>% 
# #     group_by(ecoregion, metierl4, vessellength_group) %>% 
# #     summarize(DaSFish = sum(daysatseaf))
# # fishing_eco_metier$DaSFish_ql <- TRUE
# # 
# # monitoring_eco_metier <- monitoring %>% 
# #     group_by(ecoregion, metierl4, vessellength_group) %>% 
# #     summarize(DaSMon = sum(daysatseaob))
# # monitoring_eco_metier$DaSMon_ql <- TRUE
# # 
# # fishing_monitoring_eco_metier <- merge(fishing_eco_metier, monitoring_eco_metier,
# #                                        by = c("ecoregion", "metierl4", "vessellength_group"),
# #                                        all = TRUE)
# # 
# # fishing_monitoring_eco_metier$DaSFish_ql[is.na(fishing_monitoring_eco_metier$DaSFish_ql)] <- FALSE
# # fishing_monitoring_eco_metier$DaSMon_ql[is.na(fishing_monitoring_eco_metier$DaSMon_ql)] <- FALSE
# # 
# # 
# # # QC1 - MONITORING EFFORT VS FISHING EFFORT -------------------------------
# # # More monitoring effort than fishing effort (*if dbsg was done correctly this should not be a problem)
# # # Could happen that when there are simultaneous monitoring (e.g., EM and SO)
# # # Green = no, fishing effort >= monitoring effort
# # # Red = yes, monitoring effort > fishing effort
# # total_bycatch$QC1 <- "green"
# # for(a in 1:nrow(total_bycatch)){
# #     subset <- fishing_monitoring_eco_metier[fishing_monitoring_eco_metier$ecoregion %in% 
# #                                                 total_bycatch$ecoregion[a] & 
# #                                                 fishing_monitoring_eco_metier$metierl4 %in% 
# #                                                 total_bycatch$metierl4[a], ]
# #     if(sum(subset$DaSMon, na.rm = T) >= sum(subset$DaSFish, na.rm = T)){total_bycatch$QC1[a] <- "red"}
# # }
# # 
# # 
# # # QC2 - DATA AVAILABILITY FOR TOTAL BYCATCH ESTIMATION --------------------
# # # Data availability to calculate total BPUE by ecoregion, metier lvl 4 and species
# # # Green = if there is both monitoring effort and total fishing effort
# # # Yellow = if there is NO fishing effort but YES monitoring effort
# # # Red = if there is NO monitoring effort and NO total fishing effort
# # total_bycatch$QC2 <- "green"
# # for(b in 1:nrow(total_bycatch)){
# #     subset <- fishing_monitoring_eco_metier[fishing_monitoring_eco_metier$ecoregion %in% 
# #                                                 total_bycatch$ecoregion[b] & 
# #                                                 fishing_monitoring_eco_metier$metierl4 %in% 
# #                                                 total_bycatch$metierl4[b], ]
# #     if(any(subset$DaSMon_ql == TRUE & subset$DaSFish_ql == FALSE)){total_bycatch$QC2[b] <- "yellow"}
# #     
# #     else if(
# #         all(subset$DaSMon_ql == FALSE & subset$DaSFish_ql == FALSE)){total_bycatch$QC2[b] <- "red"}
# # }
# # 
# # 
# # # QC3 - FACTORS INFLUENCING BPUE (RE RETAINED IN THE MODEL) ---------------
# # # Factors influencing BPUE (random effects retained in the model)
# # # If there is any factor influencing BPUE, do 
# # #   we have all the levels of the effect (e.g., no BPUE for one country)?
# # # Green = yes, all levels of the factor retained are available
# # # Yellow = not all levels available, but random effect is year
# # # Red = no, no all levels of the factor are available (except when the RE is year) OR
# # #   sampling protocol does not match monitoring method
# # # *tb.message reflects this
# # total_bycatch$QC3 <- "green"
# # for(c in 1:nrow(total_bycatch)){
# #     if(total_bycatch$message[c] %in% "random levels for at least one random effect not ok" &
# #        total_bycatch$model[c] %in% "n_ind ~ 1 + (1 | year)"){total_bycatch$QC3[c] <- "yellow"}
# #     if((total_bycatch$message[c] %in% "random levels for at least one random effect not ok" ||
# #         total_bycatch$message[c] %in% "samplingProtocol or monitoringMethod not ok") &  
# #        total_bycatch$model[c] != "n_ind ~ 1 + (1 | year)"){total_bycatch$QC3[c] <- "red"}
# # }
# # 
# # 
# # # Save results ------------------------------------------------------------
# # write.csv(total_bycatch, "results/total_bycatch_2021_2023_QCs.csv")
# # 
# # 
