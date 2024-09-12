# prepare report tables WGBYC 2023

# for preparation of data, see prepare_wgbyc_data

# starting data sets:
# fishing effort 
# RDB: fish_effort_rdb
# WGBYC: fish_effort

# WGBYC monitoring
# monitor

# bycatch



# table xxx effort monitor and bycatch on taxa level -----------

# I create RDB effort 

# previous table contained data from 2021 only (see prepare_report_tables.R)
# 2022-10-24 Sara asked for another selection:
# see mail 2022-10-24 12:51
# and Onenote / Bycatch / WGBYC ToR Collate data on bycatch rates:

# Elasmo & Fish:
#   fishing effort, monitoring effort, bycatch: 2021 only
# 
# Aves, reptilia, mammalia:
# fishing effort 2021
# monitoring effort 2017-2021
# bycatch 2017-2021


# variables
# - fishing effort from RDB and WGBYC
# - monitoring from WGBYC
# - bycatch on taxa level

# sum by 
# EcoRegion and Metier4


# Fishing effort RDB ------

# effort Days at sea
# Number of countries with effort data


# filter rows
# all classes (Elasmo, Fish, Aves, reptilia, mammalia):
# year 2021

fish_effort_rdb_2021 = fish_effort_rdb[Year == 2021,
                                      .(DaS_rdb = sum(DaysAtSea, na.rm = TRUE),
                                        n_country_rdb = uniqueN(Country)),
                                      keyby = .(EcoRegion, MetierL4)]

nrow(fish_effort_rdb_2021)
# 104


# write to file to Gudjon

fwrite(fish_effort_rdb_2021, "results/fish_effort_rdb_ecoreg_m4_2021.csv", sep = ";")



# Fishing effort WGBYC ------

# effort Days at sea
# Number of countries with effort data

# filter rows
# all classes (Elasmo, Fish, Aves, reptilia, mammalia): year 2021


fish_effort_wgbyc_2021 = fish_effort[Year == 2021,
                              .(DaS_wgbyc = sum(DaysAtSeaF, na.rm = TRUE),
                                n_country_wgbyc = uniqueN(Country)),
                              keyby = .(EcoRegion = Ecoregion, MetierL4)]

nrow(fish_effort_wgbyc_2021)
# 205



# Monitoring WGBYC ----

# monitor effort Days at sea
# Number of countries with monitor data

# filter rows
# two different selections:
# 1. elasmo + fish: year 2021
# 2. Aves, reptilia, mammalia: 2017-2021

# 2021
monitor_2021 = monitor[Year == 2021,
                      .(DaS_monitor = sum(DaysAtSeaOb, na.rm = TRUE),
                        n_country_monitor = uniqueN(Country)),
                      keyby = .(EcoRegion = Ecoregion, MetierL4)]

nrow(monitor_2021)
# 117

# 2017-2021
monitor_2017_2021 = monitor[Year %in% 2017:2021,
                           .(DaS_monitor = sum(DaysAtSeaOb, na.rm = TRUE),
                             n_country_monitor = uniqueN(Country)),
                           keyby = .(EcoRegion = Ecoregion, MetierL4)]

nrow(monitor_2017_2021)
# 141



# bycatch WGBYC -----


# filter rows
# two different selections:
# 1. elasmo + fish: year 2021
# 2. Aves, reptilia, mammalia: 2017-2021


# 1. elasmo + fish: year 2021
bycatch_elasmo_fish_2021 = bycatch[Year == 2021 &
                                     class %in% c("Elasmobranchii", "Fish"),
                                   .(n_individ = sum(n_individ, na.rm = TRUE),
                                     n_incident = sum(n_incident, na.rm = TRUE)),
                                   keyby = .(EcoRegion, MetierL4, class)]

# cast to wide
bycatch_elasmo_fish_2021_wd = dcast(bycatch_elasmo_fish_2021, EcoRegion + MetierL4 ~ class,
                                    value.var = c("n_individ", "n_incident"))


# 2. Aves, reptilia, mammalia: 2017-2021
bycatch_aves_rept_mamm_2017_2021 = bycatch[Year %in% 2017:2021 &
                                             class %in% c("Aves", "Reptilia", "Mammalia"),
                                           .(n_individ = sum(n_individ, na.rm = TRUE),
                                             n_incident = sum(n_incident, na.rm = TRUE)),
                                           keyby = .(EcoRegion, MetierL4, class)]

# cast to wide
bycatch_aves_rept_mamm_2017_2021_wd = dcast(bycatch_aves_rept_mamm_2017_2021, EcoRegion + MetierL4 ~ class,
                                            value.var = c("n_individ", "n_incident"))



# get all EcoRegion * MetierL4 combinations ----


# combine the all data sets to get all combinations of EcoRegion and MetierL4
# that are used in any of the tables

ecoreg_m4 = rbindlist(
  list(
    unique(fish_effort_rdb_2021[ , .(EcoRegion, MetierL4)]),
    unique(fish_effort_wgbyc_2021[ , .(EcoRegion, MetierL4)]),
    unique(monitor_2017_2021[ , .(EcoRegion, MetierL4)]),
    unique(bycatch_elasmo_fish_2021[ , .(EcoRegion, MetierL4)]),
    unique(bycatch_aves_rept_mamm_2017_2021[ , .(EcoRegion, MetierL4)])
  ))

ecoreg_m4 = unique(ecoreg_m4[!is.na(EcoRegion)])
nrow(ecoreg_m4)
# 218


setorder(ecoreg_m4, EcoRegion, MetierL4)




# join results on effort (RDB, WGBYC), monitor and bycatch -----

# fishing effort

# RDG
ecoreg_m4[fish_effort_rdb_2021, on = .(EcoRegion, MetierL4),
          `:=`(
            DaS_rdb_2021 = i.DaS_rdb,
            n_country_rdb_2021 = i.n_country_rdb)]


# WGBYC
ecoreg_m4[fish_effort_wgbyc_2021, on = .(EcoRegion, MetierL4),
          `:=`(
            DaS_wgbyc_2021 = i.DaS_wgbyc,
            n_country_wgbyc_2021 = i.n_country_wgbyc)]


# Monitor effort WGBYC, pt 1
# monitor_2021 (elasmo + fish: year 2021)
ecoreg_m4[monitor_2021, on = .(EcoRegion, MetierL4),
          `:=`(
            DaS_monitor_2021 = i.DaS_monitor,
            n_country_monitor_2021 = i.n_country_monitor)]



# Bycatch, number of individuals and incidents, pt 1 
# elasmo + fish: year 2021

ecoreg_m4[bycatch_elasmo_fish_2021_wd, on = .(EcoRegion, MetierL4),
          `:=`(
            n_individ_Elasmobranchii = i.n_individ_Elasmobranchii,
            n_incident_Elasmobranchii = i.n_incident_Elasmobranchii,
            n_individ_Fish = i.n_individ_Fish,
            n_incident_Fish = i.n_incident_Fish)]


# Monitor effort WGBYC, pt 2
# monitor_2017_2021 (Aves, reptilia, mammalia: 2017-2021)

ecoreg_m4[monitor_2017_2021, on = .(EcoRegion, MetierL4),
          `:=`(
            DaS_monitor_2017_2021 = i.DaS_monitor,
            n_country_monitor_2017_2021 = i.n_country_monitor)]



# Bycatch, number of individuals and incidents, pt 2 
# Aves, reptilia, mammalia: 2017-2021

ecoreg_m4[bycatch_aves_rept_mamm_2017_2021_wd, on = .(EcoRegion, MetierL4),
          `:=`(
            n_individ_Aves = i.n_individ_Aves,
            n_incident_Aves = i.n_incident_Aves,
            n_individ_Reptilia = i.n_individ_Reptilia,
            n_incident_Reptilia = i.n_incident_Reptilia,
            n_individ_Mammalia = i.n_individ_Mammalia,
            n_incident_Mammalia = i.n_incident_Mammalia)]

# copy ecoreg_m4
effort_monitor_bycatch_by_ecoreg_m4 = copy(ecoreg_m4)


# # replace NA in numeric columns with zero ----
# 
# # (1) columns with fishing effort and monitoring,
# 
# setnafill(effort_monitor_bycatch_by_ecoreg_m4,
#           type = "const", fill = 0,
#           cols = names(ecoreg_m4)[grep("DaS|n_co", names(ecoreg_m4))])
# 
# # (2) bycatch columns
# # discussed with Sara
# 
# # elasmo & fish
# # leave bycatch as NA, except for:
# # bycatch NA & 
# # 1. fishing effort RDB & WGBYC both 0 OR
# # 2. monitor !is.na & != 0
# 
# # elasmo 1
# effort_monitor_bycatch_by_ecoreg_m4[
#   (is.na(n_individ_Elasmobranchii) &
#     DaS_rdb_2021 == 0 &
#     DaS_wgbyc_2021 == 0),
#   n_individ_Elasmobranchii := 0
# ]
# 
# effort_monitor_bycatch_by_ecoreg_m4[is.na(n_individ_Elasmobranchii) &
#                                       ((!is.na(DaS_rdb_2021) & DaS_rdb_2021 != 0) |
#                                       (!is.na(DaS_wgbyc_2021) & DaS_wgbyc_2021 != 0))]


# write to file
fwrite(effort_monitor_bycatch_by_ecoreg_m4,
       "results/effort_monitor_bycatch_by_ecoreg_m4.csv",
       sep = ";")


######################################################
######################################################







# leftovers ------
# plot

library(ggplot2)

# effort WGBYC
ggplot(ecoreg_m4, aes(x = MetierL4, y = sum_DaS_wgbyc)) +
  facet_wrap(~EcoRegion, scales = "free_x") +
  geom_col(fill = "grey80") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# monitor WGBYC
ggplot(ecoreg_m4, aes(x = MetierL4, y = sum_DaS_monitor)) +
  facet_wrap(~EcoRegion, scales = "free_x") +
  geom_col(fill = "grey80") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# Effort data from RDB to Sara -------

# Sara wants all years (in contrast to above)

effort_rdb = fread("data/Effort data CE from RDB 2017 2021.csv",
                   na.strings = "NULL")


dim(effort_rdb)
# [1] 603498     17

# EcoRegion needs to added from lookup table "Area" - EcoRegion
# ICES_Area_to_EcoRegion.csv

area_ecoregion = fread("data/ICES_Area_to_EcoRegion.csv")

effort_rdb[area_ecoregion, on = c(Area = "ICES_Area"),
           EcoRegion := i.Ecoregion]

# several unmatched rows
effort_rdb[is.na(EcoRegion), .N]
# [1] 5632

# Katja:
# 27.2: set EcoRegion to "Barents Sea"
effort_rdb[Area == "27.2", EcoRegion := "Barents Sea"]

# leave rest as NA

# Remove NA EcoRegion
# Remove North west Atlantic & Mediterranean

effort_rdb2 = effort_rdb[!(EcoRegion %in% c("North West Atlantic", "Mediterranean Sea")) &
                           !is.na(EcoRegion)]


# metier 4 needs to be created from "FishingActivityCategoryEuropeanLvl5"
# Extract string before "_"
effort_rdb2[ , MetierL4 := sub("_.*", "", FishingActivityCategoryEuropeanLvl5)]


# select relevant columns
effort_rdb2 = effort_rdb[ , .(EcoRegion, Country, Year, MetierL4, DaysAtSea)]



# Effort RDB ------

# effort Days at sea
# Number of countries with effort data

effort_rdb_sum = effort_rdb2[ , .(sum_DaS_rdb = sum(DaysAtSea, na.rm = TRUE),
                                  n_country_rdb = uniqueN(Country)),
                              keyby = .(EcoRegion, MetierL4, Year)]

# write to file to Sara

fwrite(effort_rdb_sum, "data/fishing_effort_rdb_ecoreg_m4_year.csv", sep = ";")




# # table yyy effort monitor and bycatch on for selected species and EcoRegions -----------

# Marjorie sent a list of species and Ecoregion
# # Taxa:	EcoRegion Species (common name)
# see email 2022-09-30
# save area_species.csv
# manually added variables to table:
# vernacular name from D3 table
# scientific name


# Create table with one row per species and Metier4