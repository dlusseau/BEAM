#make tables on (1) bycatch monitoring effort, (2) number of bycatch

# Observational unit
#Ecoregion, Country, year, species,
# metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol, 

library(data.table)




# prepare monitoring effort ------




monitor = fread("A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/D2_Bycatch_monitoring_effort_Last5Years_v2.csv",
                  na.strings = "NULL")

dim(monitor)
# [1] 17495    30

str(monitor)


# filter monitor data -----

# remove LB (logbook) except for Portugal
monitor = monitor[!(MonitoringMethod == "LB" &
                      Country != "PT")]
nrow(monitor)
# [1] 15479

# remove OTH (other) except for Norway
monitor = monitor[!(MonitoringMethod == "OTH" &
                      Country != "NO")]
nrow(monitor)
# [1] 14668

# remove PO (port observers)
monitor = monitor[(MonitoringMethod != "PO")]
nrow(monitor)
# [1] 12109

# remove VO (vessel observers) for Estonia
monitor = monitor[!(MonitoringMethod == "VO" &
                      Country == "EE")]
nrow(monitor)
# [1] 11835


# check and create variables -----

monitor[ , .N, keyby = Ecoregion]
# OK

monitor[ , .N, keyby = Country]
# OK

monitor[ , .N, keyby = Year]
# OK

monitor[ , .N, keyby = MetierL4]
# oK

monitor[ , .N, keyby = MetierL5]
# OK

 
monitor[ , .N, keyby = VesselLengthRange]

# two strange levels, "F" and "NK"
# set to missing
monitor[VesselLengthRange %in% c("F", "NK"),
        VesselLengthRange := NA]
 
# create vessel length group: below and above 12 meter
# there are two length categories that overlap 12:
# : VL0015, VL0815
# for now, set those to below 12

monitor[VesselLengthRange %in% c(
  "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
  VesselLength_group := "below_12"]


monitor[VesselLengthRange %in% c(
  "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
  VesselLength_group := "above_12"]

monitor[ , .N, keyby = VesselLength_group]
# OK

monitor[ , .N, keyby = SamplingProtocol]

# fix levels
monitor[SamplingProtocol %in% c("ALL", "All", "all"),
        SamplingProtocol := "All"]

monitor[SamplingProtocol == "seabirds",
        SamplingProtocol := "Seabirds"]

# OK


monitor[ , .N, keyby = MonitoringMethod]
# OK



# sum bycatch monitoring effort ------------------

#Observational unit, group by:
#Ecoregion, Country, year,
#metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol,



monitor_sum = monitor[
  ,
  .(DaysAtSea = sum(DaysAtSeaOb, na.rm = TRUE)),
  keyby = .(Ecoregion, Country, Year,
            MetierL4, MetierL5, VesselLength_group,
            SamplingProtocol, MonitoringMethod)]


# For rows with monitoring effort, but no bycatch, add rows with zero bycatch







# prepare bycatch ------

bycatch = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/D3_Bycatch_event_Last5Years_v2.csv",
  na.strings = c("NULL", -9))


dim(bycatch)
# [1] 9126   33

str(bycatch)




# Filter bycatch data on MonitoringMethod & Country -------
# check the filter variables
bycatch[ , .N, keyby = Country]


nrow(bycatch)
# [1]  9126

# remove LB (logbook) except for Portugal
bycatch = bycatch[!(MonitoringMethod == "LB"
                    & Country != "PT")]
nrow(bycatch)
# [1] 8912

# remove OTH (other) except for Norway
bycatch = bycatch[!(MonitoringMethod == "OTH"
                    & Country != "NO")]
nrow(bycatch)
# [1] 8348

# remove PO (port observers)
bycatch = bycatch[(MonitoringMethod != "PO")]
nrow(bycatch)
# [1] 7903

# remove VO (vessel observers) for Estonia
bycatch = bycatch[!(MonitoringMethod == "VO"
                    & Country == "EE")]
nrow(bycatch)
# [1] 7865


# check and create variables -----

# add ecoregion
# read Ecoregion - ICES area data

eco_area = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2022/ToR/C_bycatch_rate/data/ICES_Area_to_EcoRegion.csv")

bycatch[ , .N, keyby = AreaCode]

bycatch[eco_area, on = c(AreaCode = "ICES_Area"),
        Ecoregion := i.Ecoregion]


# check variables

bycatch[ , .N, keyby = Ecoregion]
# OK 

bycatch[ , .N, keyby = Country]
# OK

bycatch[ , .N, keyby = Year]
# OK

bycatch[ , .N, keyby = MetierL4]
# oK

bycatch[ , .N, keyby = MetierL5]
# OK


bycatch[ , .N, keyby = VesselLengthRange]

# two strange levels, "F" and "NK"
# set to missing
bycatch[VesselLengthRange %in% c("F", "NK"),
        VesselLengthRange := NA]

# create vessel length group: below and above 12 meter
# there are two length categories that overlap 12:
# : VL0015, VL0815
# for now, set those to below 12

bycatch[VesselLengthRange %in% c(
  "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
  VesselLength_group := "below_12"]


bycatch[VesselLengthRange %in% c(
  "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
  VesselLength_group := "above_12"]

bycatch[ , .N, keyby = VesselLength_group]
# OK

bycatch[ , .N, keyby = SamplingProtocol]

# fix levels
bycatch[SamplingProtocol %in% c("ALL", "All", "all"),
        SamplingProtocol := "All"]
# OK


bycatch[ , .N, keyby = MonitoringMethod]
# OK


bycatch[ , .N, keyby = Species]
# Mix of ALL CAPS and normal caps


# new species variable with names tolower
bycatch[ , species := tolower(Species)]


# number of individuals catch

bycatch[
  ,
  n_ind := IndividualsWithPingers + IndividualsWithoutPingers]




# sum bycatch  ------------------

#Observational unit, group by:
#Ecoregion, Country, year,
#metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol, species

#Value variables: DaysAtSeaOb bycatch monitor effort DaS, number of individuals

bycatch_sum = bycatch[
  ,
  .(n_ind = sum(n_ind, na.rm = TRUE)),
  keyby = .(Ecoregion, Country, Year,
            MetierL4, MetierL5, VesselLength_group,
            SamplingProtocol, MonitoringMethod, species)]




# create list of relevant Ecoregion * species combinations  #########


# read species list

# file recieved from David 2023-09-20 15:51
annex01_species = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/species list/Annex01_WGBYC_DC_2022_Species_per_Ecoregion.csv")

# select and rename relevant column

annex01_species = annex01_species[ , .(Ecoregion, Species = `Species Scientific name`)]

str(annex01_species)
annex01_species[ , .N, keyby = Ecoregion]


# add Mediterranean species on the Priority list
# email from Paula 2023-09-22

ecoreg_species_med = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/species list/MediterraneanPrioritySpeciesToBeIncluded.csv"
)

str(ecoreg_species_med)
ecoreg_species_med[ , .N, keyby = Ecoregion]


ecoreg_species = rbindlist(
  list(annex01_species, ecoreg_species_med)
)

str(ecoreg_species)
  



# to obtain one row per (relevant) species and ecoregion
# join monitor_sum with Ecoregion * species


monitor_ecoreg_species = monitor_sum[
  ecoreg_species,
  on = .(Ecoregion),
  allow.cartesian = TRUE]

nrow(monitor_sum)
# [1] 1752
nrow(ecoreg_species)
# [1] 1238
nrow(monitor_ecoreg_species)
# [1] 184419

str(monitor_ecoreg_species)
monitor_ecoreg_species[ , .N, keyby = Ecoregion]
monitor_ecoreg_species[ , .N, keyby = Country]
monitor_ecoreg_species[ , .N, keyby = MetierL4]
monitor_ecoreg_species[ , .N, keyby = MetierL5]
monitor_ecoreg_species[ , .N, keyby = VesselLength_group]
monitor_ecoreg_species[ , .N, keyby = SamplingProtocol]
monitor_ecoreg_species[ , .N, keyby = MonitoringMethod]
monitor_ecoreg_species[ , summary(DaysAtSea)]
monitor_ecoreg_species[ , .N, Species]



# create lower caps species
monitor_ecoreg_species[ , species := tolower(Species)]

# join number of individuals bycatch from bycatch_sum
monitor_ecoreg_species[bycatch_sum,
  on = .(Ecoregion, Country, Year,
  MetierL4, MetierL5, VesselLength_group,
  SamplingProtocol, MonitoringMethod,
  species),
  n_ind := i.n_ind]


# replace non-matches (NA rows) with zero
monitor_ecoreg_species[
  is.na(n_ind),
  n_ind := 0]

str(monitor_ecoreg_species)




# check bycatch taxa that are not in 'SamplingProtocol' taxa ###############

# Sara created a species-'sampled taxa' table. Sara: Assigns all 'Cetacea' to 'mammal' -> do the same in 'SamplingProtocol' 
species_taxa = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/species list/candidate BPUE greens.csv")

species_taxa = species_taxa[ , .(Species, taxa)]

# set species tolower
species_taxa[ , species := tolower(Species)]

species_taxa[ , .N, keyby = taxa]
#      taxa  N
# 1:   Bird 31
# 2:   Fish 83
# 3: Mammal 12
# 4: Turtle  2

monitor_ecoreg_species[ , .N, keyby = SamplingProtocol]
#    SamplingProtocol      N
# 1:             <NA>     37
# 2:              All 138773
# 3:        Cetaceans   2019
# 4:             Fish  21029
# 5:          Mammals    965
# 6: ProtectedSpecies  17795
# 7:         Seabirds    765
# 8:          Turtles   1300

# rename levels in species_taxa to match those in monitor_ecoreg
species_taxa[taxa == "Bird", taxa := "Seabirds"]
species_taxa[taxa == "Mammal", taxa := "Mammals"]
species_taxa[taxa == "Turtle", taxa := "Turtles"]


# create a a copy of SamplingProtocol, taxa_monitored
monitor_ecoreg_species[ , taxa_monitored := SamplingProtocol]

# Sara assigned all 'Cetacea' to 'mammal' in species-taxa table -> do the same in 'SamplingProtocol' 
monitor_ecoreg_species[SamplingProtocol == "Cetaceans",
                       taxa_monitored := "Mammals"]


# join taxa to monitor_ecoreg_species
monitor_ecoreg_species[species_taxa, on = .(species),
                       taxa := i.taxa]

monitor_ecoreg_species[ , .N, keyby = taxa]

# create variable indicating if a row is:
# taxa_monitored is 'All' OR 'ProtectedSpecies'
# taxa == taxa_monitored
monitor_ecoreg_species[
  ,
  taxa_bycatch_monitor_ok := (taxa_monitored %in% c("All", "ProtectedSpecies")) |
    (taxa_monitored == taxa)
]

monitor_ecoreg_species[ , .N, keyby = taxa_bycatch_monitor_ok]


str(monitor_ecoreg_species)

fwrite(
  monitor_ecoreg_species,
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/output/monitor_effort_bycatch_ecoregion_species.csv")

########################################################################









# left overs -----



# read species that David needs to re-run
rerun = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/caps_species_needed.csv")

# select rows from monitor & bycatch data which needs to be rerun
monitor_ecoreg_species_rerun = monitor_ecoreg_species[
  species %in% rerun$V2
]

str(monitor_ecoreg_species_rerun)


fwrite(
  monitor_ecoreg_species_rerun,
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/output/monitor_effort_bycatch_ecoregion_species_rerun.csv")





bpue_retain = fread(
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/bpue/BPUE_retain.csv")

# set species to lower
bpue_retain[ , species := tolower(Species)]


# create indicator in bpue_retain where taxa bycatch equals taxa monitored
# taxa_bycatch_monitored_ok
# join on Ecoregion, MetierL4, Species


bpue_retain[ , taxa_bycatch_monitor_ok := TRUE]


bpue_retain[
  monitor_ecoreg_species[!(taxa_bycatch_monitor_ok)],
  on = .(Ecoregion, MetierL4, species),
  taxa_bycatch_monitor_ok :=  i.taxa_bycatch_monitor_ok]


fwrite(
  bpue_retain,
  "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/bpue/bpue_retain_taxa_issue.csv")

