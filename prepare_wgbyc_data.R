library(data.table)


# prepare fishing effort ------

fish_effort = fread("data/D1_Fishing_effort_27092022_v6.csv",
                    na.strings = "NULL")

dim(fish_effort)
# [1] 128321     23

str(fish_effort)

names(fish_effort)

fish_effort


# check variables

fish_effort[ , .N, keyby = Ecoregion]

#                       Ecoregion     N
#                          <char> <int>
#  1:                        <NA>    86
# 13:                 Iceland Sea   555
# 14:  Ionian Sea and the Central   321
# 15: Ionian Sea and the Central   1713
# 16:           Mediterranean Sea   129
# 17:         North West Atlantic   275

# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions

# Exclude 'North West Atlantic' & 'Mediterranean Sea'
# Change Icelandic Sea to Icelandic Waters
# Change 'Ionian Sea and the Central' to 'Ionian Sea and the Central Mediterranean Sea'

nrow(fish_effort)

fish_effort = fish_effort[!is.na(Ecoregion) &
                            !(Ecoregion %in% c("North West Atlantic", "Mediterranean Sea"))]
fish_effort[Ecoregion == "Iceland Sea", Ecoregion := "Icelandic Waters"]
fish_effort[startsWith(Ecoregion, "Ionian"), Ecoregion := "Ionian Sea and the Central Mediterranean Sea"]

nrow(fish_effort)
# [1] 127831


# prepare monitoring effort ------

monitor = fread("data/D2_Monitoring_Effort_27092022_v6.csv",
                  na.strings = "NULL")

dim(monitor)
# [1] 16448    30

str(monitor)


# filter monitor (D2) -----

# Rules from Carlos

# Filter According to MonitoringMethod & Country

# SQL filter from Carlos
# <> means !=

# remove LB (logbook) except for Portugal
#(table_D3.MonitoringMethod <> 'LB' or tblUpload_2.country = 'PT') 

# remove OTH (other) except for Norway
# (table_D3.MonitoringMethod <> 'OTH' or tblUpload_2.country = 'NO') 

# remove PO (port observers)
# (table_D3.MonitoringMethod <> 'PO') 

# remove VO (vessel observers) for Estonia
# ( tblUpload_2.country <> 'EE' or table_D3.MonitoringMethod <> 'VO')


nrow(monitor)
# [1] 16448

# remove LB (logbook) except for Portugal
monitor = monitor[!(MonitoringMethod == "LB" & Country != "PT")]
nrow(monitor)
# [1] 15297

# remove OTH (other) except for Norway
monitor = monitor[!(MonitoringMethod == "OTH" & Country != "NO")]
nrow(monitor)
# [1] 14622

# remove PO (port observers)
monitor = monitor[(MonitoringMethod != "PO")]
nrow(monitor)
# [1] 12554

# remove VO (vessel observers) for Estonia
monitor = monitor[!(MonitoringMethod == "VO" & Country == "EE")]
nrow(monitor)
# [1] 12280


# check Ecoregion -----


monitor[ , .N, keyby = Ecoregion]
#                       Ecoregion     N
#                           <char> <int>
#  1:                        <NA>    12
# 13:                 Iceland Sea   219
# 14:  Ionian Sea and the Central     3
# 15: Ionian Sea and the Central    425
# 16:           Mediterranean Sea    57
# 17:         North West Atlantic   183

# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions

# Exclude missing, 'North West Atlantic' & 'Mediterranean Sea'
# Change Icelandic Sea to Icelandic Waters
# Change 'Ionian Sea and the Central' to 'Ionian Sea and the Central Mediterranean Sea'

monitor = monitor[!is.na(Ecoregion) &
                    !(Ecoregion %in% c("North West Atlantic", "Mediterranean Sea"))]
monitor[Ecoregion == "Iceland Sea", Ecoregion := "Icelandic Waters"]
monitor[startsWith(Ecoregion, "Ionian"), Ecoregion := "Ionian Sea and the Central Mediterranean Sea"]


# prepare bycatch ------

bycatch = fread("data/D3_BycatchEvent_27092022_v6.csv",
                  na.strings = c("NULL", -9))

dim(bycatch)
# [1] 9472   34

str(bycatch)

bycatch


# Filter bycatch (D3) data on MonitoringMethod & Country ------

# Rules from Carlos
# SQL filter from Carlos
# <> means !=

# remove LB (logbook) except for Portugal
#(table_D3.MonitoringMethod <> 'LB' or tblUpload_2.country = 'PT') 

# remove OTH (other) except for Norway
# (table_D3.MonitoringMethod <> 'OTH' or tblUpload_2.country = 'NO') 

# remove PO (port observers)
# (table_D3.MonitoringMethod <> 'PO') 

# remove VO (vessel observers) for Estonia
# ( tblUpload_2.country <> 'EE' or table_D3.MonitoringMethod <> 'VO')


nrow(bycatch)
# [1]  9472

# remove LB (logbook) except for Portugal
bycatch = bycatch[!(MonitoringMethod == "LB" & Country != "PT")]
nrow(bycatch)
# [1] 9348

# remove OTH (other) except for Norway
bycatch = bycatch[!(MonitoringMethod == "OTH" & Country != "NO")]
nrow(bycatch)
# [1] 8796

# remove PO (port observers)
bycatch = bycatch[(MonitoringMethod != "PO")]
nrow(bycatch)
# [1] 8502

# remove VO (vessel observers) for Estonia
bycatch = bycatch[!(MonitoringMethod == "VO" & Country == "EE")]
nrow(bycatch)
# [1] 8464


# check Ecoregion -----

bycatch[ , .N, keyby = EcoRegion]
#                     EcoRegion     N
#                        <char> <int>
#                   Iceland Sea   309
# 12: Ionian Sea and the Central  137
# 13:         Mediterranean Sea    53
# 14:       North West Atlantic    68


# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions

# Exclude miossing, 'North West Atlantic' & 'Mediterranean Sea'
# Change Icelandic Sea to Icelandic Waters
# Change 'Ionian Sea and the Central' to 'Ionian Sea and the Central Mediterranean Sea'

bycatch = bycatch[!is.na(EcoRegion) &
                    !(EcoRegion %in% c("North West Atlantic", "Mediterranean Sea"))]
bycatch[EcoRegion == "Iceland Sea", EcoRegion := "Icelandic Waters"]
bycatch[startsWith(EcoRegion, "Ionian"), EcoRegion := "Ionian Sea and the Central Mediterranean Sea"]

nrow(bycatch)
# [1] 8343


# create new class in bycatch -------

# Katja:
# Bird
# Elasmo
# Reptile: Superclass == "Reptilia" & is.na(classname)
# Fish: Collapse Actinopteri + Holocephali + Myxini + Petromyzonti to fish

# first, copy classname to new variable class
bycatch[ , class := classname]

bycatch[SuperClass == "Reptilia" & is.na(classname), class := "Reptilia"]

bycatch[classname %in% c("Actinopteri", "Holocephali", "Myxini", "Petromyzonti"),
         class := "Fish"]


# sum pinger and non-pinger for individuals and incidents ------

# number of individuals
bycatch[ , n_individ := IndividualsWithPingers + IndividualsWithoutPingers]

# number of incidents
bycatch[ , n_incident := IncidentsWithPingers + IncidentsWithoutPingers]


# effort data from RDB ------

fish_effort_rdb = fread("data/Effort data CE from RDB 2017 2021.csv",
                   na.strings = "NULL")

dim(fish_effort_rdb)
# [1] 603498     17

# EcoRegion needs to added from lookup table "Area" - EcoRegion
# ICES_Area_to_EcoRegion.csv

area_ecoregion = fread("data/ICES_Area_to_EcoRegion.csv")

fish_effort_rdb[area_ecoregion, on = c(Area = "ICES_Area"),
           EcoRegion := i.Ecoregion]

# several non-matched rows
fish_effort_rdb[is.na(EcoRegion), .N]
# [1] 5632

# unmatched Areas
fish_effort_rdb[is.na(EcoRegion), sort(unique(Area))]
# [1] "21.1.A"   "21.1.B"   "21.1.C"   "21.1.D"  
# [5] "21.4.V"   "21.6.G"   "21.6.H"   "34"      
# [9] "34.1.1"   "34.1.1.3" "34.1.2"   "34.1.3"  
# [13] "34.1.3.1" "34.1.3.2" "34.2"     "34.3.1"  
# [17] "34.3.1.1" "34.3.1.2" "34.3.1.3" "34.3.3"  
# [21] "34.3.6"   "47.1"     "47.1.1"   "47.1.2"  
# [25] "47.1.3"   "47.1.4"   "47.A.1"   "87"      
# [29] "87.2.6"   "87.3.3"

# Katja:
# 27.2: set EcoRegion to "Barents Sea"
# leave rest as NA

# remove NA, North West Atlantic, Med. Sea

fish_effort_rdb[Area == "27.2", EcoRegion := "Barents Sea"]

fish_effort_rdb = fish_effort_rdb[!is.na(EcoRegion) & 
                                    !(EcoRegion %in% c("North West Atlantic", "Mediterranean Sea"))]


# metier 4 needs to be created from "FishingActivityCategoryEuropeanLvl5"
# Extract string before "_"

fish_effort_rdb[ , MetierL4 := sub("_.*", "", FishingActivityCategoryEuropeanLvl5)]

nrow(fish_effort_rdb)
# [1] 597163
# select relevant rows and columns


# save files as csv -----

fwrite(fish_effort_rdb, "data/clean/fish_effort_rdb.csv", sep = ";")
fwrite(fish_effort, "data/clean/fish_effort_wgbyc.csv", sep = ";")
fwrite(monitor, "data/clean/monitor_wgbyc.csv", sep = ";")
fwrite(bycatch, "data/clean/bycatch_wgbyc.csv", sep = ";")







# check variables ------

# check effort --------
fish_effort[ , .N, tblUploadID]
#OK

fish_effort[ , .N, keyby = Year]
# OK

fish_effort[ , .N, keyby = Month]
# 12382 NA

fish_effort[ , .N, keyby = Quarter]
# 1840 NA

fish_effort[ , .N, keyby = AreaType]
# OK

fish_effort[ , .N, keyby = AreaCode]
# OK

fish_effort[ , .N, keyby = MetierL3]
# OK

fish_effort[ , .N, keyby = MetierL4]
# OK

fish_effort[ , .N, keyby = MetierL5]
# OK

fish_effort[ , .N, keyby = MetierL6]
# <NA> 58585
# -  1141
# set "-" to NA

fish_effort[ , .N, keyby = VesselLengthRange]
# OK levels
# Set Unknown to NA?


fish_effort[ , summary(DaysAtSeaF)]
# OK?

fish_effort[ , .N, keyby = VesselsF]
# strange values
# read as integer64
# what kind of variable is this?

fish_effort[ , summary(VesselsF)]

fish_effort[ , summary(TripsF)]
# OK?
# veery high max

fish_effort[ , all.equal(TripsF, TripsF_txt)]
# equal to TripsF


fish_effort[ , summary(TotalNetsLengthF)]
# OK?
# 126104 NA

fish_effort[ , summary(TotalKmHoursF)]
# OK?

fish_effort[ , summary(NoofHaulsF)]
# 115715 NA
# very high max


fish_effort[ , summary(TotalTowTimeF)]
# very high max

fish_effort[ , .N, keyby = Country]
# OK

fish_effort[ , .N, keyby = Organization]


# check monitor ----------

monitor[ , .N, keyby = SubmissionYear]
# OK

monitor[ , .N, keyby = tblUploadID]
# OK

monitor[ , .N, keyby = Country]
# OK  

monitor[ , .N, keyby = Organization]
# OK

monitor[ , .N, keyby = Year]
# OK

monitor[ , .N, keyby = Month]
# OK

monitor[ , .N, keyby = AreaType]
# OK

monitor[ , .N, keyby = AreaCode]
# OK

monitor[ , .N, keyby = MetierL3]
# OK

monitor[ , .N, keyby = MetierL4]
# OK

monitor[ , .N, keyby = MetierL5]
# OK

monitor[ , .N, keyby = MetierL6]
# Many NA

monitor[ , .N, keyby = VesselLengthRange]
# F?

monitor[ , .N, keyby = MonitoringProgramType]
# OK

monitor[ , .N, keyby = SamplingProtocol]
# OK

monitor[ , .N, keyby = MonitoringMethod]
# OK

monitor[ , .N, keyby = TargetSpecies]

# na.string: #N/A

monitor[ , .N, keyby = PingerCharacteristics]
# OK

monitor[ , .N, keyby = OthMitigationMeasures]
# OK 

monitor[ , summary(VesselsOb)]
# OK

monitor[ , summary(TripsOb)]
# Very high max!

monitor[ , summary(DaysAtSeaOb)]
# Very high max!

monitor[ , summary(HaulsWithPingersOb)]
# OK

monitor[ , summary(TotalNetsLengthOb)]
# very high max!

monitor[ , summary(TotalKmHoursOb)]
# very high max

monitor[ , summary(NoofHaulsOb)]
# OK

monitor[ , summary(TotalTowTimeOb)]
# very high max

monitor[ , .N, keyby = Monitoring812Type]
# OK


# check bycatch ----------
bycatch[ , .N, keyby = SubmissionYear]
# OK

bycatch[ , .N, keyby = tblUploadID]
# OK

bycatch[ , .N, keyby = Country]
# OK

bycatch[ , .N, keyby = Organization]
# OK

bycatch[ , .N, keyby = RecordType]
# OK

bycatch[ , .N, keyby = Year]
# OK

bycatch[ , .N, keyby = Month]
# OK

bycatch[ , .N, keyby = AreaType]
# OK

bycatch[ , .N, keyby = AreaCode]
# OK?
# e.g. 1~5~6

bycatch[ , .N, keyby = MetierL3]
# OK

bycatch[ , .N, keyby = MetierL4]
# OK

bycatch[ , .N, keyby = MetierL5]
# OK

bycatch[ , .N, keyby = MetierL6]
# OK

bycatch[ , .N, keyby = VesselLengthRange]
# F?

bycatch[ , .N, keyby = MonitoringProgramType]
# OK

bycatch[ , .N, keyby = SamplingProtocol]
# ALL vs All vs all

bycatch[ , .N, keyby = MonitoringMethod]
# OK

bycatch[ , .N, keyby = Species]
# mixed case

bycatch[ , .N, keyby = CetaceanFlag]
# OK

bycatch[ , .N, keyby = IndividualsWithPingers]
# -9?

bycatch[ , .N, keyby = IndividualsWithoutPingers]
# -9?

bycatch[ , summary(IncidentsWithoutPingers)]
# OK

bycatch[ , summary(BycatchRateWithPingers)]
# OK

bycatch[ , summary(BycatchRateWithoutPingers)]
# OK?

bycatch[ , summary(TotalBycatchEstimate)]
# OK

bycatch[ , summary(CV)]
# OK

bycatch[ , .N, keyby = id]
# ?

bycatch[ , .N, keyby = tu_rank]
# OK

bycatch[ , .N, keyby = rank_name]
# OK

bycatch[ , .N, keyby = classname]
# OK

bycatch[ , .N, keyby = SuperClass]
# OK

bycatch[ , .N, keyby = Vernacular]
# OK


