###########################################################################################
###########################################################################################
##### BEAM pipeline
##### Sept 2024
###########################################################################################
###########################################################################################

###########################################################################################
###########################################################################################

#https://github.com/dlusseau/BEAM
##### first let's get the data 
require(icesConnect)
require(httr)
require(jsonlite)

#token manually via postman can be done with httr::POST  simply too
#postman: https://www.postman.com/ 
#POST: https://bycatch.ices.dk/api/token
# send
# Body: copy token value
#Authorization tab> Auth type = Bearer token
# paste token in token textbox (make sure to remove "")

#GET https://bycatch.ices.dk/api/GetD2_Bycatch_monitoring_effort?Year=2021

token<-"" #get manually via POST

#https://github.com/ices-eg/wg_WGBYC/tree/master

linkD1<-"https://bycatch.ices.dk/api/GetD1_Fishing_effort"
linkD2<-"https://bycatch.ices.dk/api/GetD2_Bycatch_monitoring_effort"
linkD3<-"https://bycatch.ices.dk/api/GetD3_BycatchEvent"
linkD4<-"https://bycatch.ices.dk/api/GetOverviewSubmissionTable/2023"
linkD5<-"https://bycatch.ices.dk/api/GetByCatchRoadMapListSpecies"

resp1<-ices_get_jwt(linkD1,username="",jwt=token)
D1<-content(resp1,as="text")
D1<-fromJSON(D1)
resp2<-ices_get_jwt(linkD2,username="",jwt=token)
D2<-content(resp2,as="text")
D2<-fromJSON(D2)
resp3<-ices_get_jwt(linkD3,username="",jwt=token)  
D3<-content(resp3,as="text")
D3<-fromJSON(D3)
resp4<-ices_get_jwt(linkD4,username="",jwt=token)  
D4<-content(resp4,as="text")
D4<-fromJSON(D4)
resp5<-ices_get_jwt(linkD5,username="",jwt=token)  
D5<-content(resp5,as="text")
D5<-fromJSON(D5)

write.csv(subset(D1,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D1_fishingeffort_2017_2023.csv")
write.csv(subset(D2,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D2_monitoringeffort_2017_2023.csv")
write.csv(subset(D3,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D3_bycatchevent_2017_2023.csv")
write.csv(D4,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D4_overviewsubmission.csv")
write.csv(D5,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D5_THELIST_roadmapspecies.csv")



##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
### now QC
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

### using code from Henrik

require(data.table)

D1<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D1_fishingeffort_2017_2023.csv"
D2<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D2_monitoringeffort_2017_2023.csv"
D3<-"C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D3_bycatchevent_2017_2023.csv"

##########################################################################################################################################################################################################
# prepare fishing effort ------

fish_effort = fread(D1,na.strings = "NULL")

dim(fish_effort)
str(fish_effort)
names(fish_effort)

# check variables

fish_effort[ , .N, keyby = ecoregion]

# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions
# 2024: 
# Exclude 'North West Atlantic' #NAFO area in the Labrador

nrow(fish_effort)

fish_effort = fish_effort[!is.na(ecoregion) &
                            !(ecoregion %in% c("North West Atlantic"))]
							
nrow(fish_effort)
# [1] 127831

##########################################################################################################################################################################################################
# prepare monitoring effort ------

monitor = fread(D2,na.strings = "NULL")

dim(monitor)


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

table(monitor$monitoringMethod)

# remove LB (logbook) except for Portugal
monitor = monitor[!(monitoringMethod == "LB" & country != "PT")]
nrow(monitor)

# remove OTH (other) except for Norway
monitor = monitor[!(monitoringMethod == "OTH" & country != "NO")]
nrow(monitor)


# remove PO (port observers)
monitor = monitor[(monitoringMethod != "PO")]
nrow(monitor)


# remove VO (vessel observers) for Estonia
monitor = monitor[!(monitoringMethod == "VO" & country == "EE")]
nrow(monitor)



# check Ecoregion -----
monitor[ , .N, keyby = ecoregion]

# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions
#2024
# Exclude missing, 'North West Atlantic' 

monitor = monitor[!is.na(ecoregion) &
                    !(ecoregion %in% c("North West Atlantic"))]

any(unique(monitor$ecoregion)%in%unique(fish_effort$ecoregion)==FALSE) #mismatch logical check for ecoregion names
any(unique(fish_effort$ecoregion)%in%unique(monitor$ecoregion)==FALSE) #mismatch logical check for ecoregion names

######################################################
# if TRUE then problem - buil function with a catch
######################################################

##########################################################################################################################################################################################################
# prepare bycatch ------

bycatch = fread(D3,na.strings = c("NULL", -9))

dim(bycatch)

str(bycatch)


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


# remove LB (logbook) except for Portugal
bycatch = bycatch[!(monitoringMethod == "LB" & country != "PT")]
nrow(bycatch)


# remove OTH (other) except for Norway
bycatch = bycatch[!(monitoringMethod == "OTH" & country != "NO")]
nrow(bycatch)


# remove PO (port observers)
bycatch = bycatch[(monitoringMethod != "PO")]
nrow(bycatch)


# remove VO (vessel observers) for Estonia
bycatch = bycatch[!(monitoringMethod == "VO" & country == "EE")]
nrow(bycatch)


# check Ecoregion -----

bycatch[ , .N, keyby = ecoregion]

# compare with
# https://www.ices.dk/data/Documents/Maps/ICES-Ecoregions-hybrid-statistical-areas.png

# correct Ecoregions

# Exclude missing, 'North West Atlantic' 

bycatch = bycatch[!is.na(ecoregion) &
                    !(ecoregion %in% c("North West Atlantic"))]

nrow(bycatch)

any(unique(bycatch$ecoregion)%in%unique(monitor$ecoregion)==FALSE) #mismatch logical check for ecoregion names


# create new class in bycatch -------

# Katja:
# Bird
# Elasmo
# Reptile: Superclass == "Reptilia" & is.na(classname)
# Fish: Collapse Actinopteri + Holocephali + Myxini + Petromyzonti to fish

# first, copy classname to new variable class
bycatch[ , class := classname]

bycatch[superClass == "Reptilia" & is.na(classname), class := "Reptilia"]

bycatch[classname %in% c("Actinopteri", "Holocephali", "Myxini", "Petromyzonti"),
         class := "Fish"]


# sum pinger and non-pinger for individuals and incidents ------
###
table(as.numeric(bycatch$individualsWithPingers))
### we have NAs

# number of individuals
bycatch[ , n_individ := as.numeric(individualsWithPingers) + as.numeric(individualsWithoutPingers)]

# number of incidents
bycatch[ , n_incident := as.numeric(incidentsWithPingers) + as.numeric(incidentsWithoutPingers)]



##########################################################################################################################################################################################################
# effort data from RDB ------
#######################????
# fish_effort_rdb = fread("data/Effort data CE from RDB 2017 2021.csv",
                   # na.strings = "NULL")
# dim(fish_effort_rdb)
# # [1] 603498     17
# # EcoRegion needs to added from lookup table "Area" - EcoRegion
# # ICES_Area_to_EcoRegion.csv
# area_ecoregion = fread("data/ICES_Area_to_EcoRegion.csv")
# fish_effort_rdb[area_ecoregion, on = c(Area = "ICES_Area"),
           # EcoRegion := i.Ecoregion]
# # several non-matched rows
# fish_effort_rdb[is.na(EcoRegion), .N]
# # [1] 5632
# # unmatched Areas
# fish_effort_rdb[is.na(EcoRegion), sort(unique(Area))]
# # [1] "21.1.A"   "21.1.B"   "21.1.C"   "21.1.D"  
# # [5] "21.4.V"   "21.6.G"   "21.6.H"   "34"      
# # [9] "34.1.1"   "34.1.1.3" "34.1.2"   "34.1.3"  
# # [13] "34.1.3.1" "34.1.3.2" "34.2"     "34.3.1"  
# # [17] "34.3.1.1" "34.3.1.2" "34.3.1.3" "34.3.3"  
# # [21] "34.3.6"   "47.1"     "47.1.1"   "47.1.2"  
# # [25] "47.1.3"   "47.1.4"   "47.A.1"   "87"      
# # [29] "87.2.6"   "87.3.3"
# # Katja:
# # 27.2: set EcoRegion to "Barents Sea"
# # leave rest as NA
# # remove NA, North West Atlantic, Med. Sea
# fish_effort_rdb[Area == "27.2", EcoRegion := "Barents Sea"]
# fish_effort_rdb = fish_effort_rdb[!is.na(EcoRegion) & 
                                    # !(EcoRegion %in% c("North West Atlantic", "Mediterranean Sea"))]
# # metier 4 needs to be created from "FishingActivityCategoryEuropeanLvl5"
# # Extract string before "_"
# fish_effort_rdb[ , MetierL4 := sub("_.*", "", FishingActivityCategoryEuropeanLvl5)]
# nrow(fish_effort_rdb)
# # [1] 597163
# # select relevant rows and columns









# check variables ------
###we should output this to a file
# check effort --------
fish_effort[ , .N, tblUploadID]
#OK
fish_effort[ , .N, keyby = year]
# OK
fish_effort[ , .N, keyby = month]
# 9835 NA
fish_effort$month<-as.numeric(fish_effort$month)

fish_effort[ , .N, keyby = quarter]
# 715 NA
fish_effort$quarter<-as.numeric(fish_effort$quarter)

fish_effort[ , .N, keyby = areaType]
# OK

fish_effort[ , .N, keyby = areaCode]
# OK

fish_effort[ , .N, keyby = metierL3]
# OK

fish_effort[ , .N, keyby = metierL4]
# OK

fish_effort[ , .N, keyby = metierL5]
# OK

fish_effort[ , .N, keyby = metierL6]
# <NA> zero
# -  1141
# set "-" to NA
fish_effort$metierL6[fish_effort$metierL6=="-"]<-NA

fish_effort[ , .N, keyby = vesselLengthRange]
# OK levels
# Set Unknown to NA?


fish_effort[ , summary(daysAtSeaF)]
# OK?

fish_effort[ , .N, keyby = vesselsF]
# strange values
# read as integer64
# what kind of variable is this?
fish_effort$vesselsF<-as.integer(fish_effort$vesselsF)

fish_effort[ , summary(vesselsF)]

fish_effort$tripsF<-as.numeric(fish_effort$tripsF)
fish_effort[ , summary(tripsF)]

# OK?
# veery high max

#fish_effort[ , all.equal(tripsF, tripsF_txt)]
# equal to TripsF

fish_effort$totalNetsLengthF<-as.numeric(fish_effort$totalNetsLengthF)
fish_effort[ , summary(totalNetsLengthF)]
# OK?
# 200613 NA

fish_effort$totalKmHoursF<-as.numeric(fish_effort$totalKmHoursF)
fish_effort[ , summary(totalKmHoursF)]
# OK?
# 200613 NA

fish_effort$noofHaulsF<-as.numeric(fish_effort$noofHaulsF)
fish_effort[ , summary(noofHaulsF)]
# 192574 NA
# very high max

fish_effort$totalTowTimeF<-as.numeric(fish_effort$totalTowTimeF)
fish_effort[ , summary(totalTowTimeF)]
# very high max

fish_effort[ , .N, keyby = country]
# OK
##PT-20
fish_effort$ecoregion[fish_effort$country=="PT-20",]
#correct it is Azores autonomous region iso country code

fish_effort[ , .N, keyby = organization]


# check monitor ----------

monitor[ , .N, keyby = submissionYear]
# OK

monitor[ , .N, keyby = tblUploadID]
# OK

monitor[ , .N, keyby = country]
# OK  

monitor[ , .N, keyby = organization]
# OK

monitor[ , .N, keyby = year]
# OK

monitor[ , .N, keyby = month]
# OK

monitor[ , .N, keyby = areaType]
# ICESarea
monitor$areaType[monitor$areaType=="ICESarea"]<-"ICESArea"

monitor[ , .N, keyby = areaCode]
# OK
#### !!!!!!!!!!!!!!!!!! we have levels containing tildas "~" for bycatch which cannot be assign to only one area as we move to area level estimation we need a rule to deal with those

monitor[ , .N, keyby = metierL3]
# OK

monitor[ , .N, keyby = metierL4]
# OK

monitor[ , .N, keyby = metierL5]
# OK

monitor[ , .N, keyby = metierL6]
# !!!!!!!!!!!!!!!!!NO NAs!!!!!
table(monitor$metierL6)

monitor[ , .N, keyby = vesselLengthRange]
# F?

monitor[ , .N, keyby = monitoringProgramType]
# OK

monitor[ , .N, keyby = samplingProtocol]
# OK

monitor$samplingProtocol[monitor$samplingProtocol=="seabirds"]<-"Seabirds"
monitor$samplingProtocol[monitor$samplingProtocol=="All"]<-"ALL"
monitor$samplingProtocol[monitor$samplingProtocol=="all"]<-"ALL"


monitor[ , .N, keyby = monitoringMethod]
# OK

monitor[ , .N, keyby = targetSpecies]

#multiple species possible sometimes species separated by tildas but sometimes just a space make sure to look for species rather than lok for exact match to latin name

# na.string: #N/A

monitor[ , .N, keyby = pingerCharacteristics]
# OK

monitor[ , .N, keyby = othMitigationMeasures]
# OK 
str(monitor)

table(monitor$vesselsOb)
monitor$vesselsOb<-as.integer(monitor$vesselsOb)
monitor[ , summary(vesselsOb)]
# OK

table(monitor$tripsOb)
monitor$tripsOb<-as.numeric(monitor$tripsOb)
monitor[ , summary(tripsOb)]
# Very high max!

monitor[ , summary(daysAtSeaOb)]
# Very high max!

table(monitor$haulsWithPingersOb)
monitor$haulsWithPingersOb<-as.numeric(monitor$haulsWithPingersOb)
monitor[ , summary(haulsWithPingersOb)]
# OK
monitor$totalNetsLengthOb<-as.numeric(monitor$totalNetsLengthOb)
monitor[ , summary(totalNetsLengthOb)]
# very high max!

monitor$totalKmHoursOb<-as.numeric(monitor$totalKmHoursOb)
monitor[ , summary(totalKmHoursOb)]
# very high max

monitor$noofHaulsOb<-as.numeric(monitor$noofHaulsOb)
monitor[ , summary(noofHaulsOb)]
# OK

monitor$totalTowTimeOb<-as.numeric(monitor$totalTowTimeOb)
monitor[ , summary(totalTowTimeOb)]
# very high max

monitor[ , .N, keyby = monitoring812Type]
# OK


# check bycatch ----------
bycatch[ , .N, keyby = submissionYear]
# OK

bycatch[ , .N, keyby = tblUploadID]
# OK

bycatch[ , .N, keyby = country]
# OK

bycatch[ , .N, keyby = organization]
# OK

bycatch[ , .N, keyby = recordType]
# OK

bycatch[ , .N, keyby = year]
# OK

bycatch[ , .N, keyby = month]
# OK

bycatch[ , .N, keyby = areaType]
bycatch$areaType[bycatch$areaType=="ICESarea"]<-"ICESArea"
bycatch$areaType[bycatch$areaType=="IcesArea"]<-"ICESArea"


bycatch[ , .N, keyby = areaCode]
# OK?
# e.g. 1~5~6

bycatch[ , .N, keyby = metierL3]
# OK

bycatch[ , .N, keyby = metierL4]
# OK

bycatch[ , .N, keyby = metierL5]
# OK

bycatch[ , .N, keyby = metierL6]
# OK

bycatch[ , .N, keyby = vesselLengthRange]
# F?

bycatch[ , .N, keyby = monitoringProgramType]
# OK

bycatch[ , .N, keyby = samplingProtocol]
# ALL vs All vs all
table(monitor$samplingProtocol)
bycatch$samplingProtocol[bycatch$samplingProtocol=="all"]<-"ALL"
bycatch$samplingProtocol[bycatch$samplingProtocol=="All"]<-"ALL"


bycatch[ , .N, keyby = monitoringMethod]
# OK

bycatch[ , .N, keyby = species]
# mixed case
table(bycatch$species)


bycatch[ , .N, keyby = cetaceanFlag]
# OK

bycatch[ , .N, keyby = individualsWithPingers]
# -9? NO MORE
table(bycatch$individualsWithPingers)
bycatch$individualsWithPingers<-as.numeric(bycatch$individualsWithPingers)


bycatch[ , .N, keyby = individualsWithoutPingers]
# -9? NO MORE
table(bycatch$individualsWithoutPingers)
bycatch$individualsWithoutPingers<-as.numeric(bycatch$individualsWithoutPingers)

bycatch[ , summary(incidentsWithoutPingers)]
# OK

bycatch[ , summary(bycatchRateWithPingers)]
# OK
table(bycatch$bycatchRateWithPingers)
bycatch$bycatchRateWithPingers<-as.numeric(bycatch$bycatchRateWithPingers)

bycatch[ , summary(bycatchRateWithoutPingers)]
# OK?
table(bycatch$bycatchRateWithoutPingers)
bycatch$bycatchRateWithoutPingers<-as.numeric(bycatch$bycatchRateWithoutPingers)


bycatch[ , summary(totalBycatchEstimate)]
# OK
table(bycatch$totalBycatchEstimate)
bycatch$totalBycatchEstimate<-as.numeric(bycatch$totalBycatchEstimate)


bycatch[ , summary(cv)]
# OK
table(bycatch$cv)
bycatch$cv<-as.numeric(bycatch$cv)

bycatch[ , .N, keyby = id]
# ?

bycatch[ , .N, keyby = tu_rank]
bycatch$tu_rank<-as.numeric(bycatch$tu_rank)

# OK

bycatch[ , .N, keyby = rank_name]
# OK

bycatch[ , .N, keyby = classname]
# OK

bycatch[ , .N, keyby = superClass]
# OK

bycatch[ , .N, keyby = vernacular]
# OK


##########################################################################################################################################################################################################
# save files as csv -----

#fwrite(fish_effort_rdb, "data/clean/fish_effort_rdb.csv", sep = ";")   ####check
fwrite(fish_effort, "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/fish_effort_wgbyc_2017_2024_clean.csv", sep = ";")
fwrite(monitor, "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/monitor_wgbyc_2017_2024_clean.csv", sep = ";")
fwrite(bycatch, "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/bycatch_wgbyc_2017_2024_clean.csv", sep = ";")



##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
### now prepare the data
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################

#make tables on (1) bycatch monitoring effort, (2) number of bycatch

# Observational unit
#Ecoregion, Country, year, species,
# metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol, 

require(data.table)



##########################################################################################################################################################################################################
# prepare monitoring effort ------
monitor = fread("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/monitor_wgbyc_2017_2024_clean.csv",
                  na.strings = "NA")

dim(monitor)
#checks out

str(monitor)
#better

monitor[ , .N, keyby = vesselLengthRange]

# two strange levels, "F" and "NK"
# set to missing
monitor[vesselLengthRange %in% c("F", "NK"),
        vesselLengthRange := NA]
 
# create vessel length group: below and above 12 meter
# there are two length categories that overlap 12:
# : VL0015, VL0815
# for now, set those to below 12

monitor[vesselLengthRange %in% c(
  "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
  vesselLength_group := "below_12"]


monitor[vesselLengthRange %in% c(
  "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
  vesselLength_group := "above_12"]

monitor[ , .N, keyby = vesselLength_group]
# OK

# sum bycatch monitoring effort ------------------

#Observational unit, group by:
#Ecoregion, areaCode, Country, year,
#metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol,



monitor_sum = monitor[
  ,
  .(daysAtSea = sum(daysAtSeaOb, na.rm = TRUE)),
  keyby = .(ecoregion, areaCode,country, year,
            metierL4, metierL5, vesselLength_group,
            samplingProtocol, monitoringMethod)]

###################################
# For rows with monitoring effort, but no bycatch, add rows with zero bycatch
###################################


# prepare bycatch ------

bycatch = fread(
  "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/bycatch_wgbyc_2017_2024_clean.csv",
  na.strings = c("NULL", -9,"NA"))


dim(bycatch)
# [1] 9126   33

str(bycatch)




# Filter bycatch data on MonitoringMethod & Country -------


# check and create variables -----

bycatch[ , .N, keyby = vesselLengthRange]

# two strange levels, "F" and "NK"
# set to missing
bycatch[vesselLengthRange %in% c("F", "NK"),
        vesselLengthRange := NA]

# create vessel length group: below and above 12 meter
# there are two length categories that overlap 12:
# : VL0015, VL0815
# for now, set those to below 12

bycatch[vesselLengthRange %in% c(
  "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
  vesselLength_group := "below_12"]


bycatch[vesselLengthRange %in% c(
  "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
  vesselLength_group := "above_12"]

bycatch[ , .N, keyby = vesselLength_group]
# OK


# new species variable with names tolower
bycatch[ , species := tolower(species)]


# number of individuals catch




# sum bycatch  ------------------

#Observational unit, group by:
#Ecoregion, areaCode,Country, year,
#metier 4, metier 5,  vessel length category (below or above 12m) Bycatch monitoring method, bycatch sampling protocol, species

#Value variables: DaysAtSeaOb bycatch monitor effort DaS, number of individuals

bycatch_sum = bycatch[
  ,
  .(n_ind = sum(n_individ, na.rm = TRUE)),
  keyby = .(ecoregion, areaCode,country, year,
            metierL4, metierL5, vesselLength_group,
            samplingProtocol, monitoringMethod, species)]




# create list of relevant Ecoregion * species combinations  #########


# read species list

# file recieved from David 2023-09-20 15:51
annex01_species = read.csv(
  "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D5_THELIST_roadmapspecies.csv",row.names=1,header=T)
  
mediterranean<-read.csv(
  "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/metadata/MediterraneanPrioritySpeciesToBeIncluded.csv",sep=";")
  
leftover2023<-read.csv(
  "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/metadata/caps_species_needed.csv")
  
##### now we have all the species x ecoregion set we had at the end of the show in 2023, let's merge and assure that we do not have duplicates  

  annex01_species<-as.data.table(annex01_species)
  mediterranean<-as.data.table(mediterranean)
 
 
# select and rename relevant column

annex01_species = annex01_species[ , .(ecoRegion, speciesName)]

str(annex01_species)
annex01_species[ , .N, keyby = ecoRegion]

colnames(annex01_species)<-c("Ecoregion","Species")

str(mediterranean)
mediterranean[ , .N, keyby = Ecoregion]


ecoreg_species = rbindlist(
  list(annex01_species, mediterranean)
)

str(ecoreg_species)


ecoreg_species[ , species := tolower(Species)]

###
(leftover2023$x)%in%ecoreg_species$species

###they are now all there ... let's check again after the monitor intersection
 colnames(ecoreg_species)[1]<-"ecoregion"
 
 ###problem with ecoregion levels
 ecoreg_species$ecoregion<-gsub(" ecoregion", "",ecoreg_species$ecoregion)
 
 "Bay of Biscay and the Iberian Coast"
 unique(ecoreg_species$ecoregion)
  unique(ecoreg_species$ecoregion)%in%unique(monitor_sum$ecoregion)

ecoreg_species$ecoregion[ecoreg_species$ecoregion=="Bay of Biscay and Iberian Coast"]<-"Bay of Biscay and the Iberian Coast"
 
  unique(ecoreg_species$ecoregion)%in%unique(monitor_sum$ecoregion)

##################################
### 2024 we do not have monitoring in the Arctic Ocean  
##################################


# to obtain one row per (relevant) species and ecoregion
# join monitor_sum with Ecoregion * species
monitor_sum<-(subset(monitor_sum,!is.na(ecoregion)))

monitor_ecoreg_species = monitor_sum[
  ecoreg_species,
  on = .(ecoregion),
  allow.cartesian = TRUE]
  
##should get about 267460 

nrow(monitor_sum)
# [1] 1752
nrow(ecoreg_species)
# [1] 1238
nrow(monitor_ecoreg_species)
# [1] 184419

str(monitor_ecoreg_species)
monitor_ecoreg_species[ , .N, keyby = ecoregion]
monitor_ecoreg_species[ , .N, keyby = country]
monitor_ecoreg_species[ , .N, keyby = metierL4]
monitor_ecoreg_species[ , .N, keyby = metierL5]
monitor_ecoreg_species[ , .N, keyby = vesselLength_group]
monitor_ecoreg_species[ , .N, keyby = samplingProtocol]
monitor_ecoreg_species[ , .N, keyby = monitoringMethod]
monitor_ecoreg_species[ , summary(daysAtSea)]
monitor_ecoreg_species[ , .N, species]

table_monitor_eccoregion_2017_2023<-as.data.frame(unlist(by(monitor_ecoreg_species$daysAtSea,monitor_ecoreg_species$ecoregion,sum,simplify=TRUE)))
colnames(table_monitor_eccoregion_2017_2023)<-c("daysAtSea_2017_2023")
table_monitor_eccoregion_2017_2023$daysAtSea_2017_2023[is.na(table_monitor_eccoregion_2017_2023$daysAtSea_2017_2023)]<-0
table_monitor_eccoregion_2017_2023$ecoregion<-rownames(table_monitor_eccoregion_2017_2023)

write.csv(table_monitor_eccoregion_2017_2023,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/results/table0_DaS_ecoregion_2017_2023.csv")
# we have a record that there is no monitoring in the Arctic


monitor_ecoreg_species$ecoregion[which(is.na(monitor_ecoreg_species$country))]
monitor_sum$ecoregion[which(is.na(monitor_sum$country))]
unique(monitor_sum$ecoregion)



# create lower caps species
monitor_ecoreg_species[ , species := tolower(Species)]


sum(unique(monitor_ecoreg_species$species)%in%unique(bycatch_sum$species)==FALSE)
#not all species on THE ONE LIST TO RULE THEM ALL are in the bycatch record - 81 missing

sum(unique(bycatch_sum$species)%in%unique(monitor_ecoreg_species$species)==FALSE)
#127 not on THE LIST

#let's make sure they are not just typos or mispelled
library(stringdist)
themissing<-unique(bycatch_sum$species)[which(unique(bycatch_sum$species)%in%unique(monitor_ecoreg_species$species)==FALSE)]

the_list_species<-unique(monitor_ecoreg_species$species)

## what distance makes sense?
discovery.df<-data.frame(distance=c(1:max(nchar(themissing))),NAs=0)
for (i in 1:nrow(discovery.df)) {
closematch_missing<-the_list_species[amatch(themissing, the_list_species, maxDist = discovery.df$distance[i])]
discovery.df$NAs[i]<-sum(is.na(closematch_missing))
}

plot(NAs~distance,data=discovery.df)

missing_match<-data.frame(missing=themissing,match=the_list_species[amatch(themissing, the_list_species, maxDist = 7)])
#manually rationalising only keeping those that are indeed close (homonyms or typo)

missing_match<-missing_match[c(32,96),]


##common name introduced by accident in monitor somewhere
#Stellate sturgeon
#acipenser stellatus

monitor_ecoreg_species$species[which(monitor_ecoreg_species$species=="stellate sturgeon")]<-"acipenser stellatus"

# join number of individuals bycatch from bycatch_sum
monitor_ecoreg_species[bycatch_sum,
  on = .(ecoregion, areaCode,country, year,
  metierL4, metierL5, vesselLength_group,
  samplingProtocol, monitoringMethod,
  species),
  n_ind := i.n_ind]

which(monitor_ecoreg_species$species=="nesiarchus nasutus")
sum(which(is.na(monitor_ecoreg_species$country)))

# replace non-matches (NA rows) with zero
monitor_ecoreg_species[
  is.na(n_ind),
  n_ind := 0]

str(monitor_ecoreg_species)

##remove arctic ocean which has no monitoring

monitor_ecoreg_species<-subset(monitor_ecoreg_species,!is.na(monitor_ecoreg_species$country))


# check bycatch taxa that are not in 'SamplingProtocol' taxa ###############

# # Sara created a species-'sampled taxa' table. Sara: Assigns all 'Cetacea' to 'mammal' -> do the same in 'SamplingProtocol' 
# species_taxa = fread(
  # "A:/My Documents/Projects/Bycatch/Projects/WGBYC/meeting/2023/ToR C/data/species list/candidate BPUE greens.csv")

# species_taxa = species_taxa[ , .(Species, taxa)]

# # set species tolower
# species_taxa[ , species := tolower(Species)]

# species_taxa[ , .N, keyby = taxa]
# #      taxa  N
# # 1:   Bird 31
# # 2:   Fish 83
# # 3: Mammal 12
# # 4: Turtle  2

# monitor_ecoreg_species[ , .N, keyby = SamplingProtocol]
# #    SamplingProtocol      N
# # 1:             <NA>     37
# # 2:              All 138773
# # 3:        Cetaceans   2019
# # 4:             Fish  21029
# # 5:          Mammals    965
# # 6: ProtectedSpecies  17795
# # 7:         Seabirds    765
# # 8:          Turtles   1300

# # rename levels in species_taxa to match those in monitor_ecoreg
# species_taxa[taxa == "Bird", taxa := "Seabirds"]
# species_taxa[taxa == "Mammal", taxa := "Mammals"]
# species_taxa[taxa == "Turtle", taxa := "Turtles"]


# create a a copy of SamplingProtocol, taxa_monitored
monitor_ecoreg_species[ , taxa_monitored := samplingProtocol]

unique(monitor_ecoreg_species$samplingProtocol)


# Sara assigned all 'Cetacea' to 'mammal' in species-taxa table -> do the same in 'SamplingProtocol' 
monitor_ecoreg_species[samplingProtocol == "Cetaceans",
                       taxa_monitored := "Mammals"]


#get to taxonomic cluster
devtools::install_github("ropensci/taxizedb")

library(taxizedb)
library(reshape2)

library(taxize)
sp.class<-tax_name(unique(monitor_ecoreg_species$species), get = "class")
sp.subclass<-tax_name(unique(monitor_ecoreg_species$species), get = "subclass")
sum(is.na(sp.subclass$subclass))
species.df<-sp.class
species.df$subclass<-sp.subclass$subclass
species.df$taxon<-"Fish"
species.df$taxon[species.df$class=="Aves"]<-"Seabirds"
species.df$taxon[species.df$class=="Mammalia"]<-"Mammals"
species.df$taxon[species.df$class=="Reptilia"]<-"Turtles"

write.csv(species.df,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/metadata/species2samplingProtocol.csv")

colnames(species.df)[2]<-"species"
species_taxa<-species.df[,c(2,5)]

# join taxa to monitor_ecoreg_species
monitor_ecoreg_species[species_taxa, on = .(species),
                       taxa := i.taxon]

monitor_ecoreg_species[ , .N, keyby = taxa]

# create variable indicating if a row is:
# taxa_monitored is 'All' OR 'ProtectedSpecies'
# taxa == taxa_monitored
monitor_ecoreg_species[
  ,
  taxa_bycatch_monitor_ok := (taxa_monitored %in% c("ALL","Elasmobranchs~Seabirds~Mammals", "ProtectedSpecies")) |
    (taxa_monitored == taxa)
]

monitor_ecoreg_species[ , .N, keyby = taxa_bycatch_monitor_ok]


str(monitor_ecoreg_species)

fwrite(
  monitor_ecoreg_species,
  "C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/monitor_effort_bycatch_ecoregion_area_species.csv")

##########################################################################################################################################################################################################
##########################################################################################################################################################################################################
##########################################################################################################################################################################################################






