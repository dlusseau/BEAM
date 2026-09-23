
###
all2 <- copy(all1)
colnames(all2) <- tolower(colnames(all2))
all2[, country := clean_country_code(country)]
clean_chars(all2)
add_vessel_length_groups(all2)

all2 <- all2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]
all2[metierl6 == "-", metierl6 := NA]
# grab year before current year
#all2 <- all2[year == year(Sys.time())-1] # <- might want to add a setting for this.
                                         # or at least print a message to console
                                         # saying what year of data we're using

# Save complete data before subsetting, because we need it to compute coverage
# later.


### added 23 Sept 2026 - correction to be added post benchmark and agreed with countries at wkmamby
#Divide effort by two for paired operations in some countries
all2[  metierl4 %in% c("ptm", "ptb", "spr") &  !country %in% c("dnk", "nld"),  daysatseaf := daysatseaf / 2]



fwrite(all2, "data/fishing_all_years.csv", sep = ";")

#Grab specific year 
all2 <- all2[year == max(years)]

fwrite(all2, "data/all2.csv", sep = ";")



###
obs2 <- copy(obs1)
colnames(obs2) <- tolower(colnames(obs2))
obs2[, country := clean_country_code(country)]
clean_chars(obs2)
add_vessel_length_groups(obs2)

obs2 <- obs2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]
obs2 <- obs2[!(monitoringmethod == "lb" & country != "prt")]
obs2 <- obs2[!(monitoringmethod == "oth" & country != "nor")]
obs2 <- obs2[(monitoringmethod != "po")]
obs2 <- obs2[!(monitoringmethod == "vo" & country == "est")]
fwrite(obs2, "data/obs2.csv", sep = ";")

###
bycatch2 <- copy(bycatch1)
bycatch2[, country := clean_country_code(country)]
colnames(bycatch2) <- tolower(colnames(bycatch2))


### added 14 Sept
# Depending on how the data have been downloaded and read into R, we may have
# a situation where IncidentsWithPingers and IncidentsWithoutPingers are
# interpreted by R as character, because either of these columns may be "NULL"
# (i.e. a literal NULL string). So here we just check for that and correct it
# if needed. However, if you download BEAM data through the intended workflow,
# you'll not have that problem.
if (inherits(bycatch2$incidentswithpingers, "integer") == FALSE) {
    bycatch2[, incidentswithpingers := as.integer(incidentswithpingers)]
}
if (inherits(bycatch2$incidentswithoutpingers, "integer") == FALSE) {
    bycatch2[, incidentswithpingers := as.integer(incidentswithoutpingers)]
}

clean_chars(bycatch2)
add_vessel_length_groups(bycatch2)

bycatch2 <- bycatch2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]
bycatch2 <- bycatch2[!(monitoringmethod == "lb" & country != "prt")]
bycatch2 <- bycatch2[!(monitoringmethod == "oth" & country != "nor")]
bycatch2 <- bycatch2[(monitoringmethod != "po")]
bycatch2 <- bycatch2[!(monitoringmethod == "vo" & country == "est")]
bycatch2[ , class := classname]
bycatch2[superclass == "reptilia" & is.na(classname), class := "reptilia"]
bycatch2[classname %in% c("actinopteri", "holocephali", "myxini", "petromyzonti"), class := "fish"]

bycatch2[individualswithpingers < 0, individualswithpingers := NA] #minus values = NA
bycatch2[incidentswithpingers < 0, incidentswithpingers := NA] #minus values = NA

#######
bycatch2[ , n_individ := rowSums(.SD, na.rm = TRUE),.SDcols=c("individualswithpingers","individualswithoutpingers")] #ignore NAs when summing
bycatch2[ , n_incident := rowSums(.SD, na.rm = TRUE),.SDcols=c("incidentswithpingers","incidentswithoutpingers")] #ignore NAs when summing

fwrite(bycatch2, "data/bycatch2.csv", sep = ";",na = "NA")

