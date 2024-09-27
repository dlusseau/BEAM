
###
all2 <- copy(all1)
colnames(all2) <- tolower(colnames(all2))
clean_chars(all2)
add_vessel_length_groups(all2)

all2 <- all2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]

# grab year before current year
all2 <- all2[year == year(Sys.time())-1] # <- might want to add a setting for this.
                                         # or at least print a message to console
                                         # saying what year of data we're using
all2$metierl6[all2$metierl6 == "-"] <- NA



fwrite(all2, "data/all2.csv", sep = ";")

###
obs2 <- copy(obs1)
colnames(obs2) <- tolower(colnames(obs2))
clean_chars(obs2)
add_vessel_length_groups(obs2)

obs2 <- obs2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]
obs2 <- obs2[!(monitoringmethod == "lb" & country != "pt")]
obs2 <- obs2[!(monitoringmethod == "oth" & country != "no")]
obs2 <- obs2[(monitoringmethod != "po")]
obs2 <- obs2[!(monitoringmethod == "vo" & country == "ee")]
fwrite(obs2, "data/obs2.csv", sep = ";")

###
bycatch2 <- copy(bycatch1)
colnames(bycatch2) <- tolower(colnames(bycatch2))
clean_chars(bycatch2)
add_vessel_length_groups(bycatch2)

bycatch2 <- bycatch2[!is.na(ecoregion) & !(ecoregion %in% c("north west atlantic"))]
bycatch2 <- bycatch2[!(monitoringmethod == "lb" & country != "pt")]
bycatch2 <- bycatch2[!(monitoringmethod == "oth" & country != "no")]
bycatch2 <- bycatch2[(monitoringmethod != "po")]
bycatch2 <- bycatch2[!(monitoringmethod == "vo" & country == "ee")]
bycatch2[ , class := classname]
bycatch2[superclass == "reptilia" & is.na(classname), class := "reptilia"]
bycatch2[classname %in% c("actinopteri", "holocephali", "myxini", "petromyzonti"), class := "fish"]
bycatch2[ , n_individ := individualswithpingers + individualswithoutpingers]
bycatch2[ , n_incident := incidentswithpingers + incidentswithoutpingers]
fwrite(bycatch2, "data/bycatch2.csv", sep = ";")

