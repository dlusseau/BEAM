
annotate_bpue <- function(bpue, cols) {
    #Reading total fishing effort to later compute monitoring coverage
    all_fishing_effort <- fread("data/fishing_all_years.csv")
    
    #Aggregating daysatsea across same levels used in calc_bpue
    all_fishing_effort <- all_fishing_effort[, .(
        daysatseaf = sum(daysatseaf, na.rm = TRUE)), by = cols]
    
    #Retrieve n_ind and daysatsea information (model inputs) from obs3, aggregated at the same level used in calc_bpue
    annotations <- obs3[taxon_bycatch_monitor_ok==TRUE, .(
        daysatsea = sum(daysatsea, na.rm = TRUE), n_ind = sum(n_ind, na.rm=T)), by =c("species", cols)]
    
    #Joining annotations to original bpue model
    
    bpue2 <- copy(bpue)

    bpue2[annotations, on = c("species", cols), c("daysatsea", "n_ind") := list(i.daysatsea, i.n_ind)]

     bpue2[is.na(daysatsea), daysatsea := 0] #Replace NAs by zero. There may be daysatsea in obs3 but not suitable for that taxa, hence the NAs introduced in the steps above. The number of daysatsea in bpue2 may be different than the one is obs3, because it only keeps days with suitable sampling protocol for the species taxa.
     bpue2[is.na(n_ind), n_ind := 0] #Replace NAs by zero. It means here there was no individual caught reported, when taxon_bycatch_monitor_ok = TRUE. There may be some in obs3, but they were not used in bpue estimation.
    
    bpue2[all_fishing_effort, on = cols, daysatseaf := i.daysatseaf]
    
    bpue2[,monitoring_coverage := daysatsea / daysatseaf]
    bpue2
}







