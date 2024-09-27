
beam_download <- function(years, verbose = TRUE) {

    grab_and_save <- function(url, filter = FALSE, file, report = TRUE) {
        resp <- ices_get_jwt(url, username = "", jwt = token)
        dat <- fromJSON(content(resp, as = "text"))
        if (filter == TRUE) {
            dat <- subset(dat, year %in% years)
        }
        fwrite(dat, file = file, sep = ";")
        if (verbose == TRUE) {
            cat(sprintf("Saved %s\n", file))
        }
        as.data.table(dat)
    }

    # note: these vars are all defined in the global environment 
    all1 <<- grab_and_save(url = "https://bycatch.ices.dk/api/GetD1_Fishing_effort", 
                        filter = TRUE, 
                        file = "data/all1.csv")
    
    obs1 <<- grab_and_save(url = "https://bycatch.ices.dk/api/GetD2_Bycatch_monitoring_effort", 
                        filter = TRUE, 
                        file = "data/obs1.csv")
    
    bycatch1 <<- grab_and_save(url = "https://bycatch.ices.dk/api/GetD3_BycatchEvent", 
                        filter = TRUE, 
                        file = "data/bycatch1.csv")
    
    D4 <<- grab_and_save(url = sprintf("https://bycatch.ices.dk/api/GetOverviewSubmissionTable/%s", max(years)), 
                        filter = FALSE,
                        file = "data/D4_overviewsubmission.csv")
    
    D5 <<- grab_and_save(url = "https://bycatch.ices.dk/api/GetByCatchRoadMapListSpecies",
                        filter = FALSE, 
                        file = "data/D5_THELIST_roadmapspecies.csv")

}


# token<-"" #get manually via POST
# 
# #https://github.com/ices-eg/wg_WGBYC/tree/master
# 
# linkD1<-"https://bycatch.ices.dk/api/GetD1_Fishing_effort"
# linkD2<-"https://bycatch.ices.dk/api/GetD2_Bycatch_monitoring_effort"
# linkD3<-"https://bycatch.ices.dk/api/GetD3_BycatchEvent"
# linkD4<-"https://bycatch.ices.dk/api/GetOverviewSubmissionTable/2023"
# linkD5<-"https://bycatch.ices.dk/api/GetByCatchRoadMapListSpecies"
# 
# resp1<-ices_get_jwt(linkD1,username="",jwt=token)
# D1<-content(resp1,as="text")
# D1<-fromJSON(D1)
# resp2<-ices_get_jwt(linkD2,username="",jwt=token)
# D2<-content(resp2,as="text")
# D2<-fromJSON(D2)
# resp3<-ices_get_jwt(linkD3,username="",jwt=token)  
# D3<-content(resp3,as="text")
# D3<-fromJSON(D3)
# resp4<-ices_get_jwt(linkD4,username="",jwt=token)  
# D4<-content(resp4,as="text")
# D4<-fromJSON(D4)
# resp5<-ices_get_jwt(linkD5,username="",jwt=token)  
# D5<-content(resp5,as="text")
# D5<-fromJSON(D5)
# 
# write.csv(subset(D1,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D1_fishingeffort_2017_2023.csv")
# write.csv(subset(D2,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D2_monitoringeffort_2017_2023.csv")
# write.csv(subset(D3,year>2016),file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D3_bycatchevent_2017_2023.csv")
# write.csv(D4,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D4_overviewsubmission.csv")
# write.csv(D5,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/D5_THELIST_roadmapspecies.csv"