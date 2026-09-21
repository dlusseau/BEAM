
#' Downloads data needed for BEAM from ICES bycatch database (or loads a local
#' cached version of it if available) 
#' 
#' @description
#' `beam_get_raw_data` Downloads data needed for BEAM from ICES or loads local cached data
#' 
#' @param filename Filename, one of `obs1.csv`, `all1.csv`, `bycatch1.csv` or `D4_overviewsubmission.csv`
#' @param years integer vector, one or more years to use to subset data
#' @param force_download TRUE to force downloading. Defaults to FALSE.
#' @param verbose Logical. TRUE to report progress info during downloading. Defaults to TRUE.
#' @returns A data.table with the requested data
#' @export
#' 

beam_get_raw_data <- function(filename, years, force_download = FALSE, verbose = TRUE) {
    
    ices_api_urls <- list(
        D1 = "https://bycatch.ices.dk/api/GetD1_Fishing_effort?Year=%s",
        D2 = "https://bycatch.ices.dk/api/GetD2_Bycatch_monitoring_effort?Year=%s",
        D3 = "https://bycatch.ices.dk/api/GetD3_BycatchEvent?Year=%s",
        D4 = "https://bycatch.ices.dk/api/GetOverviewSubmissionTable/%d"
    )
    
    if (file.exists(filename) & force_download == FALSE) {
        return(fread(filename))
    }
    
    filenames <- c("all1.csv", "obs1.csv", "bycatch1.csv", "D4_overviewsubmission.csv")
    
    api <- which(basename(filename) == filenames)[1]
    
    if (is.na(api) | length(api) != 1) {
        stop("BEAM: Don't know how to get this file.")
    }
    
    # Append the year(s) we want to the URL string. This is done slightly differently
    # depending on which data we're getting 
    if (api %in% 1:3) {
        url <- sprintf(ices_api_urls[[api]], paste0(years, collapse = ","))
        info <- sprintf(" (years=%s-%s)", min(years), max(years))
    } else if (api == 4) {
        url <- sprintf(ices_api_urls[[api]], max(years)+1)
        info <- sprintf(" (year=%d)", max(years)+1)
    }
    
    if (verbose == TRUE) {
        cat(sprintf("BEAM: Downloading %s %s%s", 
                    names(ices_api_urls)[api], filename, info))
    }
    
    t_start <- Sys.time()
    
    if (packageVersion("icesConnect")<"1.2") {
    resp <- icesConnect::ices_get_jwt(url, username = ices_username, jwt = ices_token, quiet = TRUE)
    } else {
    resp <- icesConnect::ices_get_jwt(url, jwt = ices_token, quiet = TRUE)  
    }
    dat <- jsonlite::fromJSON(httr::content(resp, as = "text"))
    fwrite(dat, file = filename, sep = ";")
    
    t_end <- Sys.time()
    t_elapsed <- as.integer(round(difftime(t_end,t_start,units="secs"),0))
    
    if (verbose == TRUE) {
        cat(sprintf("\rBEAM: Downloading %s complete%s in %ds, saved as %s\n", 
                    names(ices_api_urls)[api], info, t_elapsed, filename))
    }
    
    as.data.table(dat)
}
