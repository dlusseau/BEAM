
#' Calculates total bycatch for bycatch models fitted in BEAM
#' 
#' @description
#' `calc_total` calculates total bycatch based on a BPUE.
#' 
#' #' @details
#' If there are multiple rows in needle, each row is processed separately, and
#' combined into one data.table with nrow equal to the number of rows in needle.
#' 
#' @param bpue data.table with values for cols, used to subset observations from mon. Any additional columns in needle are disregarded.
#' @param cols columns over which to split the data given in obs and all. Note that if `species` is among col, it is not used when subsetting `all`.
#' @param obs data.table with observed (monitored) fishing effort data (e.g. for the last 5 years)
#' @param all data.table with all (unobserved) fishing effort data, but only for the current year
#' @param verbose Logical, passed on to predict_response. Set to FALSE if nrow(bpue) > 1.
#' @param years vector with integers indicating years of assessment. 
#' @param include.weights boolean indicator indicating whether observations from the five most recent years should have double weight as compared to older data. 
#' @returns A data.table with all columns given in cols, and additional columns showing the model formula, bpue estimates, lower and upper confidence intervals, and a logical indicating model heterogeneity in the base model (I^2). See details.
#' @seealso [calc_total()]
#' @export

calc_total <- function(bpue, cols = c("ecoregion", "metierl4", "species"), obs, allx, verbose = TRUE, years=years, include.weights=FALSE) {

    # parallelization support
    if (nrow(bpue) > 1) {
        ret <- foreach(i = 1:nrow(bpue), 
                       .export = "calc_total", # <- not 100% sure this line is needed.
                       .final = rbindlist,
                       .packages = c("data.table", "glmmTMB", "emmeans", "ggeffects")) %dopar% {
                           calc_total(bpue = bpue[i], cols = cols, obs = obs, allx = allx, verbose = FALSE, years=years, include.weights=FALSE)
                       }
        return(ret)
    }
    
    obs <- obs[bpue, on = cols]
    obs <- obs[daysatsea > 0]
    #possible_re <- c("country", "areacode", "year", "metierl5", "vessellength_group",
    #                  "samplingprotocol", "monitoringmethod")
    #obs[,(possible_re) := lapply(.SD, as.factor), .SDcols = possible_re]
    obs[,(cols) := lapply(.SD, as.factor), .SDcols = cols]
    obs[, logDAS := log(daysatsea)]
    
    # add observation weights to data(if monitoring within the last five years 
    # use full weight, if older data weight observations half as strongly in the likelihood)
    obs[, weights := ifelse(year >= (max(years)-4), 1, 0.5)] 

    ret <- bpue[, c(cols, "model"), with = FALSE]
    ret[, c("tot_mean", "tot_lwr", "tot_upr", "message", "fishing_effort") :=
            list(NA_real_, NA_real_, NA_real_, "OK", NA_real_)]

    if (is.na(bpue$model) | bpue$model == "only one") {
        return(ret)
    }
    
    
    ##change post WKBBEAM to account for the zero that pass the threshold resulting in a BPUE = 0 estimate. 
    ## note here 
    if (bpue$model == "none") {
      ret$tot_mean<-bpue$bpue   # 0 if the zero passed else NA
      return(ret)
    }
    

    form <- as.formula(bpue$model)
    re <- lme4::findbars(form) # random effects part of model formulation (if any)
    re <- sapply(re, function(x) as.character(x[[3]]))
    re.n <- length(re)

    # don't try to generate total estimates for models where either 
    # samplingProtocol or monitoringMethod is among the random effects
    if (re.n > 0 & any(c("samplingprotocol", "monitoringmethod") %in% re)) {
        ret$message <- "samplingprotocol or monitoringmethod not ok"
        return(ret)
    }
    
    # don't try to generate total estimates for models with unexplained heterogeneity (base_model_heterogeneity = TRUE but best model is n_ind~1)
    if (re.n == 0 && isTRUE(bpue$base_model_heterogeneity)) {
      ret$message <- "unable to explain heterogeneity with candidate random effects"
      return(ret)
    }
    
    
    # sum up total fishing effort per combination of all levels of the random effects
    tot <- allx[bpue, on = cols[cols != "species"], .(das_fishing = sum(daysatseaf), final_year = max(year)), by = re]
    #tot[, logDAS := log(das)] # OUTI: Days at sea is on a log scale
    tot <- tot[complete.cases(tot)]
    
    
    if (nrow(tot) == 0) {
        return(ret)
    }
    
    # don't generate estimates unless the levels of all random effects match 
    # between monitored and total effort data (NOTE!! There are cases where monitoring effort strata do not exist in the total effort stratas)
    if (any(sapply(re, function(re) !all(tot[[re]] %in% unique(obs[[re]]))))) {
        ret$message <- "levels for at least one random effect not ok"
        return(ret)
    }
    
    # Refit best model
    if(isTRUE(include.weights)){
    best <- glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = obs, weights = weights)
    } else {
    best <- glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = obs)  
    }
    
    if (re.n == 0) {
      
        tot_obs <- obs
        tot_obs <- tot_obs[year==tot$final_year,
                           .(das_monitoring = sum(daysatsea), 
                           n_ind = sum(n_ind))]
        
        # join total fishing effort and monitored effort per random factor
        tot[,names(tot_obs) := tot_obs]
        tot <- tot[,das_no_monitoring := das_fishing - das_monitoring]
        tot[, logDAS := log(das_no_monitoring)] # OUTI: Days at sea is on a log scale
        
        
        pred <- as.data.frame(emmeans(best, ~1, offset = tot$logDAS, type = "response"))
        ret[, c("tot_mean", "tot_lwr", "tot_upr") := as.list(pred[, c("response", "asymp.LCL", "asymp.UCL")] + tot$n_ind)]  # prediction for unmonitored fishing effort + observed bycatch in monitoring
        ret$fishing_effort <- sum(tot$das_fishing) # retain total fishing effort or total unmonitored fishing effort?
        
    } else {
              
        # sum monitoring effort, and number of bycaught animals, per random factor for year = THISYEAR-1
        tot_obs <- obs[,(re) := lapply(.SD, as.factor), .SDcols = re]
        tot_obs <- tot_obs[year==unique(tot$final_year),
                           .(das_monitoring = sum(daysatsea), 
                           n_ind = sum(n_ind)), 
                           by = re]
        
        tot[, (re) := lapply(.SD, as.factor), .SDcols = re] # Do we treat YEAR correctly here? or does year need to be exempted?
        # tot only has THISYEAR-1 data so if year %in% re then we only predict for THISYEAR-1
        
        # join total fishing effort and monitored effort 
        tot <- tot_obs[tot, on=re]
        tot[is.na(das_monitoring), das_monitoring := 0] # replace na with 0
        tot[is.na(n_ind), n_ind := 0] # replace na with 0
        tot <- tot[,das_no_monitoring := das_fishing - das_monitoring] #
        tot[, logDAS := log(das_no_monitoring)] # OUTI: Days at sea is on a log scale
      

        pred <- lapply(1:nrow(tot), function(i) {
          type_arg <- ifelse(packageVersion("ggeffects") >= "2.0.0", "random", "re")
           
          p <- ggpredict(model = best,
                           terms = tot[i, ..re],
                           condition = c(logDAS = tot$logDAS[i]),
                           type = type_arg,
                           interval = "confidence",
                           verbose = verbose)
         
            p <- as.data.frame(p) 
            data.table(mean = p$predicted,
                       lwr = ifelse(!is.null(p[["conf.low"]]), p[["conf.low"]], NA_real_),
                       upr = ifelse(!is.null(p[["conf.high"]]), p[["conf.high"]], NA_real_))
        }) |> rbindlist()

        ret[, c("tot_mean", "tot_lwr", "tot_upr") := as.list(colSums(pred) + sum(tot$n_ind))] # prediction for unmonitored fishing effort + observed bycatch in monitoring
        ret$fishing_effort <- sum(tot$das_fishing) # retain total fishing effort or total unmonitored fishing effort?
        
    }
    
    return(ret)
}

