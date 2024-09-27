
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
#' @returns A data.table with all columns given in cols, and additional columns showing the model formula, bpue estimates, lower and upper confidence intervals, and a logical indicating model heterogeneity in the base model (I^2). See details.
#' @seealso [calc_total()]
#' @export
#' 
#' 
#' 
# There is an error in ggpredict or one of its dependencies that causes an issue with factor
# variables that have a level "de" (it may somehow conflict with src/library/utils/R/de.R)
# a temp fix is just to substitute de with something else in all 

calc_total <- function(bpue, cols = c("ecoregion", "metierl4", "species"), obs, all, verbose = TRUE) {

    # parallelization support
    if (nrow(bpue) > 1) {
        ret <- foreach(i = 1:nrow(bpue), 
                       .export = "calc_total", # <- not 100% sure this line is needed.
                       .final = rbindlist,
                       .packages = c("data.table", "glmmTMB", "emmeans", "ggeffects")) %dopar% {
                           calc_total(bpue = bpue[i], cols = cols, obs = obs, all = all, verbose = FALSE)
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

    ret <- bpue[, c(cols, "model"), with = FALSE]
    ret[, c("mean", "lwr", "upr", "message", "effort") :=
            list(NA_real_, NA_real_, NA_real_, "OK", NA_real_)]

    if (is.na(bpue$model) | bpue$model == "none" | bpue$model == "only one") {
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
    
    # sum up total fishing effort per combination of all levels of the random effects
    tot <- all[bpue, on = cols[cols != "species"], .(das = sum(daysatseaf)), by = re]
    tot[, logDAS := log(das)]
    tot <- tot[complete.cases(tot)]
    
    if (nrow(tot) == 0) {
        return(ret)
    }
    
    # don't generate estimates unless the levels of all random effects match 
    # between monitored and total effort data
    if (any(sapply(re, function(re) !all(tot[[re]] %in% unique(obs[[re]]))))) {
        ret$message <- "levels for at least one random effect not ok"
        return(ret)
    }
    
    # temp fix, part 1:
    if ("country" %in% re) {
        obs[country == "de", country := "de_temp_fix"]
        tot[country == "de", country := "de_temp_fix"]
    }

    best <- glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = obs)

    if (re.n == 0) {

        pred <- as.data.frame(emmeans(best, ~1, offset = tot$logDAS, type = "response"))
        ret[, c("mean", "lwr", "upr") := as.list(pred[, c("response", "asymp.LCL", "asymp.UCL")])]
        ret$effort <- sum(tot$das)
        
    } else {
        
        tot[, (re) := lapply(.SD, as.factor), .SDcols = re] # Do we treat YEAR correctly here? or does year need to be exempted?
        
        pred <- lapply(1:nrow(tot), function(i) {
            p <- ggpredict(model = best,
                           terms = tot[i, ..re],
                           condition = c(logDAS = tot$logDAS[i]),
                           type = "re",
                           interval = "confidence",
                           verbose = verbose)
            
            p <- as.data.frame(p) 
            data.table(mean = p$predicted,
                       lwr = ifelse(!is.null(p[["conf.low"]]), p[["conf.low"]], NA_real_),
                       upr = ifelse(!is.null(p[["conf.high"]]), p[["conf.high"]], NA_real_))
        }) |> rbindlist()

        ret[, c("mean", "lwr", "upr") := as.list(colSums(pred))]
        ret$effort <- sum(tot$das)
        
    }
    
    # temp fix, part 2 (convert back)
    if ("country" %in% re) {
        obs[country == "de_temp_fix", country := "de"]
        tot[country == "de_temp_fix", country := "de"]
    }
    return(ret)
}


