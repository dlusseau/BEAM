
#' Calculates BPUE (bycatch per unit effort) for ICES WGBYC data
#' 
#' @description
#' `calc_bpue` takes monitored fishing effort data obtained in the ICES WGBYC
#' annual data call, and uses those data to generate a number of models of varying
#' complexity. It then compares those models using AIC, and returns the best model.
#' 
#' @details
#' If there are multiple rows in `needle`, each row is processed separately, and
#' rbind'ed into one data.table with nrow equal to the number of rows in `needle`.
#' If a parallel backend is available, the code will execute in parallel, allowing
#' faster computations.
#' 
#' @param needle data.table with values for `cols`, used to subset observations from `dat`. Any additional columns in `needle` are disregarded.
#' @param cols columns, whose unique combinations, are used to split the data given in `dat`
#' @param min_re_obs Integer specifying the minimum number of levels needed to include a term as a random effect
#' @param dat data.table with monitored fishing effort data (e.g. for the last 5 years)
#' @returns A data.table with one row for each row in `needle` and all columns given in `cols`, plus additional columns showing the final model formula, BPUE estimates, lower and upper confidence intervals, and a logical indicating model heterogeneity in the base model (I^2). See details.
#' @seealso [calc_total()]
#' @export
calc_bpue <- function(needle, cols = colnames(needle), min_re_obs = 2, dat) {

    # parallelization support
    if (nrow(needle) > 1) {
        ret <- foreach(i = 1:nrow(needle),
                       .export = "calc_bpue", # <- not 100% sure this line is needed.
                       .final = rbindlist,
                       .packages = c("data.table", "glmmTMB", "metafor", "emmeans")) %dopar% {
            calc_bpue(needle = needle[i], cols = cols, min_re_obs = min_re_obs, dat = dat)
        }
        return(ret)
    }

    dat <- dat[needle[, ..cols], on = cols, nomatch = 0]
    dat <- dat[daysatsea > 0]
    dat[, (cols) := lapply(.SD, as.factor), .SDcols = cols]
    dat[, logDAS := log(daysatsea)]
    
    ret <- needle[, ..cols]
    ret[, c("model", "bpue", "lwr", "upr", "replicates", "base_model_heterogeneity") :=
            list("none", NA_real_, NA_real_, NA_real_, nrow(dat), NA)]

    if ((nrow(dat) == 0) | (sum(dat$n_ind, na.rm = TRUE) == 0) | any(dat$n_ind < 0)) {
        return(ret)
    }
    
    # cop out when we only have 1 row of data
    if (nrow(dat) == 1) {
        bpue <- dat$n_ind / dat$daysatsea
        lwr <- bpue -1.96 * sqrt(dat$n_ind / dat$daysatsea^2)
        upr <- bpue +1.96 * sqrt(dat$n_ind / dat$daysatsea^2)
        ret[, c("bpue", "lwr", "upr", "model") := list(bpue, lwr, upr, "only one")]
        return(ret)
    }
    
    base_model <- tryCatch(glmmTMB(n_ind ~ 1, offset = logDAS, family = nbinom2, data = dat), error = function(e) e$message)
    heterogeneity <- tryCatch((rma.glmm(xi = n_ind, ti = daysatsea, measure = "IRLN", data = dat)$QEp.Wld<0.05), error = function(e) e$message)
    
    if (class(heterogeneity) != "character") {
        ret$base_model_heterogeneity <- heterogeneity
    }
    
    best <- base_model # if we have less than 5 rows of monitoring data, fit a simple model and we're done.

    # if we have more than 5 rows of monitoring data, fit all possible
    # combinations of the base model and one or more of the terms in the
    # vector below, added to the model as random effects.
    if (nrow(dat) >= 5) {

        re <- c("country", "areacode", "year", "metierl5", "vessellength_group",
                "samplingprotocol", "monitoringmethod")

        # but only consider r.e. terms where the number of unique values 
        # (i.e. levels) is greater than min_re_obs
        re <- re[sapply(re, function(x) length(unique(dat[[x]]))) >= min_re_obs]
        re.i <- do.call(CJ, replicate(length(re), c(TRUE, FALSE), simplify = FALSE))
        
        # fit all candidate models 
        candidates <- apply(re.i, 1, function(i) {
            
            if (all(i == FALSE)) {
                return(base_model)
            }
            
            re <- sprintf("(1|%s)", re[unlist(i)])
            form <- sprintf("n_ind ~ 1 + %s", paste(re, collapse = " + "))
            form <- as.formula(form)
            # wrapped this in suppressMessages just to avoid cluttering of the console when running glmmTMB, errors and warnings are still caught by tryCatch.
            suppressMessages(tryCatch(glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = dat), error = function(e) e$message, warning = function(w) w$message))
        })
        
        # pick the best one
        candidates.converged <- candidates[which(sapply(candidates, class) != "character")]
        scores <- sapply(candidates.converged, AIC)
        if (length(scores) > 1 & sum(!is.na(scores)) > 0) {
            best <- candidates.converged[[which.min(scores)]]
        }
    }
    
    if (class(best) != "character") {
        bpue.r <- as.data.frame(emmeans(best, ~1, type="response", offset = log(1))) # BPUE for one DaS we can push to total effort prediction instead next
        ret[, c("bpue", "lwr", "upr") := as.list(bpue.r[, c("response", "asymp.LCL", "asymp.UCL")])]
    }

    ret[, model := paste(format(formula(best)), collapse = "")]
    return(ret)
    
}


