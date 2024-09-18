
calc_total <- function(bpue,effort,monitor, verbose = TRUE) {


require(ggeffects)
require(data.table)
require(glmmTMB)
require(metafor)
require(emmeans)
require(glmmTMB)


    mon <- monitor[bpue, on = c("ecoregion", "metierL4", "species")] 
    mon <- mon[daysAtSea > 0]
	cols<-c("country", "areaCode","year", "metierL5", "vesselLength_group","samplingProtocol","monitoringMethod")
	mon[,(cols) := lapply(.SD, as.factor), .SDcols = cols] 
	mon$lDaS<-log(mon$daysAtSea)
	
	
    pred <- list(mean = NA_real_, lwr = NA_real_, upr = NA_real_, message = "OK", effort = NA_real_)

    if (is.na(bpue$model)) {
        return(pred)
    }
    if (bpue$model=="only one") {
        return(pred)
    }
    
    form <- as.formula(bpue$model)
    re <- lme4::findbars(form) # random effects part of model formulation (if any)
    re <- sapply(re, function(x) as.character(x[[3]]))
    re.n <- length(re)
    
    # We can't generate total bycatch estimates for models using these specific
    # random effects, so let's just stop here.
    if (re.n > 0 & any(c("samplingProtocol", "monitoringMethod") %in% re)) {
        pred$message <- "samplingProtocol or monitoringMethod not ok"
        return(as.data.table(pred))
    }
    
    tot <- effort[bpue, on = c("ecoregion", "metierL4"), .(das = (sum(daysAtSeaF))), by = re]
	tot$lDaS<-log(tot$das)
	
    # check if all the levels of the random effect match between monitored and total effort data
    # and if they're not stop here
    if (any(sapply(re, function(re) !all(tot[[re]] %in% unique(mon[[re]]))))) {
        pred$message <- "random levels for at least one random effect not ok"
        return(as.data.table(pred))
    }
    
	# sum up total fishing effort per combination of all levels of the random effects
    
    # refit the best model. optimization: we might want to save the best fitted 
    # model in calc_bpue, and readRDs it in here, rather then refit it.
    best <- glmmTMB(formula = form, offset = lDaS, family = nbinom2, data = mon)
    
	
	
    if (re.n == 0) {
		tot <- effort[bpue,on = c("ecoregion", "metierL4"), .(das = (sum(daysAtSeaF)))]
        tot$lDaS<-log(tot$das)
        pred[1:3] <- as.numeric(as.data.frame(emmeans(best, ~1, offset = tot$lDaS, type = "response"))[, c("response", "asymp.LCL", "asymp.UCL")])
        pred[5]<-sum(tot$das)
		##change here effort2 is tot note tot das is already logged there
		
    } else {
    tot<-tot[complete.cases(tot)]
	    
        pred.re <- lapply(1:nrow(tot), function(i) {
            p <- ggpredict(model=best,
							terms= tot[i, ..re],
							condition=c(lDaS = tot$lDaS[i]),
							type="re",
							interval="confidence")
	
			if (ncol(p)>2) {
            as.data.frame(p)[, c("predicted", "conf.low", "conf.high")] }
        })
		
        pred[1:3] <- as.numeric(colSums(rbindlist(pred.re))) 
        pred[5]<-sum(tot$das)
    }
        
    
    return(as.data.table(pred))
}

