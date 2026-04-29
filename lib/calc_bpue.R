
#' Calculates BPUE (bycatch per unit effort) for ICES WGBYC data
#' 
#' @description
#' `calc_bpue` takes monitored fishing effort data obtained in the ICES WGBYC
#' annual data call, and uses those data to generate a number of models of varying
#' complexity. It then compares those models using the Akaike Information Criterion (AIC), combined with a stability 
#' filter and a parsimony rule
#' @details
#' If there are multiple rows in `needle`, each row is processed separately, and
#' rbind'ed into one data.table with nrow equal to the number of rows in `needle`.
#' If a parallel backend is available, the code will execute in parallel, allowing
#' faster computations.

#' The model selection uses **AIC** with a **Parsimony Rule**:
#' 1. **Model Selection**: Valid combinations of random effects are fitted. Models with 
#'    "Singular Fit" (near-zero variance or perfect correlation) are 
#'    discarded to ensure numerical stability and prevent over-parameterization.
#' 2. **Tie-breaker**: Among stable candidates, if the difference in AIC is < 2, 
#'    the simpler model (fewer random effects) is chosen, but flagged under "alternative_models_flag".

#' @param needle data.table with values for `cols`, used to subset observations from `dat`. Any additional columns in `needle` are disregarded.
#' @param cols columns, whose unique combinations, are used to split the data given in `dat`
#' @param min_re_obs Integer specifying the minimum number of levels needed to include a term as a random effect
#' @param dat data.table with monitored fishing effort data (e.g. since 2017)
#' @param years vector with integers indicating years of assessment. 
#' @param include.weights boolean indicator indicating whether observations from the five most recent years should have double weight as compared to older data. 
#' @returns A data.table with one row for each row in `needle` and all columns given in `cols`, plus additional columns showing the final model formula, BPUE estimates, lower and upper confidence intervals, and a logical indicating model heterogeneity in the base model (I^2). See details.
#' @seealso [calc_total()]
#' @export
#' 
BEAM_progress <- function(n) {
    BEAM_pb$tick(tokens = list(step = n))
}

calc_bpue <- function(needle, cols = colnames(needle), min_re_obs = 2, dat, years, include.weights = FALSE) {

    t_start <- Sys.time()
    
    # parallelization support
    if (nrow(needle) > 1) {
        
        if (exists("BEAM_pb")) {
            BEAM_pb$terminate()
            rm("BEAM_pb", envir = .GlobalEnv)
			 }
        
        BEAM_pb <<- progress_bar$new(
            format = "calc_bpue :percent :current/:total [:bar] :elapsed | eta: :eta",
            total = nrow(needle),
            width = 60)
        
        opts <- list(progress = BEAM_progress)
        
        ret <- foreach(i = 1:nrow(needle),
                       .export = "calc_bpue", # <- not 100% sure this line is needed.
                       .final = rbindlist,
                       .packages = c("data.table", "glmmTMB", "metafor", "emmeans"),
                       .options.snow = opts) %dopar% {
            calc_bpue(needle = needle[i], cols = cols, min_re_obs = min_re_obs, dat = dat, include.weights)
                       }
        
        BEAM_pb$terminate()
        rm("BEAM_pb", envir = .GlobalEnv)
        return(ret)
    }

    dat <- dat[needle[, ..cols], on = cols, nomatch = 0]
    dat <- dat[daysatsea > 0]
    dat[, (cols) := lapply(.SD, as.factor), .SDcols = cols]
    dat[, logDAS := log(daysatsea)]
	dat[, n_ind := as.integer(n_ind)] 
    
    # add observation weights to data(if monitoring within the last five years 
    # use full weight, if older data weight observations half as strongly in the likelihood)
    dat[, weights := ifelse(year >= (max(years)-4),1,0.5)] 
    
    ret <- needle[, ..cols]
    ret[, c("model", "bpue", "lwr", "upr", "replicates", "base_model_heterogeneity", "alternative_models_flag") :=
            list("none", NA_real_, NA_real_, NA_real_, nrow(dat), NA, FALSE)]

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
    
    # fit model either including or excluding likelihood weights
    if (isTRUE(include.weights)){
      base_model <- tryCatch(glmmTMB(n_ind ~ 1, offset = logDAS, family = nbinom2, data = dat, weights = weights), error = function(e) e$message)
    } else {
      base_model <- tryCatch(glmmTMB(n_ind ~ 1, offset = logDAS, family = nbinom2, data = dat), error = function(e) e$message)
    }
    
    heterogeneity <- tryCatch((rma.glmm(xi = n_ind, ti = daysatsea, measure = "IRLN", data = dat)$QEp.Wld<0.05), error = function(e) e$message)
    
    
    if (!inherits(heterogeneity, "character")) {
        ret$base_model_heterogeneity <- heterogeneity
    }
    
    best <- base_model # if we have less than 5 rows of monitoring data, fit a simple model and we're done.

    # if we have more than 5 rows of monitoring data, fit all possible
    # combinations of the base model and one or more of the terms in the
    # vector below, added to the model as random effects.
    if (nrow(dat) >= 5) {

        re_terms <- c("country", "areacode", "year", "quarter", "metierl5", 
                "vessellength_group", "samplingprotocol", "monitoringmethod")

        # but only consider r.e. terms where the number of unique values 
        # (i.e. levels) is greater than min_re_obs. Note that this effectively
        # prevents any variables specified in the cols parameter from being included
        # as random effects in any models, since they will always have length=1.
        re_terms <- re_terms[sapply(re_terms, function(x) length(unique(dat[[x]]))) >= min_re_obs]
        re.i <- do.call(CJ, replicate(length(re_terms), c(TRUE, FALSE), simplify = FALSE)) # all possible combinations of random effect levels

		re_candidates <- apply(re.i, 1, function(i) {
        if (all(i == FALSE)) return("1")
        terms <- sprintf("(1|%s)", re_terms[unlist(i)])
        return(paste("1 +", paste(terms, collapse = " + ")))
        })							

	## MP: 09/02/2026
	## Integration of two versions of calc_bpue:
	# New model selection function written by KMB, including singularity checks and parsimony rule, 
	# New option of using weights in models, developped by TS

	  candidates <- lapply(re_candidates, function(f_str) {
      curr_formula <- as.formula(sprintf("n_ind ~ %s", f_str))
      fit <- suppressMessages(tryCatch({
		 if(isTRUE(include.weights)){
             m <- glmmTMB(formula = curr_formula, offset = logDAS, family = nbinom2, data = dat, weights = weights)
            }else{
             m <- glmmTMB(formula = curr_formula, offset = logDAS, family = nbinom2, data = dat)
            }
       
        # Check for Singularity in glmmTMB
        vc <- VarCorr(m)$cond
        is_singular <- FALSE
        for (comp in vc) {
          if (any(diag(comp) < 1e-6)) is_singular <- TRUE # Variance near zero
          if (nrow(comp) > 1) {
            cor_mat <- attr(comp, "correlation")
            if (any(abs(cor_mat[lower.tri(cor_mat)]) > 0.99)) is_singular <- TRUE # Near perfect correlation
          }
        }
        
        if (!is_singular) return(m) else return("Singular") 
      }, error = function(e) e$message))
      return(fit)
    })
    
    # Identify Valid glmmTMB objects
    valid_indices <- which(sapply(candidates, function(x) inherits(x, "glmmTMB")))
    
    if (length(valid_indices) > 0) {
      # Use AIC instead of BIC for more balanced complexity selection
      valid_aics <- sapply(candidates[valid_indices], AIC)
      valid_complexity <- lengths(regmatches(re_candidates[valid_indices], gregexpr("\\+", re_candidates[valid_indices])))
      
      sel_table <- data.table(
        Index = valid_indices,
        AIC = valid_aics,
        Complexity = valid_complexity
      )
      
      sel_table <- sel_table[is.finite(AIC)]
      
      if (nrow(sel_table) > 0) {
        min_aic <- min(sel_table$AIC)
        # Zone of indifference (delta AIC < 2)
        top_tier <- sel_table[AIC <= (min_aic + 2)]
        
        # Parsimony Rule: Choose simplest model in the top tier
        setorder(top_tier, Complexity, AIC)
        
        best_idx <- top_tier$Index[1]
        best <- candidates[[best_idx]]
        
        if (nrow(top_tier) > 1) {
          ret$alternative_models_flag <- TRUE
        }
      }
    }
  }
  
  # Final check for model validity before emmeans
  if (inherits(best, "glmmTMB")) {
    bpue.r <- tryCatch(as.data.frame(emmeans(best, ~1, type="response", offset = log(1))),
                       error = function(e) NULL)
    if (!is.null(bpue.r)) {
      ret[, c("bpue", "lwr", "upr") := as.list(bpue.r[, c("response", "asymp.LCL", "asymp.UCL")])]
    }
  }
  
  # Format the model formula as a string
  if (inherits(best, "glmmTMB")) {
    ret[, model := paste(format(formula(best)), collapse = "")]
  } else {
    ret[, model := as.character(best)]
  }

	# In cases where no heterogeneity was detected, over-write best model by base model

  if (isFALSE(ret$base_model_heterogeneity)) {
    ret[, model := paste(format(formula(base_model)), collapse = "")]
  } else {}
									

## MP: below is Legacy code - included weight option but old model selection approach - new version should integrate both but keeping the code below in case bugs occur.
#########							
        # fit all candidate models 
        #candidates <- apply(re.i, 1, function(i) {
            
           # if (all(i == FALSE)) {
                #return(base_model)
            #}
            
            #re <- sprintf("(1|%s)", re[unlist(i)])
            #form <- sprintf("n_ind ~ 1 + %s", paste(re, collapse = " + "))
            #form <- as.formula(form)
            # wrapped this in suppressMessages just to avoid cluttering of the console when running glmmTMB, errors and warnings are still caught by tryCatch.
            
            # fit model either including or excluding model weights
            #if(isTRUE(include.weights)){
            #  suppressMessages(tryCatch(glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = dat, weights = weights), error = function(e) e$message, warning = function(w) w$message))  
            #}else{
            #  suppressMessages(tryCatch(glmmTMB(formula = form, offset = logDAS, family = nbinom2, data = dat), error = function(e) e$message, warning = function(w) w$message))  
            #}
  
       # })
        
        # pick the best one
        #candidates.converged <- candidates[which(sapply(candidates, class) != "character")]
        #scores <- sapply(candidates.converged, AIC)
        #if (length(scores) > 1 & sum(!is.na(scores)) > 0) {
        #    best <- candidates.converged[[which.min(scores)]]
        #}
    #}
    
    #if (class(best) != "character") {
        #bpue.r <- as.data.frame(emmeans(best, ~1, type="response", offset = log(1))) # BPUE for one DaS we can push to total effort prediction instead next
        #ret[, c("bpue", "lwr", "upr") := as.list(bpue.r[, c("response", "asymp.LCL", "asymp.UCL")])]
    #}

    #ret[, model := paste(format(formula(best)), collapse = "")]
    ###################
							   
    t_end <- Sys.time()
    t_elapsed <- difftime(t_end,t_start,units="mins")
    cat(sprintf("calc_bpue on %d rows completed in %.01f mins\n", nrow(needle), t_elapsed))
    return(ret)
    
}




