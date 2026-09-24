## ---------------------------
##
## Script name: RMSE_bpues_V5.R
##
## Purpose of script: Calculate root mean square error (RMSE) between BPUE estimates
##                    and leave-one-out jacknife estimates, to help assess the
##                    uncertainty of any given BPUE estimate
##
## Authors: ToR C members
##
## ---------------------------

# 1) Mean BPUE + CV of estimates (emmeans)
# 2) Using LOO (=jacknife) so using ggpredict/emmeans we obtain a BPUE for the full model
#   and a BPUE for the jacknifed. Then, (BPUE_full-BPUE_jacknifed)^2 [we'll have 
#   as many values as rows of data for that unique combination of ecoregionxmetierxsps].
#   Then, sum all those rests and divide them by n-1 (cases) and make the square root

# Summary: fit full model + fit jacknifed models + ggpredict both + do (full_pred-jacknifed_pred)^2
# 

# FILTER WHEN BPUE !IS.NA

# ONLY FOR THOSE THAT WE HAVE A BPUE ESTIMATE

# pa filter, del bpue dataframe, cada fila ya tiene la combinacion unica de cada ecoxmetierxsps
# despues recuperar del full dataframe 

# CHANGELOG:
# 2025-10-15 AM Refactored code to take full advantage of data.table performance
# benefits and support parallelization.
# TODO: move calc_rmse into separate file calc_rmse.R? What do we do with the
# rest of the code in this script then?
# TODO: look into unresponsive console issue that sometimes occurs after 
# parallel execution

library(glmmTMB)
library(ggeffects)
library(emmeans)
library(data.table)

# Load data ----

#Moved data reads outside the function so it can be used both for aggregated and disaggregated data.
#full_data <- fread("data/obs3.csv")
#bpues_estimates <- fread("data/bpue1.csv")

# Preparing datasets
bpues_estimates$RMSE <- NA
# Remove NAs in BPUE estimates & 0 in obs3 DaS
bpues_estimates <- bpues_estimates[!is.na(bpue), ]
bpues_estimates <- bpues_estimates[model != "none", ]
bpues_estimates <- bpues_estimates[model != "only one", ]

# Reduce the size removing NAs in totalbycatch
# AM 15/10: changed to use data.table syntax. We should probably avoid spaces in var names...
# bpues_estimates <- bpues_estimates[!is.na(bpues_estimates$total.bycatch.2024), ]
#bpues_estimates <- bpues_estimates[!is.na(`total bycatch 2024`), ] Removing this for MAMBY

#complete_bpues$RMSE <- 0
full_data <- full_data[daysatsea > 0, ]
full_data <- full_data[!is.na(daysatsea), ]
full_data <- full_data[taxon_bycatch_monitor_ok == TRUE, ]
# Add log of monitoring (DaS)

full_data[, logDAS := log(daysatsea)]
#full_data$logDAS <- log(full_data$daysatsea)

## this is a temporary fix introduced in 2025 to deal with a new request for which data will come
#for this particular case---no anguilla in full data
#bpues_estimates<-bpues_estimates[Species!="anguilla anguilla",]
#bpues_estimates<-bpues_estimates[Species!="salmo salar",]

calc_rmse <- function(full_data, bpue) {
  
    # helper function (to avoid repeated code)
    fit_and_predict <- function(form, data, levels = NULL) {
        fit <- suppressMessages(glmmTMB(formula = form, offset = logDAS, data = data, family = nbinom2))
        re <- lme4::findbars(form) # random effects part of model formulation (if any)
        re <- sapply(re, function(x) as.character(x[[3]]))
        re.n <- length(re) # number of random effects

        if (re.n == 0) {
            pred <- as.data.table(emmeans(fit, ~1, type = "response", offset = log(1)))
            pred$predicted <- pred$response
        } else {
           # type_arg <- ifelse(packageVersion("ggeffects") >= "2.0.0", "random", "re")
          type_arg<-"random"
            pred <- ggpredict(model = fit, type =  type_arg, terms = re, verbose = FALSE)
            pred <- as.data.table(pred)
            
            if (!"facet" %in% colnames(pred)) {
                pred[, facet := "dummy"]
            }
            # if (!is.null(levels) & inherits(levels, "list")) {
            #     for (i in 1:length(levels)) {
            #         if (length(levels[[i]])) {
            #             nm <- names(levels)[i]
            #             pred[, (nm) := factor(get(nm), levels = levels[[i]])]
            #         }
            #     }
            # }
            pred <- pred[order(x, group, facet)]
        }

        list(fit = fit, re = re, re.n = re.n, pred = pred)
    }
    
    # parallelization support
    if (nrow(bpue) > 1) {
        cat("calc_rmse: setting up workers")
        
        if (exists("BEAM_pb")) {
            BEAM_pb$terminate()
            rm("BEAM_pb", envir = .GlobalEnv)
        }
        BEAM_pb <<- progress_bar$new(
            format = "calc_rmse :percent :current/:total [:bar] :elapsed | eta: :eta",
            total = nrow(bpue),
            width = 60)
            
        opts <- list(progress = BEAM_progress)
            
        ret <- foreach(i = 1:nrow(bpue),
                       .export = "calc_rmse",
                       .final = unlist,
                       .packages = c("data.table", "glmmTMB", "ggeffects", "emmeans"),
                       .options.snow = opts) %dopar% {
                         options(expressions = 500000)
                           calc_rmse(full_data, bpue = bpue[i])
                       }
            
        BEAM_pb$terminate()
        rm("BEAM_pb", envir = .GlobalEnv)
        return(ret)
    }
    
    needle <- bpue[, .(ecoregion = ecoregion, metierl4 = metierl4, species = species)] # Changed that for MAMBY
    data_subset <- full_data[needle, on = c("ecoregion", "metierl4", "species")]

    # Extract the formula of the model
    formula <- bpue[["model"]]
    #formula <- sprintf("%s + offset(logDAS)", formula)
    form <- as.formula(formula)
    
    # Fit full model
    full_model <- fit_and_predict(form, data_subset)

    #factor_levels <- list(
    #    x = levels(full_model$pred$x),
    #    group = levels(full_model$pred$group),
    #    facet = levels(full_model$pred$facet)
    #)
    
    diffs <- sapply(1:nrow(data_subset), function(j) {
        #jack <- fit_and_predict(form, data = data_subset[-j,], levels = factor_levels)
        jack <- fit_and_predict(form, data = data_subset[-j,])
        log_diff <- log1p(full_model$pred$predicted) - log1p(jack$pred$predicted)
        n <- ifelse(jack$re.n == 0, 1, full_model$re.n)
        sum(log_diff^2, na.rm = TRUE) / n
    })
    
    sqrt(sum(diffs, na.rm=TRUE)/(nrow(data_subset)-1))
}


# All species over the new obs3 and bpue_print_table
rmse <- calc_rmse(full_data = full_data, bpue = bpues_estimates)

# sometimes the above call seems to make the R console unresponsive, even after
# code execution has completed. I'm not sure why this occurs, and it only happens
# sometimes. But we can fix it by calling sink() or closeAllConnection()
# See https://stackoverflow.com/questions/26495498/rstudio-does-not-display-any-output-in-console-after-entering-code

bpues_estimates[, RMSE := rmse]

# Transform back from log scale
bpues_estimates$factor <- exp(bpues_estimates$RMSE)


#Factor=exp(RMSE) # RMSE here is the variable you have in your function but computed with the log modification above
# The lower limit of the 68% interval is: 100(1-1/factor) # in percentage
bpues_estimates$lower_factor_prop <- 100*(1-1/bpues_estimates$factor)
# The upper limit of the 68% interval is: 100(factor−1) # in percentage
bpues_estimates$upper_factor_prop <- 100*(bpues_estimates$factor-1)

# RELIABILITY LIMIT ON BOTH CONFIDENCE INTERVALS OF FACTOR
# Take the 25% limit of the lower % CI of factor (=80% similarity)
# Take the 25% limit of the upper % CIs of factor (=80% similarity)
# Which are below 25%

# Reliability
bpues_estimates$reliability_rmse <- bpues_estimates$lower_factor_prop <= 25 & bpues_estimates$upper_factor_prop <= 25

# # Plot rmse_log vs delta only for those that we are keeping
# reliable <- ggplot() +
  # geom_point(data = new.data[new.data$reliability == T, ], 
             # aes(RMSE_log, delta), col = "darkgreen") +
  # theme_bw()

# all_relvalues <- ggplot() +
  # geom_point(data = new.data, 
             # aes(RMSE_log, delta), col = "orange") +
  # theme_bw()

# # Before and after plot
# ggarrange(all_relvalues, reliable, nrow = 2)

# hist(new.data$delta[new.data$reliability == T], breaks = 200)

# SAVE SOME THINGS
fwrite(bpues_estimates, "results/bpue_reliability_rmse.csv", na="NA")

