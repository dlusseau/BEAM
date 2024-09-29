
calc_bpue_area <- function(areaCode, ml4, Species, min_re_obs = 2, ices = filtered_monitor, 
                           model_dir = "area_models/") {
  require(ggeffects)
  require(data.table)
  require(glmmTMB)
  require(metafor)
  require(emmeans)
  
  ices$lDaS <- log(ices$daysAtSea)
  
  needle <- data.table(areaCode = areaCode, metierL4 = ml4, species = Species)
  dat <- ices[needle, on = .(areaCode, metierL4, species), nomatch = NULL]
  
  cols <- c("country", "areaCode", "year", "metierL5", "vesselLength_group", "samplingProtocol", "monitoringMethod")
  dat[ , (cols) := lapply(.SD, as.factor), .SDcols = cols]
  
  print(paste("Sum of n_ind:", sum(dat$n_ind)))
  
  # Return early if no data or sum of individuals is zero
  if (nrow(dat) == 0 || sum(dat$n_ind) == 0){
    print("No individuals. Skipping model fitting.")
    return(list(bpue = NA, lwr = NA, upr = NA, model = "No data or zero bycatch", replicates = nrow(dat),
                base_model_heterogeneity = NA, bpue.cond.name = NA, bpue.cond = NA, bpue.cond.lwr = NA, bpue.cond.upr = NA))
  }
  
  # If there's only one row of data
  if (nrow(dat) == 1) {
    bpue <- dat$n_ind / dat$daysAtSea
    lwr <- bpue - 1.96 * sqrt(dat$n_ind / dat$daysAtSea^2)
    upr <- bpue + 1.96 * sqrt(dat$n_ind / dat$daysAtSea^2)
    return(list(bpue = bpue, lwr = lwr, upr = upr, model = "only one", replicates = nrow(dat),
                base_model_heterogeneity = NA, bpue.cond.name = NA, bpue.cond = NA,
                bpue.cond.lwr = NA, bpue.cond.upr = NA))
  }
  
  # Ensure model directory exists
  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }
  
  # Construct file path for model
  model_file <- file.path(model_dir, paste0("model_", areaCode, "_", ml4, "_", Species, ".rds"))
  print(paste("Model file path:", model_file))
  
  # Check if model already exists
  if (file.exists(model_file)) {
    print("Model file exists. Loading model.")
    best_model <- readRDS(model_file)
  } else {
    print("Model file does not exist. Fitting models.")
    
    # Fit base model
    base_model <- tryCatch(
      glmmTMB(n_ind ~ 1, offset = lDaS, family = nbinom2, data = dat),
      error = function(e) {
        print(paste("Error in base model fitting:", e$message))
        NULL
      }
    )
    
    if (nrow(dat) < 5) {
      best_model <- base_model
      # Fit heterogeneity base model (rma.glmm)
      heterogeneity_base <- tryCatch(
        (rma.glmm(xi = n_ind, ti = daysAtSea, measure = "IRLN", data = dat)$QEp.Wld < 0.05),
        error = function(e) e$message
      )
    } else {
      re <- c("country", "year", "metierL5", "vesselLength_group", "samplingProtocol", "monitoringMethod")
      re <- re[sapply(re, function(x) length(unique(dat[[x]])) >= min_re_obs)]
      
      if (length(re) == 0) {
        return(list(bpue = NA, lwr = NA, upr = NA, model = "No valid random effects", replicates = nrow(dat)))
      }
      
      # Create combinations of random effects to try
      re_combinations <- do.call(CJ, replicate(length(re), c(TRUE, FALSE), simplify = FALSE))
      
      # Fit models for each combination of random effects
      candidates <- lapply(1:nrow(re_combinations), function(i) {
        if (all(re_combinations[i, ] == FALSE)) return(base_model)
        re_formula <- sprintf("(1|%s)", re[unlist(re_combinations[i, ])])
        formula <- as.formula(sprintf("n_ind ~ 1 + %s", paste(re_formula, collapse = " + ")))
        tryCatch(
          glmmTMB(formula = formula, offset = lDaS, family = nbinom2, data = dat),
          error = function(e) {
            print(paste("Error fitting model with formula", formula, ":", e$message))
            NULL
          }
        )
      })
      
      # Filter out failed models and select the best one by AIC
      candidates_converged <- Filter(Negate(is.null), candidates)
      best_model <- if (length(candidates_converged) > 0) {
        candidates_converged[[which.min(sapply(candidates_converged, AIC))]]
      } else {
        base_model
      }
      
      # Fit heterogeneity base model using rma.glmm
      heterogeneity_base <- tryCatch(
        rma.glmm(xi = n_ind, ti = daysAtSea, measure = "IRLN", data = dat)$QEp.Wld < 0.05,
        error = function(e) e$message
      )
      
      # Save the best model as an RDS file
      cat("Saving best model to:", model_file, "\n")
      saveRDS(best_model, model_file)
    }
  }
  
  # Generate predictions and confidence intervals
  if (!inherits(best_model, "character")) {
    bpue_r <- as.data.frame(emmeans(best_model, ~1, type = "response", offset = log(1)))
  } else {
    bpue_r <- data.frame(response = NA, asymp.LCL = NA, asymp.UCL = NA)
  }
  
  # Check if the best model differs from the base model
  if (formula(base_model) != formula(best_model)) {
    return(list(
      bpue = bpue_r$response,
      lwr = c(bpue_r$lower.CL, bpue_r$asymp.LCL),
      upr = c(bpue_r$upper.CL, bpue_r$asymp.UCL),
      model = deparse(formula(best_model)),
      replicates = nrow(dat),
      base_model_heterogeneity = heterogeneity_base,
      bpue.cond.name = NA,
      bpue.cond = NA,
      bpue.cond.lwr = NA,
      bpue.cond.upr = NA
    ))
  } else {
    return(list(
      bpue = bpue_r$response,
      lwr = c(bpue_r$lower.CL, bpue_r$asymp.LCL),
      upr = c(bpue_r$upper.CL, bpue_r$asymp.UCL), # Use either value depending on the model
      model = deparse(formula(best_model)),
      replicates = nrow(dat),
      base_model_heterogeneity = heterogeneity_base,
      bpue.cond.name = NA,
      bpue.cond = NA,
      bpue.cond.lwr = NA,
      bpue.cond.upr = NA
    ))
  }
}


