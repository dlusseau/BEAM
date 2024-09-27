#calc_area_tot
# Read data
# bpue <- fread("~/ICES WGBYC/WGBYC2024/total_bycatch_version1.csv")
# monitor <- fread("~/ICES WGBYC/WGBYC2024/monitor_effort_bycatch_ecoregion_area_species_v2THELIST_17Sept.csv")
# effort <- fread("~/ICES WGBYC/WGBYC2024/D1_fishingeffort_2017_2023.csv")
# effort <-effort[year==2023]
# 
# non_functional_bpue <- bpue[tb.message == "random levels for at least one random effect not ok"]
# 
# # Step 2: Extract unique ecoregion, metierL4, and species combinations
# metier_species_combos <- unique(non_functional_bpue[, .(ecoregion, metierL4, species)])
# 
# # Step 3: Filter the monitor dataset based on these combinations
# filtered_monitor <- monitor[metier_species_combos, on = .(ecoregion, metierL4, species), nomatch = NULL]
# bpue <- fread("~/ICES WGBYC/WGBYC2024/BPUE_area.retain_2024_Sept19.csv")

calc_area_total <- function(bpue, effort, monitor, model_dir = "area_models/", verbose = TRUE) {
  
  require(ggeffects)
  require(data.table)
  require(glmmTMB)
  require(metafor)
  require(emmeans)
  
  # Convert bpue to a data.table if it isn't already
  if (!is.data.table(bpue)) bpue <- as.data.table(bpue)
  
  # Check for required columns
  required_columns <- c("areaCode", "metierL4", "species", "daysAtSea")
  missing_columns <- setdiff(required_columns, colnames(monitor))
  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in monitor:", paste(missing_columns, collapse = ", ")))
  }
  
  # Merge bpue data with monitor data and perform preprocessing
  if (verbose) print("Merging bpue and monitor datasets.")
  
  mon <- monitor[bpue, on = c("areaCode", "metierL4", "species"), allow.cartesian = TRUE]
  
  if (verbose) print("Inspecting merged monitor dataset.")
  print(head(mon))  # Check the first few rows of merged data
  
  mon <- mon[daysAtSea > 0]
  cols <- c("country", "year", "metierL5", "vesselLength_group", "samplingProtocol", "monitoringMethod")
  mon[, (cols) := lapply(.SD, as.factor), .SDcols = cols]
  mon$lDaS <- log(mon$daysAtSea)
  
  pred <- list(mean = NA_real_, lwr = NA_real_, upr = NA_real_, message = "OK", effort = NA_real_)
  
  # Check bpue$model for NA, "only one", or invalid cases
  if (verbose) print("Checking bpue$model.")
  
  # Ensure bpue$model is a single, non-empty character value
  if (length(bpue$model) != 1 || !is.character(bpue$model) || is.na(bpue$model) || bpue$model == "") {
    pred$message <- "bpue$model is NA, empty, or invalid"
    if (verbose) print(pred$message)
    return(as.data.table(pred))
  }
  
  if (bpue$model == "only one") {
    pred$message <- "bpue$model is 'only one'"
    if (verbose) print(pred$message)
    return(as.data.table(pred))
  }
  
  form <- as.formula(bpue$model)
  re <- lme4::findbars(form)  # Random effects part of model formulation (if any)
  re <- sapply(re, function(x) as.character(x[[3]]))
  re.n <- length(re)
  
  # Check for unsupported random effects
  if (re.n > 0 & any(c("samplingProtocol", "monitoringMethod") %in% re)) {
    pred$message <- "samplingProtocol or monitoringMethod not ok"
    if (verbose) print(pred$message)
    return(as.data.table(pred))
  }
  
  # Aggregate total effort
  if (verbose) print("Aggregating total effort.")
  tot <- effort[bpue, on = c("areaCode", "metierL4"), allow.cartesian = TRUE, .(das = sum(daysAtSeaF)), by = re]
  
  if (verbose) print("Inspecting aggregated total effort data.")
  print(head(tot))  # Check the first few rows of aggregated data
  
  tot$lDaS <- log(tot$das)
  
  # Check if all levels of random effects match between datasets
  if (any(sapply(re, function(r) !all(tot[[r]] %in% unique(mon[[r]]))))) {
    pred$message <- "random levels for at least one random effect not ok"
    if (verbose) print(pred$message)
    return(as.data.table(pred))
  }
  
  # Construct model file path
  if (verbose) print("Constructing model file path.")
  model_file <- file.path(model_dir, paste0("model_", bpue$areaCode[1], "_", bpue$metierL4[1], "_", bpue$species[1], ".rds"))
  
  # Load or fit model
  if (file.exists(model_file)) {
    if (verbose) print("Loading existing model.")
    best <- readRDS(model_file)
  } else {
    if (verbose) print("Model file not found. Re-fitting model.")
    best <- tryCatch({
      glmmTMB(formula = form, offset = lDaS, family = nbinom2, data = mon)
    }, error = function(e) {
      if (verbose) print(paste("Error fitting model:", e$message))
      NULL
    })
    if (!is.null(best)) {
      saveRDS(best, model_file) # Save the model for future use
    }
  }
  
  # Predictions
  if (is.null(best)) {
    pred$message <- "Model fitting failed"
    if (verbose) print(pred$message)
    return(as.data.table(pred))
  }
  
  if (re.n == 0) {
    tot <- effort[bpue, on = c("areaCode", "metierL4"), allow.cartesian = TRUE, .(das = sum(daysAtSeaF))]
    tot$lDaS <- log(tot$das)
    if (verbose) print("Generating predictions using emmeans.")
    pred_df <- as.data.frame(emmeans(best, ~1, offset = tot$lDaS, type = "response"))
    pred$mean <- mean(pred_df$response, na.rm = TRUE)
    pred$lwr <- mean(pred_df$asymp.LCL, na.rm = TRUE)
    pred$upr <- mean(pred_df$asymp.UCL, na.rm = TRUE)
    pred$effort <- sum(tot$das)
  } else {
    tot <- tot[complete.cases(tot)]
    if (verbose) print("Generating predictions with random effects.")
    pred_re <- lapply(1:nrow(tot), function(i) {
      p <- tryCatch({
        ggpredict(model = best,
                  terms = tot[i, ..re],
                  condition = c(lDaS = tot$lDaS[i]),
                  type = "re",
                  interval = "confidence")
      }, error = function(e) {
        if (verbose) print(paste("Error in ggpredict:", e$message))
        NULL
      })
      
      if (!is.null(p) && ncol(as.data.frame(p)) >= 3) {
        df_p <- as.data.frame(p)
        if (all(c("predicted", "conf.low", "conf.high") %in% colnames(df_p))) {
          return(df_p[, c("predicted", "conf.low", "conf.high")])
        } else {
          if (verbose) print("Columns missing in prediction result.")
          return(NULL)
        }
      } else {
        if (verbose) print("Prediction result is NULL or does not have enough columns.")
        return(NULL)
      }
    })
    
    pred_df <- do.call(rbind, pred_re)
    if (!is.null(pred_df) && nrow(pred_df) > 0) {
      pred$mean <- mean(pred_df$predicted, na.rm = TRUE)
      pred$lwr <- mean(pred_df$conf.low, na.rm = TRUE)
      pred$upr <- mean(pred_df$conf.high, na.rm = TRUE)
    }
    pred$effort <- sum(tot$das)
  }
  
  return(as.data.table(pred))
}


# 
# effort[vesselLengthRange %in% c(
#   "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
#   vesselLength_group := "below_12"]
# effort[vesselLengthRange %in% c(
#   "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
#   vesselLength_group := "above_12"]
# 
# for (i in 1:length(bpue$model)) {
#   if(bpue$model[i]>0) {
#     bpue[i,11:15]<-calc_area_total(bpue[i,],effort,filtered_monitor)
#     print(Sys.time())
#     flush.console()
#   }
#   print(Sys.time())
#   print(i)
#   flush.console()
# }
# Sys.time()-tic
# 
# a<-calc_area_total(BPUE_area.retain[],effort,filtered_monitor )
# a_filtered <- a[!is.na(mean) & !is.na(lwr) & !is.na(upr) & !is.infinite(lwr) & !is.infinite(upr)]
# write.csv(a_filtered,file("Area_tot_est_2023.csv"))
