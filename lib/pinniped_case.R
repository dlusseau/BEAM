# Pinniped-specific preparation for BEAM
#
# This file does not replace or modify BEAM model estimation.
# It only constructs three endpoints from the standard BEAM obs3/bycatch3 objects:
#   1) halichoerus grypus (pure grey seal)
#   2) phoca vitulina (pure harbour seal)
#   3) grey_harbour (grey + harbour + unresolved Phocidae/Pinnipedia)
#
# The pooled endpoint keeps one copy of the D2 monitoring denominator per BEAM
# monitoring cell and pools the D3 numerator at the same grain.

prepare_pinniped_case <- function(
    obs3_standard,
    bycatch3_standard,
    obs2,
    bycatch2,
    all2,
    results_dir = "results") {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  obs3_standard <- data.table::copy(obs3_standard)
  bycatch3_standard <- data.table::copy(bycatch3_standard)
  obs2 <- data.table::copy(obs2)
  bycatch2 <- data.table::copy(bycatch2)
  all2 <- data.table::copy(all2)

  # Exact single ICES area codes only. Composite codes such as
  # 27.3.d.24~27.3.d.25 are audited rather than split between areas.
  ices_area_keys <- unique(
    obs2[
      areatype == "icesarea" &
        !is.na(areacode) &
        !grepl("~", areacode, fixed = TRUE),
      .(ecoregion, areacode)
    ]
  )

  combined_area_monitoring <- obs2[
    areatype == "icesarea" &
      !is.na(areacode) &
      grepl("~", areacode, fixed = TRUE)
  ]
  combined_area_bycatch <- bycatch2[
    areatype == "icesarea" &
      !is.na(areacode) &
      grepl("~", areacode, fixed = TRUE)
  ]
  combined_area_fishing <- all2[
    areatype == "icesarea" &
      !is.na(areacode) &
      grepl("~", areacode, fixed = TRUE)
  ]

  data.table::fwrite(
    combined_area_monitoring,
    file.path(results_dir, "pinniped_case_combined_area_monitoring_audit.csv"),
    sep = ";", na = "NA"
  )
  data.table::fwrite(
    combined_area_bycatch,
    file.path(results_dir, "pinniped_case_combined_area_bycatch_audit.csv"),
    sep = ";", na = "NA"
  )
  data.table::fwrite(
    combined_area_fishing,
    file.path(results_dir, "pinniped_case_combined_area_fishing_audit.csv"),
    sep = ";", na = "NA"
  )

  pure_species <- c("halichoerus grypus", "phoca vitulina")
  unresolved_species <- c("phocidae", "pinnipedia")
  combined_sources <- c(pure_species, unresolved_species)

  # Pure species rows remain the standard BEAM rows.
  obs3_pure <- obs3_standard[
    species %chin% pure_species
  ][
    ices_area_keys,
    on = .(ecoregion, areacode),
    nomatch = 0
  ]

  # Unresolved pinnipeds are included in the combined endpoint but are also
  # exported separately so their contribution remains visible.
  unresolved_detail <- bycatch3_standard[
    species %chin% unresolved_species
  ][
    ices_area_keys,
    on = .(ecoregion, areacode),
    nomatch = 0
  ]

  unresolved_area_summary <- unresolved_detail[
    , .(
      n_records = .N,
      n_unclassified = sum(n_ind, na.rm = TRUE)
    ),
    by = .(ecoregion, areacode, year, metierl4, species)
  ][order(year, ecoregion, areacode, metierl4, species)]

  unresolved_totals <- unresolved_detail[
    , .(
      n_records = .N,
      n_unclassified = sum(n_ind, na.rm = TRUE)
    ),
    by = species
  ]

  data.table::fwrite(
    unresolved_detail,
    file.path(results_dir, "pinniped_case_unclassified_detail.csv"),
    sep = ";", na = "NA"
  )
  data.table::fwrite(
    unresolved_area_summary,
    file.path(results_dir, "pinniped_case_unclassified_area_summary.csv"),
    sep = ";", na = "NA"
  )
  data.table::fwrite(
    unresolved_totals,
    file.path(results_dir, "pinniped_case_unclassified_totals.csv"),
    sep = ";", na = "NA"
  )

  # This is the aggregation grain used to construct a valid pooled monitoring
  # row before handing the data back to the standard BEAM modelling functions.
  monitor_grain <- c(
    "ecoregion",
    "areacode",
    "country",
    "year",
    "quarter",
    "metierl4",
    "metierl5",
    "vessellength_group",
    "samplingprotocol",
    "monitoringmethod"
  )

  # Standard generate_the_list.R gives grey and harbour rows copies of the
  # same D2 monitoring denominator. Check that before deduplicating them.
  monitoring_denominator_qc <- obs3_pure[
    , .(
      n_distinct_daysatsea = uniqueN(daysatsea),
      min_daysatsea = suppressWarnings(min(daysatsea, na.rm = TRUE)),
      max_daysatsea = suppressWarnings(max(daysatsea, na.rm = TRUE))
    ),
    by = monitor_grain
  ][n_distinct_daysatsea > 1]

  if (nrow(monitoring_denominator_qc) > 0) {
    stop(
      "Grey and harbour rows do not have identical D2 monitoring denominators in ",
      nrow(monitoring_denominator_qc),
      " monitoring cells. Inspect $monitoring_denominator_qc before pooling."
    )
  }

  pooled_monitoring <- obs3_pure[
    , .(
      daysatsea = max(daysatsea, na.rm = TRUE),
      taxa_monitored = taxa_monitored[1],
      taxon_bycatch_monitor_ok = any(
        taxon_bycatch_monitor_ok %in% TRUE,
        na.rm = TRUE
      )
    ),
    by = monitor_grain
  ]

  pooled_monitoring[, `:=`(
    species = "grey_harbour",
    taxon = "mammals"
  )]

  # Pool the D3 numerator at exactly the same grain.
  pooled_bycatch <- bycatch3_standard[
    species %chin% combined_sources
  ][
    ices_area_keys,
    on = .(ecoregion, areacode),
    nomatch = 0
  ][
    , .(n_ind = sum(n_ind, na.rm = TRUE)),
    by = monitor_grain
  ]

  # Positive D3 cells without a corresponding D2 denominator are audited and
  # are not silently forced into BPUE estimation.
  unmatched_pooled_bycatch <- pooled_bycatch[
    !pooled_monitoring,
    on = monitor_grain
  ][n_ind > 0]

  data.table::fwrite(
    unmatched_pooled_bycatch,
    file.path(results_dir, "pinniped_case_unmatched_pooled_bycatch.csv"),
    sep = ";", na = "NA"
  )

  pooled_obs3 <- data.table::copy(pooled_monitoring)
  pooled_obs3[, n_ind := 0]
  pooled_obs3[
    pooled_bycatch,
    on = monitor_grain,
    n_ind := i.n_ind
  ]
  pooled_obs3[, n_ind := round(n_ind, 0)]

  obs3_pinniped <- data.table::rbindlist(
    list(obs3_pure, pooled_obs3),
    use.names = TRUE,
    fill = TRUE
  )

  # D1 effort remains species independent and is used once for each endpoint.
  all2_pinniped <- all2[
    areatype == "icesarea" &
      !is.na(areacode) &
      !grepl("~", areacode, fixed = TRUE)
  ]

  source_counts <- bycatch3_standard[
    species %chin% combined_sources
  ][
    ices_area_keys,
    on = .(ecoregion, areacode),
    nomatch = 0
  ][
    , .(n = sum(n_ind, na.rm = TRUE)),
    by = species
  ]

  count_qc <- data.table::data.table(
    source_count = source_counts[, sum(n, na.rm = TRUE)],
    pooled_count = pooled_obs3[, sum(n_ind, na.rm = TRUE)]
  )

  list(
    obs3 = obs3_pinniped,
    all2 = all2_pinniped,
    unresolved_detail = unresolved_detail,
    unresolved_area_summary = unresolved_area_summary,
    unresolved_totals = unresolved_totals,
    unmatched_pooled_bycatch = unmatched_pooled_bycatch,
    monitoring_denominator_qc = monitoring_denominator_qc,
    count_qc = count_qc,
    ices_area_keys = ices_area_keys
  )
}


# Add descriptive diagnostics to BEAM BPUE output without changing any model,
# BPUE, zero-assessment, or raising decision.
add_pinniped_diagnostics <- function(
    bpue,
    obs3_pinniped,
    cols = c("ecoregion", "areacode", "metierl4", "species")) {

  bpue <- data.table::copy(bpue)
  obs3_pinniped <- data.table::copy(obs3_pinniped)

  input_summary <- obs3_pinniped[
    taxon_bycatch_monitor_ok %in% TRUE,
    .(
      observed_n = sum(n_ind, na.rm = TRUE),
      monitored_days = sum(daysatsea, na.rm = TRUE),
      n_monitoring_cells = .N,
      years_with_monitoring = uniqueN(year[!is.na(year)])
    ),
    by = cols
  ]

  out <- merge(bpue, input_summary, by = cols, all.x = TRUE)

  if (!"bpue" %in% names(out)) out[, bpue := NA_real_]
  if (!"model" %in% names(out)) out[, model := NA_character_]

  out[, result_status := data.table::fcase(
    observed_n == 0 & !is.na(bpue) & bpue == 0,
      "zero observed; BPUE=0 returned by BEAM",
    observed_n == 0 & (is.na(bpue) | model == "none"),
      "zero observed; no BPUE model retained",
    observed_n > 0 & !is.na(bpue),
      "BPUE estimated",
    observed_n > 0 & is.na(bpue),
      "positive bycatch; BPUE unavailable",
    default = "review"
  )]

  out[]
}


# Add a descriptive raising status after calc_total(). This does not modify the
# value returned by calc_total().
add_pinniped_total_status <- function(results) {
  results <- data.table::copy(results)

  if (!"tot_mean" %in% names(results)) results[, tot_mean := NA_real_]
  if (!"observed_n" %in% names(results)) results[, observed_n := NA_real_]

  results[, total_status := data.table::fcase(
    !is.na(tot_mean),
      "raised estimate available",
    is.na(tot_mean) & observed_n == 0,
      "no raised estimate; zero-observation case",
    is.na(tot_mean) & observed_n > 0,
      "no raised estimate despite positive bycatch",
    default = "no raised estimate"
  )]

  results[]
}
