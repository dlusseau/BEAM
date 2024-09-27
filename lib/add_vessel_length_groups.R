
add_vessel_length_groups <- function(x) {
    
    levels <- list(
        "below_12" = c("vl0006", "vl0008", "vl0010", "vl0015", "vl0608", "vl0612", "vl0810", "vl0815", "vl1012",
                       "VL0006", "VL0008", "VL0010", "VL0015", "VL0608", "VL0612", "VL0810", "VL0815", "VL1012"),
        "above_12" = c("vl1215", "vl1218", "vl1518", "vl15XX", "vl1824", "vl2440", "vl40xx", "vl40XX",
                      "VL1215", "VL1218", "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX")
    )
    
    x[grepl("^(f|nk)$", vessellengthrange, ignore.case = TRUE), vessellengthrange := NA_character_]
    x[, vessellength_group := factor(vessellengthrange, 
                                     levels = unlist(levels), 
                                     labels = rep(names(levels), lengths(levels)))]
    
}
