################################################################################
# Baltic Sea areacode-level BPUE for seals and harbour porpoise
################################################################################
# ============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

# ============================================================================
# 1. Annotate BPUE
# ============================================================================

# For cases where ecoregion-level estimation and prediction is considered 
# inappropriate, we can go to the areacode level.

obs3_area_BS_seals <- subset(obs3,(ecoregion=="baltic sea") &
                       (taxon=="mammals"))

# obs3_area_BS_seals <- subset(obs3,(ecoregion=="baltic sea") &
#                        (species =="pusa hispida" | species =="phocoena phocoena"))


#### SPLIT THE DATA WHERE ~IS USED FOR AREAS
# Count how many ICES areas each record contains
obs3_area_BS_seals <- obs3_area_BS_seals %>%
  mutate(
    n_areas = str_count(areacode, fixed("~")) + 1,
    daysatsea = daysatsea / n_areas
  ) %>%
  separate_rows(areacode, sep = "~") %>%
  select(-n_areas)

setDT(obs3_area_BS_seals)

obs3_area_BS_seals$daysatsea = round(obs3_area_BS_seals$daysatsea, 0)

### continue

#### Split also all2


#### SPLIT THE DATA WHERE ~IS USED FOR AREAS
# Count how many ICES areas each record contains

all2 <- all2 %>%
  mutate(
    n_areas = str_count(areacode, fixed("~")) + 1,
    daysatseaf = daysatseaf / n_areas
  ) %>%
  separate_rows(areacode, sep = "~") %>%
  select(-n_areas)


setDT(all2)

all2$daysatseaf = round(all2$daysatseaf, 0)


### cotninue to the estimate


area_m4_spec_BS_seals <- unique(obs3_area_BS_seals[, .(ecoregion, areacode, metierl4, species)])



bpue1_area_BS_seals <- calc_bpue(needle = area_m4_spec_BS_seals, 
                        cols = c("ecoregion", "areacode", "metierl4", "species"),
                        years = years, #added this
                        dat = obs3_area_BS_seals[taxon_bycatch_monitor_ok==TRUE])



# Optional export
# fwrite(
#   bpue1_area_BS_seals,
#   file = "data/bpue1_area_BS_seals.csv",
#   sep = ";",
#   na = "NA"
# )

bpue2_area_BS_seals <- annotate_bpue(bpue1_area_BS_seals, cols = c("ecoregion", "areacode", "metierl4"))

fwrite(bpue2_area_BS_seals, file = "data/bpue2_area_BS_seals.csv", sep = ";",na="NA")


#### zero assessment

#bpue1 <- fread("data/bpue1.csv")
#bpue2 <- fread("data/bpue2.csv")
source("lib/zeros_assessment.R")
zero_assessments(bpue2=bpue2_area_BS_seals,bpue1=bpue1_area_BS_seals)


# ============================================================================
# 2. MAKE into Bycatch estimate
# ============================================================================

total_bycatch_area_BS_seals <- vector(mode = "list", length = nrow(bpue1_area_BS_seals))

for (i in 1:nrow(bpue1_area_BS_seals)) {
    cat("\rProcessing row ", i, "/", nrow(bpue1_area_BS_seals), sep="")

    total_bycatch_area_BS_seals[[i]] <- calc_total(
        bpue = bpue1_area_BS_seals[i],
        cols = c("ecoregion", "areacode", "metierl4", "species"),
        obs = obs3_area_BS_seals[taxon_bycatch_monitor_ok==TRUE],
        allx = all2,
        years = years,
        verbose = FALSE
    )
}

total_bycatch_area_BS_seals <- rbindlist(total_bycatch_area_BS_seals)
#fwrite(total_bycatch_area_BS_seals, file = "data/tot_byc_area_BS_seals.csv", sep = ";",na="NA")

total_bycatch_area_BS_seals

#### Reliability test


library(glmmTMB)
library(ggeffects)
library(emmeans)
library(data.table)

##bpue_estimates
full_data <-obs3_area_BS_seals
bpues_estimates <-bpue1_area_BS_seals

source("lib/reliability_estimation_wkmamby.R")

tot_reliability <- total_bycatch_area_BS_seals

tot_reliability[bpues_estimates,
                on = .(ecoregion, metierl4, species),
                reliability_rmse := i.reliability_rmse]

tot_reliability$tot_lwr_log <- log10(tot_reliability$tot_lwr+1)
tot_reliability$tot_upr_log<-log10(tot_reliability$tot_upr+1)
tot_reliability$delta<-tot_reliability$tot_upr_log-tot_reliability$tot_lwr_log

tot_reliability$reliability_CI <- tot_reliability$delta < 2

monitoring_cutoff <- bpue2_area_BS_seals
monitoring_cutoff$monitoring_coverage <- monitoring_cutoff$daysatsea/monitoring_cutoff$daysatseaf

tot_reliability[monitoring_cutoff,
                on = .(ecoregion, metierl4, species),
                monitoring_coverage := i.monitoring_coverage]

tot_reliability$reliability_coverage <- tot_reliability$monitoring_coverage > 0.001

tot_reliability$overall_reliability <- tot_reliability$reliability_CI == TRUE & tot_reliability$reliability_rmse == TRUE & tot_reliability$reliability_coverage ==TRUE

#Need to bring in the zeros ?

tot_reliability[tot_mean == 0, overall_reliability := TRUE]

fwrite(tot_reliability, "results/tot_reliability_BS_seals.csv")


# ============================================================================
# 2. Define Baltic Sea areas
# ============================================================================

baltic_areas <- c(
  "27.3.c.22",
  "27.3.b.23",
  "27.3.d.24",
  "27.3.d.25",
  "27.3.d.26",
  "27.3.d.27",
  "27.3.d.28.1",
  "27.3.d.28.2",
  "27.3.d.29",
  "27.3.d.30",
  "27.3.d.31",
  "27.3.d.32"
)

species_to_map <- c(
  "pusa hispida",
  "phocoena phocoena",
  "halichoerus grypus",
  "phoca vitulina"
)

# ============================================================================
# ============================================================================
# 3. Create species x metier x area total-bycatch table
# ============================================================================

total_bycatch_table <- tot_reliability[overall_reliability == TRUE] %>%
  as.data.frame() %>%

  # Keep only the two target species
  filter(
    species %in% species_to_map
  ) %>%

  # Keep only Baltic Sea areas used in the maps
  filter(
    areacode %in% baltic_areas
  ) %>%

  # Keep only variables needed for tables and maps
  select(
    species,
    metierl4,
    areacode,
    tot_mean,
    tot_lwr,
    tot_upr,
    fishing_effort,
    message
  ) %>%

  # Combine possible repeated species-metier-area rows
  group_by(
    species,
    metierl4,
    areacode
  ) %>%

  summarise(
    tot_mean = if (
      all(is.na(tot_mean))
    ) {
      NA_real_
    } else {
      sum(tot_mean, na.rm = TRUE)
    },

    tot_lwr = if (
      all(is.na(tot_lwr))
    ) {
      NA_real_
    } else {
      sum(tot_lwr, na.rm = TRUE)
    },

    tot_upr = if (
      all(is.na(tot_upr))
    ) {
      NA_real_
    } else {
      sum(tot_upr, na.rm = TRUE)
    },

    fishing_effort = if (
      all(is.na(fishing_effort))
    ) {
      NA_real_
    } else {
      sum(fishing_effort, na.rm = TRUE)
    },

    message = paste(
      unique(message[!is.na(message)]),
      collapse = "; "
    ),

    .groups = "drop"
  ) %>%

  # Calculate absolute uncertainty interval width
  mutate(
    variation = tot_upr - tot_lwr
  ) %>%

  # Create every species-metier-area combination
  complete(
    species = species_to_map,
    metierl4,
    areacode = baltic_areas
  ) %>%

  arrange(
    species,
    metierl4,
    match(areacode, baltic_areas)
  )

# Inspect the complete table
print(
  total_bycatch_table,
  n = Inf
)

# Inspect only rows containing total-bycatch estimates
total_bycatch_table %>%
  filter(!is.na(tot_mean)) %>%
  print(n = Inf)

# Inspect messages other than OK
total_bycatch_table %>%
  filter(
    !is.na(message),
    message != "OK"
  ) %>%
  print(n = Inf)

# Optional export
# write.csv(
#   total_bycatch_table,
#   "total_bycatch_table_Baltic_Sea_species.csv",
#   row.names = FALSE
# )

# ============================================================================
# 4. Optional wide tables for inspection
# ============================================================================

total_bycatch_table_wide <- total_bycatch_table %>%
  select(
    species,
    metierl4,
    areacode,
    tot_mean
  ) %>%
  pivot_wider(
    names_from = areacode,
    values_from = tot_mean
  )

total_bycatch_variation_wide <- total_bycatch_table %>%
  select(
    species,
    metierl4,
    areacode,
    variation
  ) %>%
  pivot_wider(
    names_from = areacode,
    values_from = variation
  )

print(total_bycatch_table_wide)
print(total_bycatch_variation_wide)

# ============================================================================
# 5. Read and prepare Baltic Sea ICES polygons
# ============================================================================

ices <- st_read(
  paste0(
    "/Users/janiluke/Library/CloudStorage/",
    "OneDrive-Luonnonvarakeskus/",
    "VMS data/ICES_areas/",
    "ICES_Areas_20160601_cut_dense_3857.shp"
  ),
  quiet = TRUE
)

# Save the current spherical-geometry setting
previous_s2_setting <- sf_use_s2()

# Use planar geometry for repairing projected EPSG:3857 polygons
sf_use_s2(FALSE)

ices_baltic <- ices %>%

  # Retain only the required Baltic areas
  filter(
    Area_Full %in% baltic_areas
  ) %>%

  # Keep only the code and geometry
  select(
    areacode = Area_Full
  )

# Repair invalid geometries
ices_baltic <- st_make_valid(
  ices_baltic
)

cat(
  "Invalid ICES geometries after repair:",
  sum(!st_is_valid(ices_baltic)),
  "\n"
)

# Simplify in EPSG:3857; tolerance is in metres
ices_baltic <- st_simplify(
  ices_baltic,
  dTolerance = 1000,
  preserveTopology = TRUE
)

# Transform to longitude-latitude after repair and simplification
ices_baltic <- st_transform(
  ices_baltic,
  4326
)

# Restore the previous s2 setting
sf_use_s2(previous_s2_setting)

# ============================================================================
# 6. Join results to ICES polygons
# ============================================================================
# ============================================================================
# START REPLACEMENT: 6. Join total-bycatch estimates to ICES polygons
# ============================================================================

mapdat_total_bycatch <- ices_baltic %>%
  left_join(
    total_bycatch_table,
    by = "areacode"
  )

# Check values entering the maps
map_values_total_bycatch <- mapdat_total_bycatch %>%
  st_drop_geometry() %>%
  filter(!is.na(tot_mean)) %>%
  select(
    species,
    metierl4,
    areacode,
    tot_mean,
    tot_lwr,
    tot_upr,
    variation,
    fishing_effort,
    message
  ) %>%
  arrange(
    species,
    metierl4,
    match(areacode, baltic_areas)
  )


  map_values_total_bycatch

# ============================================================================
# 7. Natural Earth land background
# ============================================================================

land_baltic <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Exact cropping is unnecessary because coord_sf() defines the visible extent.
# Keeping the original layer also avoids geometry-intersection problems.

# ============================================================================
# 8. Helper function for map creation
# ============================================================================

create_species_map <- function(
    spatial_data,
    species_name,
    value_column,
    legend_title,
    plot_title,
    palette_option = "plasma",
    use_sqrt = TRUE
) {

  species_data <- spatial_data %>%
  filter(
    species == species_name
  ) %>%

  # Remove metier maps where the selected variable has only NA values
  group_by(
    metierl4
  ) %>%
  filter(
    any(!is.na(.data[[value_column]]))
  ) %>%
  ungroup()

  # Select the fill transformation
fill_scale <- scale_fill_gradientn(
  colours = c(
    "#edc9c0",
    "#d6846b",
    "#6b1208"
  ),
  na.value = "white",
  name = legend_title,
  trans = if (use_sqrt) "sqrt" else "identity"
)


  ggplot() +

    # Land background
    geom_sf(
      data = land_baltic,
      fill = "grey90",
      colour = "grey50",
      linewidth = 0.2
    ) +

    # ICES area results
    geom_sf(
      data = species_data,
      aes(
        fill = .data[[value_column]]
      ),
      colour = "black",
      linewidth = 0.3
    ) +

    # One panel for each metier
    facet_wrap(
      ~metierl4,
      ncol = 2
    ) +

    # Baltic Sea extent
    coord_sf(
      xlim = c(9, 33),
      ylim = c(53, 67),
      expand = FALSE,
      crs = st_crs(4326)
    ) +

    fill_scale +

    labs(
      title = plot_title
    ) +

    theme_bw() +

    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),

      strip.background = element_rect(
        fill = "grey95",
        colour = "grey60"
      ),

      strip.text = element_text(
        face = "bold"
      ),

      legend.position = c(0.82, 0.25) #or legend.position = "bottom"
    )
}

# ============================================================================
# 9. Create four graphs
# ============================================================================

# Graph 1: Pusa hispida BPUE
p_pusa_byc <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "pusa hispida",
  value_column = "tot_mean",
  legend_title = "Estimated total Bycatch",
  plot_title = "Pusa hispida: estimated total bycatch",
  palette_option = "plasma",
  use_sqrt = TRUE
)

# Graph 2: Pusa hispida uncertainty interval width
p_pusa_variation <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "pusa hispida",
  value_column = "variation",
  legend_title = "95 % CI width",
  plot_title = "Pusa hispida: uncertainty interval width",
  palette_option = "viridis",
  use_sqrt = TRUE
)

# Graph 3: Phocoena phocoena BPUE
p_phocoena_byc <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "phocoena phocoena",
  value_column = "tot_mean",
  legend_title = "Estimated total Bycatch",
  plot_title = "Phocoena phocoena: estimated total bycatch",
  palette_option = "plasma",
  use_sqrt = TRUE
)

# Graph 4: Phocoena phocoena uncertainty interval width
p_phocoena_variation <- create_species_map(
  spatial_data = mapdat,
  species_name = "phocoena phocoena",
  value_column = "variation",
  legend_title = "95 % CI width",
  plot_title = "Phocoena phocoena: uncertainty interval width",
  palette_option = "viridis",
  use_sqrt = TRUE
)

# Graph 5: Phoca vitulina BPUE
p_phoca_vitulina_byc <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "phoca vitulina",
  value_column = "tot_mean",
  legend_title = "Estimated total Bycatch",
  plot_title = "Phoca vitulina: estimated total bycatch",
  palette_option = "plasma",
  use_sqrt = TRUE
)

# Graph 6: Phoca vitulina uncertainty interval width
p_phoca_vitulina_variation <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "phoca vitulina",
  value_column = "variation",
  legend_title = "95 % CI width",
  plot_title = "Phoca vitulina: uncertainty interval width",
  palette_option = "viridis",
  use_sqrt = TRUE
)

# Graph 7: Halichoerus grypus BPUE
p_halichoerus_grypus_byc <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "halichoerus grypus",
  value_column = "tot_mean",
  legend_title = "Estimated total Bycatch",
  plot_title = "Halichoerus grypus: estimated total bycatch",
  palette_option = "plasma",
  use_sqrt = TRUE
)

# Graph 8: Halichoerus grypus uncertainty interval width
p_halichoerus_grypus_variation <- create_species_map(
  spatial_data = mapdat_total_bycatch,
  species_name = "halichoerus grypus",
  value_column = "variation",
  legend_title = "95 % CI width",
  plot_title = "Halichoerus grypus: uncertainty interval width",
  palette_option = "viridis",
  use_sqrt = TRUE
)

# ============================================================================
# 10. Display the four graphs
# ============================================================================

# print(p_pusa_byc)
# print(p_pusa_variation)
print(p_phocoena_byc)
# print(p_phocoena_variation)
print(p_phoca_vitulina_byc)
# print(p_phoca_vitulina_variation)
print(p_halichoerus_grypus_byc)
# print(p_halichoerus_grypus_variation)


# ============================================================================
# 11. Save total-bycatch maps
# ============================================================================
# ggsave(
#   filename = "total_bycatch_pusa_hispida.png",
#   plot = p_pusa_byc,
#   width = 12,
#   height = 8,
#   dpi = 300,
#   bg = "white"
# )

ggsave(
  filename = "total_bycatch_phoca_vitulina.png",
  plot = p_phoca_vitulina_byc,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "total_bycatch_phocoena_phocoena.png",
  plot = p_phocoena_byc,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ggsave(
ggsave(
  filename = "total_bycatch_halichoerus_grypus.png",
  plot = p_halichoerus_grypus_byc,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)


#### SAVE FINAL TABLE FOR REPORT

fwrite(tot_reliability, "tot_reliability_BS_mammals.csv", sep = ";", na = "NA")
