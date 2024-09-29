# BEAM vX.Y

# imports all R packages that are needed for BEAM
library(data.table) # data wrangling
library(reshape2) # data wrangling
library(icesConnect) # downloading data from ICES servers
library(httr) # parsing downloaded data above
library(jsonlite) # parsing downloaded data
library(glmmTMB) # model fitting
library(metafor) # testing heterogeneity 
library(ggeffects) # BPUE estimation
library(emmeans) # BPUE estimation
library(doParallel) # parallel computation
library(stringdist) # data checking and cleaning (string manipulation)
#library(taxize) # data cleaning
library(taxizedb) # data cleaning, devtools::install_github("ropensci/taxizedb")

setwd("C:/Users/a21997/OneDrive - Havforskningsinstituttet/Prosjekter/git/BEAM")

# beam_lib.R imports all BEAM internal functionality
source("lib/beam_lib.R")

# uncomment to enable parallel execution
# pick a reasonable number of workers for your system
registerDoParallel(10)

# BEAM WORKFLOW

################################################################################
# STEP 1: Download WGBYC data from ICES server
# This gives us all1, obs1 and bycatch1 (and D4 and D5) in the global env (and saved on disk)
#beam_download(years = 2017:2023, token = "access_token_goes_here")
obs1 <- fread("C:/Users/a21997/OneDrive - Havforskningsinstituttet/Møter/2024-09 ICES WGBYC 2024/Data/D2_monitoringeffort_2017_2023.csv")
all1 <- fread("C:/Users/a21997/OneDrive - Havforskningsinstituttet/Møter/2024-09 ICES WGBYC 2024/Data/D1_fishingeffort_2017_2023.csv")
bycatch1 <- fread("C:/Users/a21997/OneDrive - Havforskningsinstituttet/Møter/2024-09 ICES WGBYC 2024/Data/D3_bycatchevent_2017_2023.csv")
D5 <- fread("C:/Users/a21997/OneDrive - Havforskningsinstituttet/Møter/2024-09 ICES WGBYC 2024/Data/D5_THELIST_roadmapspecies.csv")

################################################################################
# STEP 2: data cleaning and preparation
# gives us all2, obs2 and bycatch2 in global env (and saved on disk)
source("lib/clean_data.R")

################################################################################
# STEP 3: additional data cleaning and preparation
# generates obs3.csv (plus some additional files)
source("lib/generate_the_list.R")

################################################################################
# STEP 4: calculate BPUE estimates for ecoregion * metierl4 * species combos,
# based on the list created above.
eco_m4_spec <- unique(obs3[, .(ecoregion, metierl4, species)])
system.time(bpue1 <- calc_bpue(eco_m4_spec, dat = obs3)) # 237 secs per 1000 row with 10 workers
fwrite(bpue1, file = "data/bpue1.csv", sep = ";")

################################################################################

# Are there any other intermediary steps here?

################################################################################
# STEP 5: estimate total bycatch based on BPUE table
tot1 <- vector(mode = "list", length = nrow(bpue1))
for (i in 2313:nrow(bpue1)) {
    cat("\rProcessing row", i)
    tot1[[i]] <- calc_total(bpue1[i], obs = obs3, all = all2, verbose = FALSE)
}
tot1 <- rbindlist(tot1)
tot1 <- calc_total(bpue1[1:1000], obs = obs3, all = all2, cols = c("ecoregion", "metierl4", "species"))
fwrite(tot1, file = "data/tot1.csv", sep = ";")

################################################################################
# For any ecoregion * metierl4 * species combos where we can't predict total
# bycatch due to mismatches between r.e. levels in monitored and total effort data,
# try refitting those, but on the area * metierl4 * species level instead and 
# predict using those finer-scale models.
# Step 2: Extract unique ecoregion, metierL4, and species combinations
source("lib/add_finer_scale_models_when_needed.R")
bpue2 <- rbindlist(list(bpue1, bpue_fine_scale))
tot2 <- rbindlist(list(tot1, tot_fine_scale))

################################################################################
# final QC, assign color codes to models according to rules set out in script
source("lib/QCs_TotalBycatch.R")

################################################################################
# any final steps, e.g. final formatting, table generation and visualisation
fwrite(total_bycatch, "final_beam_table.csv")


######################################################


