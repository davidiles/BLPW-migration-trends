#------------------------------------------------
# Script to combine two "western" strata ("West" and "Central" into a single large western stratum)
#------------------------------------------------

#------------------------------------------------
# Load/install packages
#------------------------------------------------
my_packs <- c('tidyverse','RColorBrewer','sp','sf','rgdal','raster','janitor')
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

library(ggsflabel) # download from devtools::install_github("yutannihilation/ggsflabel")

# Projection
lcc <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

#------------------------------------------------
# Function to read/format breeding strata
#------------------------------------------------

setwd("~/iles_ECCC/Landbirds/BLPW-migration-trends-2022_2strata/analysis/0_data/strata")

strata_crs = "+proj=longlat +datum=WGS84 +no_defs "

strata_sf <- read_sf(dsn = "BLPW_Strata2.shp") %>%
  st_set_crs(st_crs(crs(strata_crs)))

# Combine western strata into a single one
sf::sf_use_s2(FALSE)
strata_eastern = strata_sf[3,] %>%
  mutate(Strata = "East")
strata_western = strata_sf[1:2,] %>%
  summarise(geometry = st_union(geometry)) %>%
  add_column(Strata = "West")

strata_sf = dplyr::bind_rows(strata_eastern,strata_western)
strata_sf$Strata <- factor(strata_sf$Strata)
plot(strata_sf)

# For context
BCR_poly <- st_read("../bcr/BCR_Terrestrial_master.shp") %>%
  subset(., !(PROVINCE_S %in% c("HAWAIIAN ISLANDS")) & COUNTRY %in% c("USA","CANADA")) %>%
  st_transform(st_crs(strata_sf))

strata <- list(strata_sf = strata_sf,BCR_poly = BCR_poly)

save(strata, file = "strata_East_West.RData")
