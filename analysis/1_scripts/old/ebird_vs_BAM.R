# --------------------------------------------
# Compare ebird and BAM
# --------------------------------------------

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse','sf','cowplot','bbsBayes','ebirdst','terra') 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)
library(ggsflabel) # custom library download from devtools::install_github("yutannihilation/ggsflabel")

#------------------------------------------------
# Clear working memory
#------------------------------------------------

rm(list=ls())
setwd("D:/Working_Files/1_Projects/Landbirds/BLPW-migration-trends/analysis/1_scripts")
setwd("../../")

# ------------------------------------------------
# Set ebirdst properties
# ------------------------------------------------

# Download eBird relative abundance map
# Need to set access key (only once) before downloading ranges
# ebirdst::set_ebirdst_access_key("ntqm1ha68fov",overwrite=TRUE)

# Ensure path is correctly set
#usethis::edit_r_environ()

# This should read:
# EBIRDST_KEY='ntqm1ha68fov'
# EBIRDST_DATA_DIR='D:/Working_Files/1_Projects/Landbirds/BLPW-migration-trends/analysis/0_data/eBird/'

# Download data
ebirdst_download("Blackpoll Warbler", pattern = "abundance_seasonal_mean_hr") 
path <- get_species_path("Blackpoll Warbler")

#------------------------------------------------
# Load strata shapefile; initially created using combine_strata.R
#------------------------------------------------

load("analysis/0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
nstrata <- nrow(strata_sf)
strata_sf$stratum_number <- 1:nstrata

# ------------------------------------------------
# Load BCR polygons
# ------------------------------------------------
bcr <- st_read("analysis/0_data/bcr/BCR_Terrestrial_master.shp") %>%
  subset(COUNTRY == "CANADA") %>%
  st_make_valid() %>%
  st_union

# ------------------------------------------------
# CALCULATE RELATIVE ABUNDANCES WITHIN EACH STRATUM BASED ON EBIRD
# ------------------------------------------------

# Load high resolution relative abundance ebirdraster
ebirdrast <- rast(paste0(path,"/seasonal/bkpwar_abundance_seasonal_mean_hr_2021.tif"))
ebirdrast <- ebirdrast[["breeding"]]

# Extract relative abundance in each stratum
ebirdrast <- crop(ebirdrast, st_transform(bcr, st_crs(ebirdrast))) %>%
  crop(.,st_transform(strata_sf, st_crs(ebirdrast)))

ebirdrast_df <- terra::extract(ebirdrast,st_transform(strata_sf,st_crs(ebirdrast))) %>%
  group_by(ID) %>%
  summarize(RelAbund = sum(breeding, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID]) 

# ------------------------------------------------
# CALCULATE RELATIVE ABUNDANCES WITHIN EACH STRATUM BASED ON BAM
# ------------------------------------------------

bamrast <- rast("analysis/0_data/bam_density_raster/pred-BLPW-CAN-Mean.tif")

bamrast <- crop(bamrast, st_transform(bcr, st_crs(bamrast)))
bamrast_df <- terra::extract(bamrast,st_transform(strata_sf,st_crs(bamrast))) %>%
  group_by(ID) %>%
  summarize(RelAbund = sum(`pred-BLPW-CAN-Mean`, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID])

ebirdrast <- project(ebirdrast,bamrast)

plot(bamrast)
lines(st_transform(strata_sf, st_crs(bamrast)))
plot(ebirdrast)
lines(st_transform(strata_sf, st_crs(bamrast)))
# ------------------------------------------------
# COMPARE RELATIVE ABUNDANCES
# ------------------------------------------------

ebird = ebirdrast_df$RelAbund/sum(ebirdrast_df$RelAbund)
BAM = bamrast_df$RelAbund/sum(bamrast_df$RelAbund)

relabund <- rbind(data.frame(Stratum = c("East","West"),
                             relabund = ebirdrast_df$RelAbund/sum(ebirdrast_df$RelAbund),
                             Source = "eBird"),
                  data.frame(Stratum = c("East","West"),
                             relabund = bamrast_df$RelAbund/sum(bamrast_df$RelAbund),
                             Source = "BAM"))


ggplot(relabund)+
  geom_bar(aes(x = Stratum, y = relabund, fill = Source), stat = "identity", position = position_dodge())+
  theme_bw()+
  ylab("Relative Abundance")+
  ggtitle("Relative Abundance in East and West Strata")
