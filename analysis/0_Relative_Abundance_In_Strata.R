# *******************************************************************
# *******************************************************************
# To calculate range-wide trends, calculate independent measures of relative
# abundance in each stratum.
# *******************************************************************
# *******************************************************************

my_packs <- c('tidyverse','sf','ebirdst','terra','scales','ggpubr') 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# ------------------------------------------------
# Set graphical theme
# ------------------------------------------------

# Theme for plotting
CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA), 
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'transparent'),
                            panel.grid.minor = element_line(colour = 'transparent'),
                            
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",fill = NA),
                            axis.line = element_line(colour = "black"),
                            
                            strip.text = element_text(size = 12, colour = "black"),
                            
                            strip.background = element_rect(colour = "black",
                                                            fill = "gray95",
                                                            linetype = "solid"),
                            
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            
                            panel.background = element_rect(fill = "white"))

# Colors for plotting
strata_colours <- c("#016b9b","#D18A80")

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BLPW-migration-trends/analysis")

# ------------------------------------------------
# Object to store relative abundance estimates
# ------------------------------------------------

relabund <- data.frame()

#------------------------------------------------
# Load BBS strata shapefile
#------------------------------------------------

dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

#------------------------------------------------
# Load migration strata shapefile
#------------------------------------------------

load("0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf %>% st_transform(st_crs(BBS_strata_boundaries))
strata_sf$stratum_number <- 1:nrow(strata_sf)
strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))

# ------------------------------------------------
# USING EBIRD
# ------------------------------------------------

# Download relative abundance raster
#ebirdst_download_status("bkpwar",pattern = "abundance_seasonal_mean")

# Load high resolution relative abundance raster
blpw_ebird <- load_raster("bkpwar", 
                          product = "abundance", 
                          period = "seasonal",
                          metric = "mean",
                          resolution = "3km")

blpw_ebird <- blpw_ebird[["breeding"]]

# Extract relative abundance in each stratum
blpw_ebird <- crop(blpw_ebird, st_transform(strata_sf, st_crs(blpw_ebird)))
blpw_ebird_df <- terra::extract(blpw_ebird,strata_sf) %>%
  group_by(ID) %>%
  summarize(Sum = sum(breeding, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID])

relabund <- rbind(relabund, data.frame(Source = "eBird",
                                       Stratum = c("East","West"),
                                       Sum = blpw_ebird_df$Sum))

# ~~~
# Map of relative abundance
# ~~~
blpw_ebird[blpw_ebird <= 0.0001] <- NA
blpw_ebird <- mask(blpw_ebird, st_transform(strata_sf, st_crs(blpw_ebird)))

eBird_Map <- ggplot(data = strata_sf) +
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  tidyterra::geom_spatraster(data = blpw_ebird)+
  geom_sf(data = strata_sf, aes(col = Strata), fill = "transparent", linewidth = 1)+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5))+
  ggtitle("")+
  coord_sf(xlim = c(-3500000,2800000),
           ylim=c(500000,4500000))+
  scale_fill_gradientn(colors = c("transparent","lightblue","dodgerblue","darkblue"),
                       na.value = "transparent",
                       trans="log10",
                       name = "Relative Abundance",
                       labels = comma)+
  scale_color_manual(values=strata_colours)+
  ggtitle("Relative Abundance (eBird 2022)")
eBird_Map

# png(file ="2_output/Results_Appendix/eBird_relabund.png", units = "in", width = 8, height = 6, res = 600)
# eBird_Map
# dev.off()

# ------------------------------------------------
# USING BAM
# ------------------------------------------------

# Load high resolution relative abundance raster
blpw_bam <- terra::rast("0_data/bam_density_raster/pred-BLPW-CAN-Mean.tif")

# Extract relative abundance in each stratum
blpw_bam <- crop(blpw_bam, st_transform(strata_sf, st_crs(blpw_bam)))
blpw_bam_df <- terra::extract(blpw_bam,strata_sf)
colnames(blpw_bam_df)[2] <- "Density"
blpw_bam_df <- blpw_bam_df %>%
  group_by(ID) %>%
  summarize(Sum = sum(Density, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID])

# Eastern stratum was 1.45 times larger than western stratum in 2022
relabund <- rbind(relabund, data.frame(Source = "BAM",
                                       Stratum = c("East","West"),
                                       Sum = blpw_bam_df$Sum))

# ~~~
# Map of relative abundance
# ~~~
#blpw_bam[blpw_bam <= 0.001] <- NA
blpw_bam <- mask(blpw_bam, st_transform(strata_sf, st_crs(blpw_bam)))

BAM_Map <- ggplot(data = strata_sf) +
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  tidyterra::geom_spatraster(data = blpw_bam)+
  geom_sf(data = strata_sf, aes(col = Strata), fill = "transparent", linewidth = 1)+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5))+
  ggtitle("")+
  coord_sf(xlim = c(-3500000,2800000),
           ylim=c(500000,4500000))+
  scale_fill_gradientn(colors = c("transparent","lightblue","dodgerblue","darkblue"),
                       na.value = "transparent",
                       trans="log10",
                       name = "Relative Abundance",
                       labels = comma)+
  scale_color_manual(values=strata_colours)+
  ggtitle("Relative Abundance (BAM 2011)")
BAM_Map

# png(file ="2_output/Results_Appendix/BAM_relabund.png", units = "in", width = 8, height = 6, res = 600)
# BAM_Map
# dev.off()

# ------------------------------------------------
# Combine maps
# ------------------------------------------------

maps_combined <- ggpubr::ggarrange(eBird_Map, BAM_Map, nrow = 2, labels = c("(a)","(b)"))

png(file ="2_output/Results_Appendix/Relative_Abundance_Maps.png", units = "in", width = 8, height = 10, res = 600)
maps_combined
dev.off()

# ------------------------------------------------
# Save relative abundance estimates
# ------------------------------------------------

write.csv(relabund, "2_output/Results_MainText/relabund_eBird_BAM.csv", row.names = FALSE)
