# *******************************************************************
# *******************************************************************
# To calculate range-wide trends, calculate independent measures of relative
# abundance in each stratum.
# *******************************************************************
# *******************************************************************

my_packs <- c('tidyverse','sf','ebirdst','terra','scales','ggpubr') 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

library(bbsBayes)
library(ggrepel)

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

BBS_strata_boundaries <- sf::read_sf("C:/Users/IlesD/AppData/Local/R/win-library/4.3/bbsBayes/maps/BBS_USGS_strata.shp")

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
#blpw_ebird[blpw_ebird <= 0.0001] <- NA
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

# ------------------------------------------------
# USING BAM
# NOTE THAT BAM LACKS RELATIVE ABUNDANCE PREDICTIONS FOR AK.  
# USE BBS TO IMPUTE IT.
# ------------------------------------------------

# Load high resolution relative abundance raster from BAM
blpw_bam <- terra::rast("0_data/bam_density_raster/pred-BLPW-CAN-Mean.tif")
blpw_bam <- crop(blpw_bam, st_transform(strata_sf, st_crs(blpw_bam)))

# Convert from males/ha to males/km2
values(blpw_bam) <- values(blpw_bam)*100

# Stratum boundaries
BBS_strata_boundaries <- sf::read_sf("C:/Users/IlesD/AppData/Local/R/win-library/4.3/bbsBayes/maps/BBS_USGS_strata.shp") %>%
  st_intersection(st_union(strata_sf))

# Load BBS model
load("1_output/BBS/BBS_analysis.RData")

# Mean relative abundance in each region
stratum_indices <- generate_indices(jags_mod = jags_mod,
                                    jags_data = jags_data,
                                    regions = c("stratum"))$data_summary %>%
  dplyr::rename(ST_12 = Region) %>%
  group_by(ST_12) %>%
  summarize(Index = mean(Index,na.rm = TRUE))

BBS_strata_boundaries$area_km2 <- as.numeric(st_area(BBS_strata_boundaries))/1000000
BBS_strata_boundaries <- full_join(BBS_strata_boundaries, stratum_indices)

# bam estimate is total birds (birds per ha * 100 ha/km2)
comparison <- terra::extract(blpw_bam,BBS_strata_boundaries, bind = TRUE, fun = sum, na.rm = TRUE) %>%
  dplyr::rename(bam = pred.BLPW.CAN.Mean)  %>% 
  as.data.frame() %>%
  dplyr::select(ST_12,PROVSTATE,COUNTRY,area_km2,Index,bam)

comparison$Index[comparison$Index==0] <- NA

# Predictions of BAM's relative abundance in the missing strata (Alaska)
comparison$bam_area <- comparison$bam/comparison$area_km2
bam_lm <- lm(log(bam_area)~log(Index),data = subset(comparison, COUNTRY != "US"))

# Extract relative abundance in each stratum, fill in using prediction
blpw_bam_df <- terra::extract(blpw_bam,strata_sf)
colnames(blpw_bam_df)[2] <- "Density"
blpw_bam_df <- blpw_bam_df %>%
  group_by(ID) %>%
  summarize(Sum = sum(Density, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID])

# Add contribution from Alaska
comparison$pred <- exp(predict(bam_lm, newdata = comparison)) * comparison$area_km2

AK_contribution <- sum(subset(comparison,PROVSTATE == "AK")$pred,na.rm = TRUE)

blpw_bam_df$Sum[blpw_bam_df$Stratum == "West"] <- blpw_bam_df$Sum[blpw_bam_df$Stratum == "West"] + AK_contribution

# Western stratum was 2.43 times larger than eastern stratum in 2011
relabund <- rbind(relabund, data.frame(Source = "BAM",
                                       Stratum = c("East","West"),
                                       Sum = blpw_bam_df$Sum))



# Plot comparison between BBS and BAM
BAM_BBS_comparison_plot <- ggplot(comparison)+
  #geom_vline(data = subset(comparison, PROVSTATE == "AK"), aes(xintercept = Index), col = "dodgerblue", linetype = 2)+
  geom_smooth(data = subset(comparison, COUNTRY != "US"), aes(x = Index, y = bam/area_km2), method = "lm")+
  
  geom_point(data = subset(comparison, COUNTRY != "US"), aes(x = Index, y = bam/area_km2))+
  
  geom_point(data = subset(comparison, COUNTRY == "US" & PROVSTATE == "AK"), aes(x = Index, y = pred/area_km2), col = "dodgerblue", size = 3)+
  
  geom_text_repel(data = subset(comparison, COUNTRY != "US"),aes(x = Index, y = bam/area_km2, label = ST_12), size = 2, hjust = 1, vjust = 0)+
  geom_text_repel(data = subset(comparison, COUNTRY == "US" & PROVSTATE == "AK"),aes(x = Index, y = pred/area_km2, label = ST_12), size = 2, col = "dodgerblue", hjust = 1, vjust = 1)+
  
  
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  xlab("BBS (birds/route)")+
  ylab("BAM (birds/area)")+
  theme_bw()
BAM_BBS_comparison_plot

png(file ="1_output/Results_Appendix/BAM_BBS_comparison_plot.png", units = "in", width = 5, height = 5, res = 600)
BAM_BBS_comparison_plot
dev.off()

# ~~~
# Map of relative abundance
# ~~~

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
                       name = "males/km^2",
                       labels = comma)+
  scale_color_manual(values=strata_colours)+
  ggtitle("Relative Abundance (BAM 2011)")
BAM_Map


# ------------------------------------------------
# Combine maps
# ------------------------------------------------

maps_combined <- ggpubr::ggarrange(BAM_Map, eBird_Map, nrow = 2, labels = c("(a)","(b)"))

png(file ="1_output/Results_Appendix/Relative_Abundance_Maps.png", units = "in", width = 8, height = 10, res = 600)
maps_combined
dev.off()

# ------------------------------------------------
# Save relative abundance estimates
# ------------------------------------------------

write.csv(relabund, "1_output/Results_MainText/relabund_eBird_BAM.csv", row.names = FALSE)
