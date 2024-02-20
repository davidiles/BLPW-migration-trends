# *******************************************************************
# *******************************************************************
# Summary of Main Results for Publication
# Note that figures in appendices are produced in earlier model-specific scripts (e.g., goodness-of-fit evaluation)
# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse','sf','cowplot','bbsBayes','ebirdst','terra','scales') 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)
library(ggsflabel) # custom library download from devtools::install_github("yutannihilation/ggsflabel")

#------------------------------------------------
# Clear working memory
#------------------------------------------------

rm(list=ls())

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BLPW-migration-trends/analysis")

# *******************************************************************
# *******************************************************************
# Set custom theme and load additional shapefiles
# *******************************************************************
# *******************************************************************

`%!in%` <- Negate(`%in%`)

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

#------------------------------------------------
# Load migration strata shapefile
#------------------------------------------------

load("0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
strata_sf$stratum_number <- 1:nrow(strata_sf)
strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))

#------------------------------------------------
# Load BBS strata shapefile
#------------------------------------------------

dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

# *******************************************************************
# *******************************************************************
# LOAD FITTED MODELS 
# *******************************************************************
# *******************************************************************

# ---------------------------------
# Pre-breeding migration
# ---------------------------------

load("1_output/Spring/wksp_Spring.RData")
count_df_Spring <- count_df
station_data_summarized_sf_Spring <- station_data_summarized_sf
jags_data_Spring <- jags_data
out_Spring <- out

# ---------------------------------
# Post-breeding migration
# ---------------------------------

load("1_output/Fall/wksp_Fall.RData")
count_df_Fall <- count_df
station_data_summarized_sf_Fall <- station_data_summarized_sf
jags_data_Fall <- jags_data
out_Fall <- out

# ---------------------------------
# BBS
# ---------------------------------

load(file = "1_output/BBS/BBS_analysis.RData")
jags_data_BBS <- jags_data
out_BBS <- jags_mod
strat_data_BBS <- strat_data

rm(out)

# *******************************************************************
# *******************************************************************
# PREPARE FIGURE 1: MAP OF SURVEY LOCATIONS
# *******************************************************************
# *******************************************************************

# Convert to AEA 
strata_sf <- strata_sf %>% st_transform(st_crs(BBS_strata_boundaries))

all_stations <- c(station_data_summarized_sf_Spring$station,
                  station_data_summarized_sf_Fall$station) %>% unique()
both_seasons <- all_stations[all_stations %in% station_data_summarized_sf_Spring$station &
                               all_stations %in% station_data_summarized_sf_Fall$station]
Spring_only <- station_data_summarized_sf_Spring$station[station_data_summarized_sf_Spring$station %!in% both_seasons]
Fall_only <- station_data_summarized_sf_Fall$station[station_data_summarized_sf_Fall$station %!in% both_seasons]


all_stations <- bind_rows(station_data_summarized_sf_Fall,station_data_summarized_sf_Spring) %>%
  dplyr::select(station) %>% unique()

all_stations$Season <- "Both"
all_stations$Season[all_stations$station %in% Spring_only] <- "Spring"
all_stations$Season[all_stations$station %in% Fall_only] <- "Fall"
all_stations$Season <- factor(all_stations$Season, levels = c("Spring","Both","Fall"),
                              labels = c("Pre-breeding","Both","Post-breeding"))

Season_Cols <- RColorBrewer::brewer.pal(10,"PRGn")[c(3,8)]
map_combined <- ggplot(data = strata_sf) +
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7, col = "transparent")+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5))+
  geom_sf(data = all_stations,size = 1.5, alpha = 1, shape = 19, col = "black")+
  geom_sf(data = all_stations,aes(col = Season),size = 1, alpha = 1, shape = 19)+
  
  geom_sf_label_repel(data = all_stations,
                      aes(label = station, col = Season),
                      max.overlaps = 50,
                      force = 1.5,
                      force_pull = 0.1,
                      size = 3, 
                      alpha = 1,
                      min.segment.length = 0,
                      fontface = "bold",
                      show_guide = FALSE)+
  scale_color_manual(values=c(Season_Cols[1],"black",Season_Cols[2]))+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = "none") +
  ggtitle("")+
  coord_sf(xlim = c(-3500000,2800000),
           ylim=c(-500000,4200000))
map_combined

png(file ="1_output/Results_MainText/Figure1.png", units = "in", width = 8, height = 6, res = 600)
map_combined
dev.off()

# *******************************************************************
# *******************************************************************
# Load annual indices from migration monitoring, prepared by previous scripts
# *******************************************************************
# *******************************************************************

strata_colours <- c("#016b9b","#D18A80","black")

Spring_results <- readRDS("1_output/Spring/results_Spring.rds")
Spring_indices <- Spring_results$indices_summarized

Fall_results <- readRDS("1_output/Fall/results_Fall.rds")
Fall_indices <- Fall_results$indices_summarized

BBS_results <- readRDS("1_output/BBS/results_BBS.rds")
BBS_indices <- BBS_results$indices_summarized

ylim <- range(c(Spring_indices[,c("indices_rescaled_q025","indices_rescaled_q975")],Fall_indices[,c("indices_rescaled_q025","indices_rescaled_q975")]))

# ~~~~~~~~~~~~
# Pre-breeding migration indices
# ~~~~~~~~~~~~


Spring_indices$Stratum <- factor(Spring_indices$Stratum, levels = c("West","East","Continental"))

# Plot indices for pre-breeding migration
Fig2_Spring <- ggplot(data = Spring_indices, 
                      aes(x = Year,
                          y = indices_rescaled_q500, 
                          ymin = indices_rescaled_q025, 
                          ymax = indices_rescaled_q975,
                          col = Stratum),
)+
  geom_errorbar(width = 0)+
  geom_point()+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=ylim)+
  scale_color_manual(values=strata_colours,guide="none")+
  ylab("Population Index")+
  facet_grid(.~Stratum)

Fig2_Spring

# ~~~~~~~~~~~~
# Pre-breeding migration indices
# ~~~~~~~~~~~~

Fall_indices$Stratum <- factor(Fall_indices$Stratum, levels = c("West","East","Continental"))

Fig2_Fall <- ggplot(data = Fall_indices, 
                    aes(x = Year,
                        y = indices_rescaled_q500, 
                        ymin = indices_rescaled_q025, 
                        ymax = indices_rescaled_q975,
                        col = Stratum),
)+
  geom_errorbar(width = 0)+
  geom_point()+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=ylim)+
  scale_color_manual(values=strata_colours,guide="none")+
  ylab("Population Index")+
  facet_grid(.~Stratum)

Fig2_Fall

# ~~~~~~~~~~~~
# BBS
# ~~~~~~~~~~~~

BBS_indices$Stratum <- factor(BBS_indices$Stratum, levels = c("West","East","Continental"))

Fig2_BBS <- ggplot(data = BBS_indices, 
                    aes(x = Year,
                        y = indices_rescaled_q500, 
                        ymin = indices_rescaled_q025, 
                        ymax = indices_rescaled_q975,
                        col = Stratum),
)+
  geom_errorbar(width = 0)+
  geom_point()+
  scale_y_continuous(trans="log10")+
  coord_cartesian(ylim=ylim)+
  scale_color_manual(values=strata_colours,guide="none")+
  ylab("Population Index")+
  facet_grid(.~Stratum)

Fig2_BBS
