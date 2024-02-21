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
# Summarize and plot estimates of annual indices from each program
#  - rescale annual indices based on relative abundance
# *******************************************************************
# *******************************************************************

# ---------------------------------
# LOAD RELATIVE ABUNDANCE ESTIMATES
# ---------------------------------

relabund <- read.csv("1_output/Results_MainText/relabund_eBird_BAM.csv")

source_of_estimate <- "BAM"
relabund <- subset(relabund, Source == source_of_estimate)
relabund$Sum <- relabund$Sum/relabund$Sum[1] # West stratum will have relative abundance of 1

# Year to which relative abundance will be "pinned"
if (source_of_estimate == "eBird") relabund_year <- 2018
if (source_of_estimate == "BAM") relabund_year <- 2011

# ---------------------------------
# Load indices from migration monitoring
# ---------------------------------

# Full posterior sample
Spring_indices <- readRDS("1_output/Spring/results_Spring.rds")$indices_df %>% mutate(Season = "Spring")
Fall_indices <- readRDS("1_output/Fall/results_Fall.rds")$indices_df %>% mutate(Season = "Fall")

mig_indices <- bind_rows(Spring_indices,Fall_indices) %>% subset(Stratum %in% c("East","West"))
mig_indices$row <- 1:nrow(mig_indices)

# ---------------------------------
# Rescale annual indices based on relative abundance in defined year (year of relative abundance map)
# ---------------------------------

# Calculate annual means of each index (used for rescaling based on relabund_year above)
annual_means <- mig_indices %>%
  group_by(Season, Year,year_number,Stratum,stratum_number) %>%
  summarize(index_mean = mean(indices))

relabund_mean <- subset(annual_means, Year %in% relabund_year)

# Create new column containing the mean index in the year of the relative abundance estimate
mig_indices$mean_index_in_relabund_year <- NA
mig_indices$rescaling_factor <- NA
mig_indices$indices_rescaled2 <- NA

for (season in unique(mig_indices$Season)){
  
  for (stratum in unique(mig_indices$Stratum)){
    
    # Relevant rows in dataframe
    rows <- subset(mig_indices, Stratum == stratum & Season == season)$row
    
    # Mean index in year of relative abundance estimate
    mig_indices$mean_index_in_relabund_year[rows] <- subset(annual_means, Year == relabund_year & Stratum == stratum & Season == season)$index_mean
    
    # Factor by which to rescale estimates
    mig_indices$rescaling_factor[rows] <- subset(relabund,Stratum == stratum)$Sum
    
  }
}

# Rescale indices
mig_indices$indices_rescaled <- mig_indices$indices/mig_indices$mean_index_in_relabund_year * mig_indices$rescaling_factor

# ---------------------------------
# Calculate range-wide indices by summing rescaled indices across strata for each sample from posterior, in each year
# ---------------------------------

continental_indices <- mig_indices %>%
  group_by(sample,Year,Season) %>%
  summarize(indices_rescaled = sum(indices_rescaled)) %>%
  mutate(Stratum = "Continental")

# append to mig_indices
mig_indices <- bind_rows(mig_indices, continental_indices)

# ---------------------------------
# Load BBS indices
# ---------------------------------

BBS_indices <- readRDS("1_output/BBS/results_BBS.rds")$BBS_indices_df %>% 
  mutate(Source = "Breeding Bird Survey") %>% 
  dplyr::rename(indices_rescaled = indices)

# ---------------------------------
# Combine all estimates
# ---------------------------------

mig_indices$Source[mig_indices$Season == "Fall"] <- "Post-breeding migration"
mig_indices$Source[mig_indices$Season == "Spring"] <- "Pre-breeding migration"

all_indices <- bind_rows(mig_indices,BBS_indices)
all_indices$Source <- factor(all_indices$Source, levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"))
all_indices$Stratum <- factor(all_indices$Stratum, levels = c("West","East","Continental"))

# ---------------------------------
# Calculate trends and percent change for each sample from the posterior
# ---------------------------------

all_indices <- all_indices %>%
  dplyr::select(sample,Stratum,Source,Year,indices_rescaled) %>%
  group_by(sample,Stratum,Source) %>%
  mutate(log_change = log(indices_rescaled/indices_rescaled[Year == 1998]),
         percent_change = 100*(indices_rescaled-indices_rescaled[Year == 1998])/indices_rescaled[Year == 1998],
         trend = 100*((indices_rescaled/indices_rescaled[Year == 1998])^(1/(2018-1998))-1))

# Convert log change to percent change using 100*(exp(log_change)-1)

# ---------------------------------
# Summarize annual indices (for plotting trajectories)
# ---------------------------------

analysis_summary <- all_indices %>%
  group_by(Year,Stratum,Source) %>%
  summarize(log_change_q025 = quantile(log_change,0.025) ,
            log_change_q500 = quantile(log_change,0.500),
            log_change_q975 = quantile(log_change,0.975),
            
            percent_change_q025 = quantile(percent_change,0.025) %>% round(),
            percent_change_q500 = quantile(percent_change,0.500) %>% round(),
            percent_change_q975 = quantile(percent_change,0.975) %>% round(),
            
            trend_q025 = quantile(trend,0.025) %>% round(1),
            trend_q500 = quantile(trend,0.500) %>% round(1),
            trend_q975 = quantile(trend,0.975) %>% round(1),
            
            prob_positive = round(mean(trend > 0),2)
  )

# Summary in final year (represents overall trend from 1998 to 2008)

summary_2018 <- subset(analysis_summary, Year == 2018) %>%
  dplyr::select(-log_change_q025,-log_change_q500,-log_change_q975) %>%
  as.data.frame()

summary_2018

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

strata_colours <- c("#016b9b","#D18A80","black")


# ------------
# Prepare y axis scale (conversion between log-scale change and percent change) 
# Convert log-scale change to percent change using: 100*(exp(log_change)-1)
# ------------

y_axis_breaks <- c(0,50,100,200,300)
log_breaks <- log(y_axis_breaks/100+1)
log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
y_axis_breaks <- 100*(exp(log_breaks)-1)

y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
y_axis_labels <- paste0(y_axis_labels," %")

# ------------
# Text labels for each panel
# ------------

text_labels <- expand.grid(Year = 1998,
                           log_y = max(log_breaks),
                           Stratum = factor(c("West","East","Continental"),levels = c("West","East","Continental")),
                           Source = factor(c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"),levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey")))
text_labels$panel <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")

# ------------
# Plot itself
# ------------

Fig2 <- ggplot()+
  geom_ribbon(data = analysis_summary,
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975, 
                  fill = Stratum, 
                  col = Stratum),alpha = 0.2, col = "transparent")+
  geom_line(data = analysis_summary,
            aes(x = Year, 
                y = log_change_q500, 
                col = Stratum),linewidth = 2)+
  geom_text(data = text_labels, aes(x = Year, y = log_y, label = panel),hjust=0, size = 6,vjust=1)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_grid(Source~Stratum)+
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")


png(file = paste0("1_output/Results_MainText/Figure2_",source_of_estimate,".png"), units = "in", width = 7, height = 7, res = 600)
Fig2
dev.off()

