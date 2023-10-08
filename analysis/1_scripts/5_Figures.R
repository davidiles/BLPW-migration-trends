# *******************************************************************
# *******************************************************************
# Figures / Main Results for publication
# *******************************************************************
# *******************************************************************

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

# ------------------------------------------------
# Set working directory
# ------------------------------------------------

stub <- function() {}
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

dirname <- thisPath()

setwd(dirname)
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

# *******************************************************************
# *******************************************************************
# Load fitted models for Spring and Fall
# *******************************************************************
# *******************************************************************

# Spring
load("analysis/2_output/Spring/analysis_Spring.RData")
count_df_Spring <- count_df
station_data_summarized_sf_Spring <- station_data_summarized_sf
jags_data_Spring <- jags_data
out_Spring <- out

# Fall
load("analysis/2_output/Fall/analysis_Fall.RData")
count_df_Fall <- count_df
station_data_summarized_sf_Fall <- station_data_summarized_sf
jags_data_Fall <- jags_data
out_Fall <- out

# BBS analysis
load(file = "analysis/2_output/BBS/BBS_analysis.RData")
jags_data_BBS <- jags_data
out_BBS <- jags_mod
strat_data_BBS <- strat_data

# Only retain relevant objects
rm(list=ls()[! ls() %in% c("count_df_Spring","station_data_summarized_sf_Spring","jags_data_Spring","out_Spring",
                           "count_df_Fall", "station_data_summarized_sf_Fall","jags_data_Fall","out_Fall",
                           "jags_data_BBS","out_BBS","strat_data_BBS",
                           "dirname","path")])


# *******************************************************************
# *******************************************************************
# Set custom theme and load additional shapefiles
# *******************************************************************
# *******************************************************************

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

load("analysis/0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
strata_sf$stratum_number <- 1:nrow(strata_sf)
strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))

#------------------------------------------------
# Load BBS strata shapefile
#------------------------------------------------

dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

#------------------------------------------------
# Define custom 'east' and 'west' strata for BBS analysis
#------------------------------------------------

stratum_indices <- generate_indices(jags_mod = out_BBS,
                                    jags_data = jags_data_BBS,
                                    regions = c("stratum"))

stratum_trends <- generate_trends(indices = stratum_indices,
                                  Min_year = 2008,
                                  Max_year = 2018)

# Only include strata that were included in analysis
BBS_strata_boundaries <- subset(BBS_strata_boundaries, ST_12 %in% stratum_trends$Strata_included)

# Overlap between migration strata and BBS strata
strata_sf <- st_transform(strata_sf, st_crs(BBS_strata_boundaries))
BBS_mig_strata = st_intersects(BBS_strata_boundaries,strata_sf) %>% as.numeric()
BBS_mig_strata[is.na(BBS_mig_strata)] <- 1
BBS_mig_strata <- c("East","West")[BBS_mig_strata] %>% factor(levels = c("West","East"))
BBS_strata_boundaries$EW <- BBS_mig_strata  %>% factor(levels = c("West","East"))

# Plot BBS strata grouped into east and west
ggplot(BBS_strata_boundaries, aes(fill = EW))+
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(col = "black")+
  scale_fill_manual(values = strata_colours, name = "Stratum Grouping")+
  ggtitle("BBS analytical strata")

st_comp_regions <- get_composite_regions(strata_type = "bbs_usgs")
st_comp_regions$EW <- "Other"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "East")$ST_12)] <- "East"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "West")$ST_12)] <- "West"

EW_indices <- generate_indices(jags_mod = out_BBS,
                               jags_data = jags_data_BBS,
                               alt_region_names = st_comp_regions,
                               regions = "EW")

# *******************************************************************
# *******************************************************************
# MAIN RESULTS
# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Figure 1: Maps of migration monitoring stations
#------------------------------------------------

# Convert to AEA 
strata_sf <- strata_sf %>% st_transform(st_crs(BBS_strata_boundaries))

map_Spring <- ggplot(data = strata_sf) +
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7, col = "transparent")+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5))+
  geom_sf(data = station_data_summarized_sf_Spring ,size = 1, alpha = 1, shape = 19,col = "black")+
  geom_sf_label_repel(data = station_data_summarized_sf_Spring,
                      aes(label = station),
                      max.overlaps = 50,
                      force = 1.5,
                      force_pull = 0.1,
                      size = 1.5, 
                      alpha = 0.8,
                      min.segment.length = 0)+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = "none") +
  ggtitle("Pre-breeding ('spring') migration")

map_Fall <- ggplot(data = strata_sf) +
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7, col = "transparent")+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5))+
  geom_sf(data = station_data_summarized_sf_Fall ,size = 1, alpha = 1, shape = 19,col = "black")+
  geom_sf_label_repel(data = station_data_summarized_sf_Fall,
                      aes(label = station),
                      max.overlaps = 50,
                      force = 1.5,
                      force_pull = 0.1,
                      size = 1.5, 
                      alpha = 0.8,
                      min.segment.length = 0)+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = "none") +
  ggtitle("Post-breeding ('fall') migration")


Figure1 <- plot_grid(map_Spring,map_Fall,nrow=1,labels=c("(a)","(b)"))

png(file ="analysis/2_output/Figures_Main_Text/Figure1.png", units = "in", width = 10, height = 5, res = 600)
Figure1
dev.off()

# ------------------------------------------------
# CALCULATE CONTINENTAL TRENDS, USING EBIRD FOR RELATIVE ABUNDANCES WITHIN EACH STRATUM
# ------------------------------------------------

# Load high resolution relative abundance raster
rast <- rast(paste0(path,"/seasonal/bkpwar_abundance_seasonal_mean_hr_2021.tif"))
rast <- rast[["breeding"]]

# Extract relative abundance in each stratum
strata_sf <- st_transform(strata_sf, st_crs(rast))
rast <- crop(rast, strata_sf)
rast_df <- terra::extract(rast,strata_sf) %>%
  group_by(ID) %>%
  summarize(RelAbund = sum(breeding, na.rm = TRUE)) %>%
  mutate(Stratum = c("East","West")[ID])

# Assume that ebird map represents relative abundance in final year.  Rescale everything X to this final year
X_Spring <- out_Spring$sims.list$X
X_Fall <- out_Fall$sims.list$X

for (j in 1:jags_data_Spring$nstrata){
  for (t in 1:jags_data_Spring$nyear){
    
    X_Spring[,j,t] <- X_Spring[,j,t]/X_Spring[,j,jags_data_Spring$nyear] * rast_df$RelAbund[j]
    X_Fall[,j,t] <- X_Fall[,j,t]/X_Fall[,j,jags_data_Fall$nyear] * rast_df$RelAbund[j]
    
    
  }
}

# Estimated total abundance each year
XSum_Spring <- apply(X_Spring,c(1,3),sum)
XSum_Fall <- apply(X_Fall,c(1,3),sum)

(rast_df$RelAbund[1]/rast_df$RelAbund[2]) %>% round(2) # Eastern stratum was 1.3 times larger than western stratum in 2021

#------------------------------------------------
# Figure 2: Regional population trajectories relative to 2000
#------------------------------------------------

regional_change_since_2000 <- data.frame()

for (j in 1:jags_data_Spring$nstrata){
  for (t in 1:jags_data_Spring$nyear){
    
    # ------------
    # Spring estimates
    # ------------
    
    log_change <- log(X_Spring[,j,t]/X_Spring[,j,1])
    
    regional_change_since_2000 <- rbind(regional_change_since_2000,
                               data.frame(stratum_number = j,
                                          year_number = t,
                                          Source = "Pre-breeding migration",
                                          
                                          log_change_q025 = quantile(log_change,0.025),
                                          log_change_q50 = quantile(log_change,0.5),
                                          log_change_q975 = quantile(log_change,0.975),
                                          
                                          prob_large_decrease = mean(log_change<log(0.5)),
                                          prob_large_increase = mean(log_change>log(2)),
                                          
                                          prob_decline = mean(log_change<0)
                               ))
    
    # ------------
    # Fall estimates
    # ------------
    
    log_change <- log(X_Fall[,j,t]/X_Fall[,j,1])
    
    regional_change_since_2000 <- rbind(regional_change_since_2000,
                               data.frame(stratum_number = j,
                                          year_number = t,
                                          Source = "Post-breeding migration",
                                          
                                          log_change_q025 = quantile(log_change,0.025),
                                          log_change_q50 = quantile(log_change,0.5),
                                          log_change_q975 = quantile(log_change,0.975),
                                          
                                          prob_large_decrease = mean(log_change<log(0.5)),
                                          prob_large_increase = mean(log_change>log(2)),
                                          
                                          prob_decline = mean(log_change<0)
                               ))
    
    # ------------
    # BBS estimates
    # ------------
    
    if (j==1) X <- EW_indices$samples$EW_East
    if (j==2) X <- EW_indices$samples$EW_West
    
    log_change <- log(X[,t]/X[,1])
    
    regional_change_since_2000 <- rbind(regional_change_since_2000,
                               data.frame(stratum_number = j,
                                          year_number = t,
                                          Source = "Breeding Bird Survey",
                                          
                                          log_change_q025 = quantile(log_change,0.025),
                                          log_change_q50 = quantile(log_change,0.5),
                                          log_change_q975 = quantile(log_change,0.975),
                                          
                                          prob_large_decrease = mean(log_change<log(0.5)),
                                          prob_large_increase = mean(log_change>log(2)),
                                          
                                          prob_decline = mean(log_change<0)
                               ))
    
    
  }
}

regional_change_since_2000$Year = (2000:2018)[regional_change_since_2000$year_number]
regional_change_since_2000$Stratum = c("East","West")[regional_change_since_2000$stratum_number] %>% factor(levels = c("West","East"))
regional_change_since_2000$Source <- factor(regional_change_since_2000$Source, levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"))

# ------------
# Prepare y axis scale (conversion between log-scale change and percent change) Convert log-scale change to percent change using: 100*(exp(log_change)-1)
# ------------

y_axis_breaks <- c(0,100,300,900)
log_breaks <- log(y_axis_breaks/100+1)
log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
y_axis_breaks <- 100*(exp(log_breaks)-1)

y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
y_axis_labels <- paste0(y_axis_labels," %")

# ------------
# Text labels for each panel
# ------------

text_labels <- expand.grid(Year = 2000,
                           log_y = max(log_breaks),
                           Stratum = factor(c("West","East"),levels = c("West","East")),
                           Source = factor(c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"),levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey")))
text_labels$panel <- c("(a)","(b)","(c)","(d)","(e)","(f)")
Figure2 <- ggplot()+
  geom_ribbon(data = regional_change_since_2000,
              aes(x = Year, 
                  y = log_change_q50, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975, 
                  fill = Stratum, 
                  col = Stratum),alpha = 0.2, col = "transparent")+
  geom_line(data = regional_change_since_2000,
            aes(x = Year, 
                y = log_change_q50, 
                ymin = log_change_q025, 
                ymax = log_change_q975, 
                fill = Stratum, 
                col = Stratum),linewidth = 2)+
  geom_text(data = text_labels, aes(x = Year, y = log_y, label = panel),hjust=0, size = 6,vjust=1)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_grid(Source~Stratum)+
  ylab("Percent change relative to 2000")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")

Figure2

png(file ="analysis/2_output/Figures_Main_Text/Figure2.png", units = "in", width = 6, height = 8, res = 600)
Figure2
dev.off()

#------------------------------------------------
# Figure 3: Estimates of continental change relative to 2000
#------------------------------------------------

BBS_indices <- generate_indices(jags_mod = out_BBS,
                            jags_data = jags_data_BBS,
                            regions = "continental")

continental_change_since_2000 <- data.frame()

for (t in 1:jags_data_Spring$nyear){
  
  # ------------
  # Spring estimates
  # ------------
  
  log_change <- log(XSum_Spring[,t]/XSum_Spring[,1])
  
  continental_change_since_2000 <- rbind(continental_change_since_2000,
                             data.frame(year_number = t,
                                        Source = "Pre-breeding migration",
                                        
                                        log_change_q025 = quantile(log_change,0.025),
                                        log_change_q50 = quantile(log_change,0.5),
                                        log_change_q975 = quantile(log_change,0.975),
                                        prob_large_decrease = mean(log_change<log(0.5)),
                                        prob_large_increase = mean(log_change>log(2)),
                                        
                                        prob_decline = mean(log_change<0)
                             ))
  
  # ------------
  # Fall estimates
  # ------------
  
  log_change <- log(XSum_Fall[,t]/XSum_Fall[,1])
  
  continental_change_since_2000 <- rbind(continental_change_since_2000,
                             data.frame(year_number = t,
                                        Source = "Post-breeding migration",
                                        
                                        log_change_q025 = quantile(log_change,0.025),
                                        log_change_q50 = quantile(log_change,0.5),
                                        log_change_q975 = quantile(log_change,0.975),
                                        prob_large_decrease = mean(log_change<log(0.5)),
                                        prob_large_increase = mean(log_change>log(2)),
                                        
                                        prob_decline = mean(log_change<0)
                             ))
  
  # ------------
  # BBS estimates
  # ------------
  
  log_change <- log(BBS_indices$samples$continental_Continental[,t]/BBS_indices$samples$continental_Continental[,1])
  
  continental_change_since_2000 <- rbind(continental_change_since_2000,
                             data.frame(year_number = t,
                                        Source = "Breeding Bird Survey",
                                        
                                        log_change_q025 = quantile(log_change,0.025),
                                        log_change_q50 = quantile(log_change,0.5),
                                        log_change_q975 = quantile(log_change,0.975),
                                        prob_large_decrease = mean(log_change<log(0.5)),
                                        prob_large_increase = mean(log_change>log(2)),
                                        
                                        prob_decline = mean(log_change<0)
                             ))
  
  
}

continental_change_since_2000$Year = (2000:2018)[continental_change_since_2000$year_number]
continental_change_since_2000$Source <- factor(continental_change_since_2000$Source, levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"))

# ------------
# Prepare y axis scale (conversion between log-scale change and percent change) Convert log-scale change to percent change using: 100*(exp(log_change)-1)
# ------------

y_axis_breaks <- c(0,100,300,900)
log_breaks <- log(y_axis_breaks/100+1)
log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
y_axis_breaks <- 100*(exp(log_breaks)-1)

y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
y_axis_labels <- paste0(y_axis_labels," %")

Figure3 <- ggplot(continental_change_since_2000,
                  aes(x = Year, 
                      y = log_change_q50, 
                      ymin = log_change_q025, 
                      ymax = log_change_q975))+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 1)+
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_grid(.~Source)+
  ylab("Percent change relative to 2000")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")

Figure3

png(file ="analysis/2_output/Figures_Main_Text/Figure3.png", units = "in", width = 8, height = 4, res = 600)
Figure3
dev.off()

#-------------------------------------------------------------------
# Summarize change/trend estimates from each program, for each region
#-------------------------------------------------------------------
# 
# trend_fn <- function(N2,N1,year_interval){
#   percent_change = (N2-N1)/N1 * 100
#   trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
#   return(trend)
# }
# 
# trend_df <- data.frame()
# for (s in 1:jags_data$nstrata){
#   for (start_year in 1:jags_data$nyear){
#     for (end_year in 1:jags_data$nyear){
#       if (end_year <= start_year) next
# 
#       trend_est <- trend_fn(N2 = out$sims.list$X[,s,end_year], N1 = out$sims.list$X[,s,start_year], year_interval = end_year - start_year)
# 
#       trend_df <- rbind(trend_df, data.frame(stratum_number = s,
#                                              start_year = start_year,
#                                              end_year = end_year,
#                                              trend_est_q50 = quantile(trend_est,0.5),
#                                              trend_est_q025 = quantile(trend_est,0.025),
#                                              trend_est_q975 = quantile(trend_est,0.975),
#                                              prob_decline = round(mean(trend_est<0),2)))
#     }
#   }
# 
# }
# trend_df <- trend_df %>% arrange(start_year)
# trend_df$start_year_abs <- year_vec[trend_df$start_year]
# trend_df$end_year_abs <- year_vec[trend_df$end_year]
# trend_df$Stratum <- c("East","West")[trend_df$stratum_number]
# 
# 
# trend_df$'Trend Period' <- paste0(trend_df$start_year_abs," to ", trend_df$end_year_abs)
# trend_df$'Trend Length (Yrs)' <- trend_df$end_year_abs - trend_df$start_year_abs
# trend_df$'Trend Estimate' <- paste0(round(trend_df$trend_est_q50,1)," (",round(trend_df$trend_est_q025,1)," to ",round(trend_df$trend_est_q975,1),")")
# trend_df$'P(Decline)' <- trend_df$prob_decline
# trend_df$'Width of 95% CRI' <- round(trend_df$trend_est_q975 - trend_df$trend_est_q025,2)
# 
# trend_summary <- trend_df %>%
#   subset(end_year_abs == max(trend_df$end_year_abs) & trend_df$'Trend Length' %in% c(10,max(trend_df$'Trend Length'))) %>%
#   dplyr::select(Stratum,'Trend Period','Trend Length (Yrs)',"Trend Estimate",'Width of 95% CRI','P(Decline)') %>%
#   arrange('Trend Length (Yrs)')
# trend_summary
# 
# write.csv(trend_summary, file = paste0(output_directory,"/tables/",focal_season,"_trend_summary.csv"),row.names = FALSE)


#------------------------------------------------
# Figure 4: Station-level annual indices (overlaid with estimated proportion of birds from each stratum)
#------------------------------------------------

#-------
# Spring
#-------

stations_Spring <- count_df_Spring[,c("station","station_number")] %>% unique() %>% arrange(station_number)

# Station-level indices
station_indices_Spring <- out_Spring$sims.list$T_star %>% 
  reshape2::melt() %>%
  rename(sample = Var1, station_number = Var2, year_number = Var3, index = value) %>%
  mutate(Station = stations_Spring$station[station_number],
         Year = (2000:2018)[year_number]) %>%
  group_by(Station,Year) %>%
  summarize(index_q50 = quantile(index,0.5),
            index_q025 = quantile(index,0.025),
            index_q975 = quantile(index,0.975))

# Annual composition (mean)
station_composition_Spring <- out_Spring$mean$station_composition %>%
  reshape2::melt() %>%
  rename(stratum_number = Var1, station_number = Var2, year_number = Var3, composition = value) %>%
  mutate(Stratum = factor(c("East","West")[stratum_number],levels = c("West","East")),
         Station = stations_Spring$station[station_number],
         Year = (2000:2018)[year_number]) 

# Only include station-year combinations where count data were collected
counts_collected_Spring <- count_df_Spring %>%
  mutate(Station = stations_Spring$station[station_number],
         Year = (2000:2018)[year_number]) %>%
  group_by(Station,Year) %>%
  summarize(counts_collected = sum(count)>0)

station_indices_Spring <- station_indices_Spring %>% full_join(counts_collected_Spring) %>% na.omit()
station_composition_Spring <- station_composition_Spring %>% full_join(counts_collected_Spring) %>% na.omit()

# Plot
plot_station_composition_Spring <- ggplot() +
  geom_bar(data = station_composition_Spring,
           aes(x = Year, 
               y = composition,  
               fill = Stratum),
           alpha = 0.7,
           stat = "identity")+
  geom_point(data = station_indices_Spring,
             aes(x = Year, 
                 y = index_q50),
             size = 0.5)+
  geom_errorbar(data = station_indices_Spring,
                aes(x = Year, 
                    ymin = index_q025,
                    ymax = index_q975),
                width = 0,
                alpha = 0.5)+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~., scales = "free_y", ncol = 3)+
  ggtitle("Esimates of annual station composition\n\nPre-breeding ('Spring') migration")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Index of abundance")

png(file ="analysis/2_output/Figures_Main_Text/Figure4_Spring.png", units = "in", width = 8, height = 8, res = 600)
plot_station_composition_Spring
dev.off()

# Mean annual proportion of birds from each stratum, at each station
station_composition_sum_Spring <- station_composition_Spring %>% 
  group_by(Station,Year) %>%
  summarize(Sum = sum(composition)) %>%
  full_join(station_composition_Spring) %>%
  mutate(proportion = composition/Sum) %>%
  group_by(Station,Stratum) %>%
  summarize(Mean_Proportion = mean(proportion))

#-------
# Fall
#-------

stations_Fall <- count_df_Fall[,c("station","station_number")] %>% unique() %>% arrange(station_number)

# Station-level indices
station_indices_Fall <- out_Fall$sims.list$T_star %>% 
  reshape2::melt() %>%
  rename(sample = Var1, station_number = Var2, year_number = Var3, index = value) %>%
  mutate(Station = stations_Fall$station[station_number],
         Year = (2000:2018)[year_number]) %>%
  group_by(Station,Year) %>%
  summarize(index_q50 = quantile(index,0.5),
            index_q025 = quantile(index,0.025),
            index_q975 = quantile(index,0.975))

# Annual composition (mean)
station_composition_Fall <- out_Fall$mean$station_composition %>%
  reshape2::melt() %>%
  rename(stratum_number = Var1, station_number = Var2, year_number = Var3, composition = value) %>%
  mutate(Stratum = factor(c("East","West")[stratum_number],levels = c("West","East")),
         Station = stations_Fall$station[station_number],
         Year = (2000:2018)[year_number]) 

# Only include station-year combinations where count data were collected
counts_collected_Fall <- count_df_Fall %>%
  mutate(Station = stations_Fall$station[station_number],
         Year = (2000:2018)[year_number]) %>%
  group_by(Station,Year) %>%
  summarize(counts_collected = sum(count)>0)

station_indices_Fall <- station_indices_Fall %>% full_join(counts_collected_Fall) %>% na.omit()
station_composition_Fall <- station_composition_Fall %>% full_join(counts_collected_Fall) %>% na.omit()

# Plot
plot_station_composition_Fall <- ggplot() +
  geom_bar(data = station_composition_Fall,
           aes(x = Year, 
               y = composition,  
               fill = Stratum),
           alpha = 0.7,
           stat = "identity")+
  geom_point(data = station_indices_Fall,
             aes(x = Year, 
                 y = index_q50),
             size = 0.5)+
  geom_errorbar(data = station_indices_Fall,
                aes(x = Year, 
                    ymin = index_q025,
                    ymax = index_q975),
                width = 0,
                alpha = 0.5)+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~., scales = "free_y", ncol = 3)+
  ggtitle("Esimates of annual station composition\n\nPost-breeding ('fall') migration")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Index of abundance")

png(file ="analysis/2_output/Figures_Main_Text/Figure4_Fall.png", units = "in", width = 8, height = 8, res = 600)
plot_station_composition_Fall
dev.off()



# # *******************************************************************
# # *******************************************************************
# # BONUS FIGURES: CHANGE RELATIVE TO 2008
# # *******************************************************************
# # *******************************************************************
# year_vec <- 2000:2018
# start_year <- which(year_vec == 2008)
# 
# change_since_2008 <- data.frame()
# 
# for (t in start_year:jags_data_Spring$nyear){
#   
#   # ------------
#   # Spring estimates
#   # ------------
#   
#   log_change <- log(XSum_Spring[,t]/XSum_Spring[,start_year])
#   
#   change_since_2008 <- rbind(change_since_2008,
#                              data.frame(year_number = t,
#                                         Source = "Pre-breeding migration",
#                                         
#                                         log_change_q025 = quantile(log_change,0.025),
#                                         log_change_q50 = quantile(log_change,0.5),
#                                         log_change_q975 = quantile(log_change,0.975),
#                                         
#                                         prob_decline = mean(log_change<0)
#                              ))
#   
#   # ------------
#   # Fall estimates
#   # ------------
#   
#   log_change <- log(XSum_Fall[,t]/XSum_Fall[,start_year])
#   
#   change_since_2008 <- rbind(change_since_2008,
#                              data.frame(year_number = t,
#                                         Source = "Post-breeding migration",
#                                         
#                                         log_change_q025 = quantile(log_change,0.025),
#                                         log_change_q50 = quantile(log_change,0.5),
#                                         log_change_q975 = quantile(log_change,0.975),
#                                         
#                                         prob_decline = mean(log_change<0)
#                              ))
#   
#   # ------------
#   # BBS estimates
#   # ------------
#   
#   log_change <- log(BBS_indices$samples$continental_Continental[,t]/BBS_indices$samples$continental_Continental[,start_year])
#   
#   change_since_2008 <- rbind(change_since_2008,
#                              data.frame(year_number = t,
#                                         Source = "Breeding Bird Survey",
#                                         
#                                         log_change_q025 = quantile(log_change,0.025),
#                                         log_change_q50 = quantile(log_change,0.5),
#                                         log_change_q975 = quantile(log_change,0.975),
#                                         
#                                         prob_decline = mean(log_change<0)
#                              ))
#   
#   
# }
# 
# 
# change_since_2008$Year = year_vec[change_since_2008$year_number]
# change_since_2008$Source <- factor(change_since_2008$Source, levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"))
# 
# # ------------
# # Prepare y axis scale (conversion between log-scale change and percent change) Convert log-scale change to percent change using: 100*(exp(log_change)-1)
# # ------------
# 
# y_axis_breaks <- c(0,100,300,900)
# log_breaks <- log(y_axis_breaks/100+1)
# log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
# y_axis_breaks <- 100*(exp(log_breaks)-1)
# 
# y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
# y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
# y_axis_labels <- paste0(y_axis_labels," %")
# 
# Figure_SXX <- ggplot(change_since_2008,
#                   aes(x = Year, 
#                       y = log_change_q50, 
#                       ymin = log_change_q025, 
#                       ymax = log_change_q975))+
#   geom_ribbon(alpha = 0.2, col = "transparent")+
#   geom_line(linewidth = 1)+
#   scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
#   coord_cartesian(ylim=range(log_breaks))+
#   facet_grid(.~Source)+
#   ylab("Percent change relative to 2008")+
#   xlab("Year")+
#   geom_hline(yintercept = 0, linetype = 2)+
#   ggtitle("")
# 
# Figure_SXX
# 
# png(file ="analysis/2_output/Figures_Main_Text/Figure_SXX.png", units = "in", width = 8, height = 4, res = 600)
# Figure_SXX
# dev.off()
