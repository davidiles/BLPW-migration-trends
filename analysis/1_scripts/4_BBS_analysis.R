# *******************************************************************
# *******************************************************************
# Use bbsBayes to estimate continental trends of Blackpoll Warbler, 
# and compare to those from migration monitoring
# *******************************************************************
# *******************************************************************

library(bbsBayes)
library(tidyverse)
library(sf)
library(sp)

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

#------------------------------------------------
# Graphical themes
#------------------------------------------------

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

# --------------------------------------------
# Load migration analysis shapefile
# --------------------------------------------
# Load fitted model
load(paste0("analysis/2_output/Spring/analysis_Spring.RData"))


# ------------------------------------------------
# Set working directory
# ------------------------------------------------

strat_data <- stratify(by = "bbs_usgs")

species = "Blackpoll Warbler"
species_file = gsub(species,pattern = " ",replacement = "_")

jags_data <- prepare_jags_data(strat_data = strat_data,
                               species_to_run = species,
                               model = "firstdiff",
                               heavy_tailed = FALSE,
                               min_year = 2000,
                               max_year = 2018)

nsamp <- 2000
nb <- 50000
nt <- 100
ni <- nb + nsamp*nt

jags_mod <- run_model(jags_data = jags_data,
                      parameters_to_save = c("n","n3"),
                      n_iter = ni,
                      n_burnin = nb,
                      n_thin = nt,
                      parallel = TRUE)

max(unlist(jags_mod$Rhat),na.rm = TRUE)

save.image("analysis/2_output/BBS/BBS_analysis.RData")

# Load model output
load(file = "analysis/2_output/BBS/BBS_analysis.RData")

# --------------------------------------------
# Estimate continental trends
# --------------------------------------------

indices <- generate_indices(jags_mod = jags_mod,
                            jags_data = jags_data,
                            regions = "continental")

trends <- generate_trends(indices = indices,
                          Min_year = 2000,
                          Max_year = 2018)


tp = plot_indices(indices = indices,
                  species = "Blackpoll Warbler",
                  add_observed_means = TRUE)

png(file = "analysis/2_output/BBS/figures/BBS_continental.png", units = "in", width = 8, height = 6, res = 600)
tp
dev.off()

# --------------------------------------------
# Routes included in analysis
# --------------------------------------------

route_locations <- subset(strat_data$route_strat, 
                          rt.uni %in% jags_data$route & Year >= 2000) %>%
  dplyr::select(Latitude,Longitude) %>%
  unique() %>%
  st_as_sf(coords = c("Longitude","Latitude"),
           crs = CRS("+init=epsg:4326"))

# --------------------------------------------
# Load BBS strata
# --------------------------------------------

dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

stratum_indices <- generate_indices(jags_mod = jags_mod,
                                    jags_data = jags_data,
                                    regions = c("stratum"))

stratum_trends <- generate_trends(indices = stratum_indices,
                                  Min_year = 2008,
                                  Max_year = 2018)

# Only include strata that were included in analysis
BBS_strata_boundaries <- subset(BBS_strata_boundaries, ST_12 %in% stratum_trends$Strata_included)

# --------------------------------------------
# Change projection for migration strata
# --------------------------------------------

strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))
BCR_poly <- strata$BCR_poly

strata_sf <- strata_sf %>% st_transform(st_crs(BBS_strata_boundaries))
BCR_poly <- BCR_poly %>% st_transform(st_crs(BBS_strata_boundaries))
xlim <- c(-162, -56)
ylim <- c(35, 70)

strata_fig <- ggplot() +
  geom_sf(data = BCR_poly, fill = "gray90", col = "transparent")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7, col = "transparent")+
  #xlab("Longitude")+
  #ylab("Latitude")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray85", linetype = "dashed", size = 0.5))+
  #theme(axis.title=element_blank(),
  #      axis.text=element_blank(),
  #      axis.ticks=element_blank())+
  geom_sf(data = BBS_strata_boundaries, fill = "black", col = "gray25",
          alpha = 0.1,
          linewidth = 0.3)+
  geom_sf(data = route_locations, size = 0.3)+
  scale_fill_manual(values = strata_colours, name = "Migration\nAnalysis\nStrata") +
  ggtitle("North American Breeding Bird Survey\nanalytical strata and route locations\ncontributing data to analysis")
strata_fig

png(file = "analysis/2_output/BBS/figures/BBS_data.png", units = "in", width = 5, height = 4, res = 600)
strata_fig 
dev.off()

# --------------------------------------------
# Estimate trends within custom strata
# --------------------------------------------

# Overlap between migration strata and BBS strata
BBS_mig_strata = st_intersects(BBS_strata_boundaries,strata_sf) %>% as.numeric()
BBS_mig_strata[is.na(BBS_mig_strata)] <- 1
BBS_mig_strata <- c("East","West")[BBS_mig_strata] %>% factor(levels = c("West","East"))
BBS_strata_boundaries$East_West <- BBS_mig_strata  %>% factor(levels = c("West","East"))

ggplot(BBS_strata_boundaries, aes(fill = East_West))+
  geom_sf(col = "black")+
  scale_fill_manual(values = strata_colours, name = "Stratum")+
  ggtitle("BBS Strata")

st_comp_regions <- get_composite_regions(strata_type = "bbs_usgs")
st_comp_regions$East_West <- "Other"
st_comp_regions$East_West[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, East_West == "East")$ST_12)] <- "East"
st_comp_regions$East_West[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, East_West == "West")$ST_12)] <- "West"


east_west_indices <- generate_indices(jags_mod = jags_mod,
                                      jags_data = jags_data,
                                      alt_region_names = st_comp_regions,
                                      regions = "East_West")

east_west_trends_2000_2018 <- generate_trends(indices = east_west_indices,
                                              Min_year = 2000,
                                              Max_year = 2018)

east_west_trends_2008_2018 <- generate_trends(indices = east_west_indices,
                                              Min_year = 2008,
                                              Max_year = 2018)


# Plot trajectories
tp = plot_indices(indices = east_west_indices,
                  species = "Blackpoll Warbler",
                  add_observed_means = TRUE,
                  min_year = 2000,
                  max_year = 2018)

library(cowplot)

east_west_trajectory_plot <- plot_grid(tp$West, tp$East, nrow=1, ncol = 2)

png(file = "analysis/2_output/BBS/figures/BBS_East_West_Trajectory.png", units = "in", width = 12, height = 6, res = 600)
east_west_trajectory_plot 
dev.off()

# **************************************************************************************
# **************************************************************************************
# Comparisons between trend estimates from spring migration, fall migration, and BBS
# **************************************************************************************
# **************************************************************************************

trend_fn <- function(N2,N1,year_interval){
  percent_change = (N2-N1)/N1 * 100
  trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
  return(trend)
}

# Load fitted model for spring
load("analysis/2_output/Spring/analysis_Spring.RData")
out_Spring <- out

# Load fitted model for fall
load("analysis/2_output/Fall/analysis_Fall_alternative.RData")
out_Fall <- out
rm(out)

# --------------------------------------------
# Trends from 2000 to 2018
# --------------------------------------------

# Trend samples (10-year)
yr_start <- which(year_vec == 2000)
yr_end <- which(year_vec == 2018)

BBS_East_10yr <- data.frame(Source = "BBS",
                            Stratum = "East",
                            Start_Year = yr_start,
                            End_Year = yr_end,
                            Trend_Length = yr_end - yr_start,
                            trend_samples = trend_fn(N2 = east_west_indices$samples$East_West_East[,yr_end],
                                                     N1 = east_west_indices$samples$East_West_East[,yr_start],
                                                     year_interval = yr_end-yr_start))

BBS_West_10yr <- data.frame(Source = "BBS",
                            Stratum = "West",
                            Start_Year = yr_start,
                            End_Year = yr_end,
                            Trend_Length = yr_end - yr_start,
                            trend_samples = trend_fn(N2 = east_west_indices$samples$East_West_West[,yr_end],
                                                     N1 = east_west_indices$samples$East_West_West[,yr_start],
                                                     year_interval = yr_end-yr_start))

Spring_East_10yr <- data.frame(Source = "Pre-breeding Migration",
                               Stratum = "East",
                               Start_Year = yr_start,
                               End_Year = yr_end,
                               Trend_Length = yr_end - yr_start,
                               trend_samples = trend_fn(N2 = out_Spring$sims.list$X[,1,yr_end],
                                                        N1 = out_Spring$sims.list$X[,1,yr_start],
                                                        year_interval = yr_end-yr_start))

Spring_West_10yr <-  data.frame(Source = "Pre-breeding Migration",
                                Stratum = "West",
                                Start_Year = yr_start,
                                End_Year = yr_end,
                                Trend_Length = yr_end - yr_start,
                                trend_samples = trend_fn(N2 = out_Spring$sims.list$X[,2,yr_end],
                                                         N1 = out_Spring$sims.list$X[,2,yr_start],
                                                         year_interval = yr_end-yr_start))

Fall_East_10yr <-  data.frame(Source = "Post-breeding Migration",
                              Stratum = "East",
                              Start_Year = yr_start,
                              End_Year = yr_end,
                              Trend_Length = yr_end - yr_start,
                              trend_samples = trend_fn(N2 = out_Fall$sims.list$X[,1,yr_end],
                                                       N1 = out_Fall$sims.list$X[,1,yr_start],
                                                       year_interval = yr_end-yr_start))

Fall_West_10yr <- data.frame(Source = "Post-breeding Migration",
                             Stratum = "West",
                             Start_Year = yr_start,
                             End_Year = yr_end,
                             Trend_Length = yr_end - yr_start,
                             trend_samples = trend_fn(N2 = out_Fall$sims.list$X[,2,yr_end],
                                                      N1 = out_Fall$sims.list$X[,2,yr_start],
                                                      year_interval = yr_end-yr_start))

trend_estimates <- rbind(BBS_East_10yr,
                         BBS_West_10yr,
                         Spring_East_10yr,
                         Spring_West_10yr,
                         Fall_East_10yr,
                         Fall_West_10yr )

trend_estimates$Stratum <- factor(trend_estimates$Stratum, levels = c("West","East"))

trend_estimates$Source <- factor(trend_estimates$Source, levels = c("BBS","Pre-breeding Migration","Post-breeding Migration"))
plot1 <- ggplot(data = trend_estimates, aes(x = Source, y = trend_samples, col = Stratum, fill = Stratum))+
  geom_hline(yintercept=0,linetype=2)+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975), alpha = 0.7,size = 1)+
  facet_grid(.~Stratum)+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = FALSE) +
  scale_color_manual(values = strata_colours, name = "Stratum", guide = FALSE) +
  theme_bw()+
  coord_cartesian(ylim=c(-20,20))+
  xlab("Source of Trend Estimate")+
  ylab("Trend Estimate\n\n(Posterior distribution)")+
  ggtitle(paste0("Trend estimates from ",year_vec[yr_start], " to ",year_vec[yr_end]))
plot1

# --------------------------------------------
# Trends from 2008 to 2018
# --------------------------------------------

# Trend samples (10-year)
yr_start <- which(year_vec == 2008)
yr_end <- which(year_vec == 2018)

BBS_East_10yr <- data.frame(Source = "BBS",
                            Stratum = "East",
                            Start_Year = yr_start,
                            End_Year = yr_end,
                            Trend_Length = yr_end - yr_start,
                            trend_samples = trend_fn(N2 = east_west_indices$samples$East_West_East[,yr_end],
                                                     N1 = east_west_indices$samples$East_West_East[,yr_start],
                                                     year_interval = yr_end-yr_start))

BBS_West_10yr <- data.frame(Source = "BBS",
                            Stratum = "West",
                            Start_Year = yr_start,
                            End_Year = yr_end,
                            Trend_Length = yr_end - yr_start,
                            trend_samples = trend_fn(N2 = east_west_indices$samples$East_West_West[,yr_end],
                                                     N1 = east_west_indices$samples$East_West_West[,yr_start],
                                                     year_interval = yr_end-yr_start))

Spring_East_10yr <- data.frame(Source = "Pre-breeding Migration",
                               Stratum = "East",
                               Start_Year = yr_start,
                               End_Year = yr_end,
                               Trend_Length = yr_end - yr_start,
                               trend_samples = trend_fn(N2 = out_Spring$sims.list$X[,1,yr_end],
                                                        N1 = out_Spring$sims.list$X[,1,yr_start],
                                                        year_interval = yr_end-yr_start))

Spring_West_10yr <-  data.frame(Source = "Pre-breeding Migration",
                                Stratum = "West",
                                Start_Year = yr_start,
                                End_Year = yr_end,
                                Trend_Length = yr_end - yr_start,
                                trend_samples = trend_fn(N2 = out_Spring$sims.list$X[,2,yr_end],
                                                         N1 = out_Spring$sims.list$X[,2,yr_start],
                                                         year_interval = yr_end-yr_start))

Fall_East_10yr <-  data.frame(Source = "Post-breeding Migration",
                              Stratum = "East",
                              Start_Year = yr_start,
                              End_Year = yr_end,
                              Trend_Length = yr_end - yr_start,
                              trend_samples = trend_fn(N2 = out_Fall$sims.list$X[,1,yr_end],
                                                       N1 = out_Fall$sims.list$X[,1,yr_start],
                                                       year_interval = yr_end-yr_start))

Fall_West_10yr <- data.frame(Source = "Post-breeding Migration",
                             Stratum = "West",
                             Start_Year = yr_start,
                             End_Year = yr_end,
                             Trend_Length = yr_end - yr_start,
                             trend_samples = trend_fn(N2 = out_Fall$sims.list$X[,2,yr_end],
                                                      N1 = out_Fall$sims.list$X[,2,yr_start],
                                                      year_interval = yr_end-yr_start))

trend_estimates <- rbind(BBS_East_10yr,
                         BBS_West_10yr,
                         Spring_East_10yr,
                         Spring_West_10yr,
                         Fall_East_10yr,
                         Fall_West_10yr )

trend_estimates$Stratum <- factor(trend_estimates$Stratum, levels = c("West","East"))

trend_estimates$Source <- factor(trend_estimates$Source, levels = c("BBS","Pre-breeding Migration","Post-breeding Migration"))
plot2 <- ggplot(data = trend_estimates, aes(x = Source, y = trend_samples, col = Stratum, fill = Stratum))+
  geom_hline(yintercept=0,linetype=2)+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975), alpha = 0.7,size = 1)+
  facet_grid(.~Stratum)+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = FALSE) +
  scale_color_manual(values = strata_colours, name = "Stratum", guide = FALSE) +
  theme_bw()+
  coord_cartesian(ylim=c(-20,20))+
  xlab("Source of Trend Estimate")+
  ylab("Trend Estimate\n\n(Posterior distribution)")+
  ggtitle(paste0("Trend estimates from ",year_vec[yr_start], " to ",year_vec[yr_end]))
plot2

# --------------------------------------------
# Calculate and plot change relative to 2000
# --------------------------------------------
BBS_change_2000 <- data.frame()

for (y in 1:jags_data$nyear){
  percent_change_East <- 100*(east_west_indices$samples$East_West_East[,y]-east_west_indices$samples$East_West_East[,1])/east_west_indices$samples$East_West_East[,1] 
  BBS_change_2000 <- rbind(BBS_change_2000,
                           data.frame(Year = year_vec[y],
                                      Stratum = "East",
                                      percent_change_q50 = quantile(percent_change_East,0.5),
                                      percent_change_q025 = quantile(percent_change_East,0.025),
                                      percent_change_q975 = quantile(percent_change_East,0.975)))
  
  percent_change_West <- 100*(east_west_indices$samples$East_West_West[,y]-east_west_indices$samples$East_West_West[,1])/east_west_indices$samples$East_West_West[,1] 
  BBS_change_2000 <- rbind(BBS_change_2000,
                           data.frame(Year = year_vec[y],
                                      Stratum = "West",
                                      percent_change_q50 = quantile(percent_change_West,0.5),
                                      percent_change_q025 = quantile(percent_change_West,0.025),
                                      percent_change_q975 = quantile(percent_change_West,0.975)))
  
  
}
BBS_change_2000$Stratum <- factor(BBS_change_2000$Stratum, levels = c("West","East"))

regional_trajectory_plot_2000 <- ggplot(data = BBS_change_2000, aes(x = Year, 
                                   y = percent_change_q50, 
                                   ymin = percent_change_q025,
                                   ymax = percent_change_q975,
                                   fill = Stratum,
                                   col = Stratum))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 2)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  facet_grid(Stratum~.)+
  ylab("Percent change relative to 2000")+
  xlab("Year")+
  ggtitle("Regional trajectory\n\nBBS")
regional_trajectory_plot_2000


BBS_change_2008 <- data.frame()

for (y in which(year_vec == 2008):jags_data$nyear){
  percent_change_East <- 100*(east_west_indices$samples$East_West_East[,y]-east_west_indices$samples$East_West_East[,which(year_vec == 2008)])/east_west_indices$samples$East_West_East[,which(year_vec == 2008)] 
  BBS_change_2008 <- rbind(BBS_change_2008,
                           data.frame(Year = year_vec[y],
                                      Stratum = "East",
                                      percent_change_q50 = quantile(percent_change_East,0.5),
                                      percent_change_q025 = quantile(percent_change_East,0.025),
                                      percent_change_q975 = quantile(percent_change_East,0.975)))
  
  percent_change_West <- 100*(east_west_indices$samples$East_West_West[,y]-east_west_indices$samples$East_West_West[,which(year_vec == 2008)])/east_west_indices$samples$East_West_West[,which(year_vec == 2008)] 
  BBS_change_2008 <- rbind(BBS_change_2008,
                           data.frame(Year = year_vec[y],
                                      Stratum = "West",
                                      percent_change_q50 = quantile(percent_change_West,0.5),
                                      percent_change_q025 = quantile(percent_change_West,0.025),
                                      percent_change_q975 = quantile(percent_change_West,0.975)))
  
  
}
BBS_change_2008$Stratum <- factor(BBS_change_2008$Stratum, levels = c("West","East"))

regional_trajectory_plot_2008 <- ggplot(data = BBS_change_2008, aes(x = Year, 
                                                                    y = percent_change_q50, 
                                                                    ymin = percent_change_q025,
                                                                    ymax = percent_change_q975,
                                                                    fill = Stratum,
                                                                    col = Stratum))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 2)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  facet_grid(Stratum~.)+
  ylab("Percent change relative to 2008")+
  xlab("Year")+
  ggtitle("Regional trajectory\n\nBBS")
regional_trajectory_plot_2008

# 
# bbs_dir <- rappdirs::app_dir("bbsBayes")$data()
# stratum_shapes <- sf::read_sf(paste0(rappdirs::app_dir("bbsBayes")$data(),"/maps/BBS_USGS_strata.shp"))
# stratum_trends_sf <- full_join(stratum_shapes,stratum_trends,by = c("ST_12" = "Region")) #%>% na.omit()
# 
# # Plot
# stratum_trends_sf$Trend_discrete <- cut(stratum_trends_sf$Trend/100,breaks = c(-Inf,-0.05,-0.025,0,0.025,0.05,Inf))
# stratum_trends_sf$Trend_discrete <- factor(stratum_trends_sf$Trend_discrete, levels = rev(levels(stratum_trends_sf$Trend_discrete)))
# colpal <- RColorBrewer::brewer.pal(length(levels(stratum_trends_sf$Trend_discrete)),"RdYlBu")
# 
# stratum_trends_sf <- stratum_trends_sf %>% st_transform(lcc)
# route_locations <- route_locations %>% st_transform(lcc)
# 
# coords <- st_coordinates(route_locations)
# xlim <- expand_range(range(coords[,1]), mul = 0.05)
# ylim <- expand_range(range(coords[,2]), mul = 0.05)
# BBS_trends_map <- ggplot(data = stratum_trends_sf) +
#   geom_sf(data = stratum_trends_sf, fill = "gray97", col = "gray90")+
#   #geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85")+
#   geom_sf(data = na.omit(stratum_trends_sf),aes(fill = Trend_discrete))+
#   geom_sf(data = route_locations, size = 0.25, shape = 19)+
#   coord_sf(xlim = xlim, ylim = ylim)+
#   xlab("Longitude")+
#   ylab("Latitude")+
#   scale_fill_manual(values = rev(colpal), 
#                     name = "10-year trend",
#                     labels = rev(c("< -0.05",
#                                    "-0.025 to -0.05",
#                                    "0 to -0.025",
#                                    "0 to 0.025",
#                                    "0.025 to 0.05",
#                                    ">0.05")), drop = FALSE
#   ) +
#   theme(panel.background = element_rect(fill = "transparent", colour = "black"),
#         panel.border = element_rect(fill = "transparent", color = "black"),
#         panel.grid.major = element_line(color = "gray90", linetype = "dashed", size = 0.5))+
#   ggtitle("BBS trend results")
# 
# BBS_trends_map
#
# # --------------------------------------------
# # Prepare table of output
# # --------------------------------------------
# 
# national_indices <- generate_indices(jags_mod = jags_mod,
#                                      jags_data = jags_data,
#                                      regions = c("national"))
# 
# national_trend_2000_2018 <- generate_trends(indices = national_indices,Min_year = 2000,Max_year = 2018)  %>%  subset(Region == "CA")
# national_trend_2008_2018 <- generate_trends(indices = national_indices,Min_year = 2008,Max_year = 2018)  %>%  subset(Region == "CA")
# 
# # Estimate of percent change between endpoints
# years <- seq(national_indices$startyear,national_indices$startyear+dim(national_indices$samples$national_CA)[2]-1)
# indices_2000 <- national_indices$samples$national_CA[,which(years == 2000)]
# indices_2008 <- national_indices$samples$national_CA[,which(years == 2008)]
# indices_2018 <- national_indices$samples$national_CA[,which(years == 2018)]
# 
# pchange_2000_2018 <- (indices_2018 - indices_2000)/indices_2000 * 100
# mean(pchange_2000_2018  <= -30)
# 
# pchange_2008_2018 <- (indices_2018 - indices_2008)/indices_2008 * 100
# mean(pchange_2008_2018  <= -30)
# 
# trend_table <- rbind(national_trend_2000_2018,
#                      national_trend_2008_2018)[,-6]
# 
# trend_table$prob_increase <- c(mean(pchange_2000_2018  > 0),mean(pchange_2008_2018  > 0))
# trend_table$prob_30_decline <- c(mean(pchange_2000_2018  <= -30),mean(pchange_2008_2018  <= -30))
# trend_table$prob_50_decline <- c(mean(pchange_2000_2018  <= -50),mean(pchange_2008_2018  <= -50))
# 
# trend_table <- trend_table %>% add_column(stratum_name = "National", season = "BBS") %>%
#   dplyr::select(stratum_name,season,Start_year, End_year,
#                 Trend,Trend_Q0.025,Trend_Q0.975,
#                 Percent_Change,Percent_Change_Q0.025,Percent_Change_Q0.975,
#                 prob_increase,
#                 prob_30_decline,
#                 prob_50_decline)
# 
# write.csv(trend_table, file = "BBS_trend_estimates.csv", row.names = FALSE)
