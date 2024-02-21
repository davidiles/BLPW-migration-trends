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

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BLPW-migration-trends/analysis")

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

# ------------------------------------------------
# Conduct analysis
# ------------------------------------------------

year_vec <- seq(1998,2018)

# strat_data <- stratify(by = "bbs_usgs")
# 
# species = "Blackpoll Warbler"
# species_file = gsub(species,pattern = " ",replacement = "_")
# 
# jags_data <- prepare_jags_data(strat_data = strat_data,
#                                species_to_run = species,
#                                model = "firstdiff",
#                                heavy_tailed = FALSE,
#                                min_year = min(year_vec),
#                                max_year = max(year_vec))
# 
# nsamp <- 2000
# nb <- 50000
# nt <- 100
# ni <- nb + nsamp*nt
# 
# jags_mod <- run_model(jags_data = jags_data,
#                       parameters_to_save = c("n","n3"),
#                       n_iter = ni,
#                       n_burnin = nb,
#                       n_thin = nt,
#                       parallel = TRUE)
# 
# max(unlist(jags_mod$Rhat),na.rm = TRUE)
# 
# save.image("1_output/BBS/BBS_analysis.RData")

# ------------------------------------------------
# Load fitted model
# ------------------------------------------------

load("1_output/BBS/BBS_analysis.RData")

# ------------------------------------------------
# Load migration strata
# ------------------------------------------------

load("0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
strata_sf$stratum_number <- 1:nrow(strata_sf)
strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))

#------------------------------------------------
# Load BBS strata shapefile
#------------------------------------------------

dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

# ------------------------------------------------
# Define custom 'east' and 'west' strata for BBS analysis; calculate indices in these strata
# ------------------------------------------------

stratum_indices <- generate_indices(jags_mod = jags_mod,
                                    jags_data = jags_data,
                                    regions = c("stratum"))

stratum_trends <- generate_trends(indices = stratum_indices,
                                  Min_year = min(year_vec),
                                  Max_year = max(year_vec))

# Only include strata that were included in analysis
BBS_strata_boundaries <- subset(BBS_strata_boundaries, ST_12 %in% stratum_trends$Strata_included)

# Overlap between migration strata and BBS strata
strata_sf <- st_transform(strata_sf, st_crs(BBS_strata_boundaries))
BBS_mig_strata = st_intersects(BBS_strata_boundaries,strata_sf) %>% as.numeric()
BBS_mig_strata[is.na(BBS_mig_strata)] <- 1
BBS_mig_strata <- c("East","West")[BBS_mig_strata] %>% factor(levels = c("West","East"))
BBS_strata_boundaries$EW <- BBS_mig_strata  %>% factor(levels = c("West","East"))

st_comp_regions <- get_composite_regions(strata_type = "bbs_usgs")
st_comp_regions$EW <- "Other"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "East")$ST_12)] <- "East"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "West")$ST_12)] <- "West"

BBS_EW <- generate_indices(jags_mod = jags_mod,
                           jags_data = jags_data,
                           alt_region_names = st_comp_regions,
                           regions = "EW")

BBS_East <- BBS_EW$samples$EW_East %>% 
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "East",
         Year = year_vec[year_number])

BBS_West <- BBS_EW$samples$EW_West %>% 
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "West",
         Year = year_vec[year_number])

# ------------------------------------------------
# Calculate continental indices
# ------------------------------------------------

BBS_Continental <- generate_indices(jags_mod = jags_mod,
                                    jags_data = jags_data,
                                    regions = "continental")$samples[[1]] %>%
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "Continental",
         Year = year_vec[year_number])

# ------------------------------------------------
# Merge into single table
# ------------------------------------------------

BBS_indices_df <- bind_rows(BBS_East,BBS_West,BBS_Continental)

BBS_indices_summarized <- BBS_indices_df %>%
  group_by(Stratum,Year) %>%
  summarize(
    
    # On original scale
    indices_q025 = quantile(indices,0.025),
    indices_q500 = quantile(indices,0.500),
    indices_q975 = quantile(indices,0.975))

# ------------------------------------------------
# Calculate trends
# ------------------------------------------------

yend = 2018 # End year
y0 = 1998   # Start year

BBS_trend_df <- BBS_indices_df %>%
  group_by(Stratum,sample) %>%
  mutate(trend = 100*((indices[Year == yend]/indices[Year == y0])^(1/(yend-y0))-1),
         perc_change_1998 = 100*(indices[Year == 2018] - indices[Year == 1998])/indices[Year == 1998],
         perc_change_2008 = 100*(indices[Year == 2018] - indices[Year == 2008])/indices[Year == 2008]
  )

BBS_trend_summarized <- BBS_trend_df %>%
  group_by(Stratum) %>%
  summarize(trend_q500 = quantile(trend,0.500) %>% round(1),
            trend_q025 = quantile(trend,0.025)  %>% round(1),
            trend_q975 = quantile(trend,0.975)  %>% round(1),
            
            prob_positive = mean(trend>0)  %>% round(2),
            
            perc_change_1998_q500 = quantile(perc_change_1998,0.500) %>% round(0),
            perc_change_1998_q025 = quantile(perc_change_1998,0.025)  %>% round(0),
            perc_change_1998_q975 = quantile(perc_change_1998,0.975)  %>% round(0),
            
            perc_change_2008_q500 = quantile(perc_change_2008,0.500) %>% round(0),
            perc_change_2008_q025 = quantile(perc_change_2008,0.025)  %>% round(0),
            perc_change_2008_q975 = quantile(perc_change_2008,0.975)  %>% round(0),
  )

BBS_trend_summarized 

BBS_results <- list(BBS_indices_df = BBS_indices_df,
                    BBS_indices_summarized = BBS_indices_summarized,
                    BBS_trend_df = BBS_trend_df,
                    BBS_trend_summarized = BBS_trend_summarized)

BBS_results$trend_summarized %>% as.data.frame()

saveRDS(BBS_results,"1_output/BBS/results_BBS.rds")




# ############################################
# Additional Plots
# ############################################



# --------------------------------------------
# Estimate continental trends
# --------------------------------------------

indices <- generate_indices(jags_mod = jags_mod,
                            jags_data = jags_data,
                            regions = "continental")

trends <- generate_trends(indices = indices,
                          Min_year = min(year_vec),
                          Max_year = max(year_vec))


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

year_vec <- seq(1998,2018)
# Define custom 'east' and 'west' strata for BBS analysis; calculate indices in these strata
stratum_indices <- generate_indices(jags_mod = out_BBS,
                                    jags_data = jags_data_BBS,
                                    regions = c("stratum"))

stratum_trends <- generate_trends(indices = stratum_indices,
                                  Min_year = 1998,
                                  Max_year = 2018)

# Only include strata that were included in analysis
BBS_strata_boundaries <- subset(BBS_strata_boundaries, ST_12 %in% stratum_trends$Strata_included)

# Overlap between migration strata and BBS strata
strata_sf <- st_transform(strata_sf, st_crs(BBS_strata_boundaries))
BBS_mig_strata = st_intersects(BBS_strata_boundaries,strata_sf) %>% as.numeric()
BBS_mig_strata[is.na(BBS_mig_strata)] <- 1
BBS_mig_strata <- c("East","West")[BBS_mig_strata] %>% factor(levels = c("West","East"))
BBS_strata_boundaries$EW <- BBS_mig_strata  %>% factor(levels = c("West","East"))

st_comp_regions <- get_composite_regions(strata_type = "bbs_usgs")
st_comp_regions$EW <- "Other"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "East")$ST_12)] <- "East"
st_comp_regions$EW[which(st_comp_regions$region %in% subset(BBS_strata_boundaries, EW == "West")$ST_12)] <- "West"

BBS_EW <- generate_indices(jags_mod = out_BBS,
                           jags_data = jags_data_BBS,
                           alt_region_names = st_comp_regions,
                           regions = "EW")

BBS_East <- BBS_EW$samples$EW_East %>% 
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "East",
         Year = year_vec[year_number])

BBS_West <- BBS_EW$samples$EW_West %>% 
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "West",
         Year = year_vec[year_number])

BBS_Continental <- generate_indices(jags_mod = out_BBS,
                                    jags_data = jags_data_BBS,
                                    regions = "continental")$samples[[1]] %>%
  reshape2::melt() %>% 
  dplyr::rename(sample = Var1, year_number = Var2, indices = value) %>%
  mutate(Stratum = "Continental",
         Year = year_vec[year_number])

BBS_indices_df <- bind_rows(BBS_East,BBS_West,BBS_Continental)

BBS_indices_summarized <- BBS_indices_df %>%
  group_by(Stratum,Year) %>%
  summarize(
    
    # On original scale
    indices_q025 = quantile(indices,0.025),
    indices_q500 = quantile(indices,0.500),
    indices_q975 = quantile(indices,0.975))

# --------------------------------------------
# Calculate and plot change relative to 1998
# --------------------------------------------

BBS_change_1998 <- data.frame()

for (y in 1:jags_data$nyear){
  percent_change_East <- 100*(east_west_indices$samples$East_West_East[,y]-east_west_indices$samples$East_West_East[,1])/east_west_indices$samples$East_West_East[,1] 
  BBS_change_1998 <- rbind(BBS_change_1998,
                           data.frame(Year = year_vec[y],
                                      Stratum = "East",
                                      percent_change_q50 = quantile(percent_change_East,0.5),
                                      percent_change_q025 = quantile(percent_change_East,0.025),
                                      percent_change_q975 = quantile(percent_change_East,0.975)))
  
  percent_change_West <- 100*(east_west_indices$samples$East_West_West[,y]-east_west_indices$samples$East_West_West[,1])/east_west_indices$samples$East_West_West[,1] 
  BBS_change_1998 <- rbind(BBS_change_1998,
                           data.frame(Year = year_vec[y],
                                      Stratum = "West",
                                      percent_change_q50 = quantile(percent_change_West,0.5),
                                      percent_change_q025 = quantile(percent_change_West,0.025),
                                      percent_change_q975 = quantile(percent_change_West,0.975)))
  
  
}
BBS_change_1998$Stratum <- factor(BBS_change_1998$Stratum, levels = c("West","East"))

regional_trajectory_plot_1998 <- ggplot(data = BBS_change_1998, aes(x = Year, 
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
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  ggtitle("Regional trajectory\n\nBBS")
regional_trajectory_plot_1998


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