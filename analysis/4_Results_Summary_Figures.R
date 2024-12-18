# *******************************************************************
# *******************************************************************
# Summary of Main Results for Publication
# Note that figures in appendices are produced in earlier model-specific scripts (e.g., goodness-of-fit evaluation)
# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse','sf','cowplot','bbsBayes','ebirdst','terra','scales','gridExtra') 
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

#dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
#BBS_strata_boundaries <- sf::read_sf(dsn = dsn)
BBS_strata_boundaries <- sf::read_sf("C:/Users/IlesD/AppData/Local/R/win-library/4.3/bbsBayes/maps/BBS_USGS_strata.shp")

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
  geom_sf(data = strata$BCR_poly, fill = "gray96", col = "gray85", linewidth = 0.75)+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7, col = "transparent")+
  xlab("")+
  ylab("")+
  geom_sf(data = all_stations,size = 2.5, alpha = 1, shape = 19, col = "black")+
  geom_sf(data = all_stations,aes(col = Season),size = 2, alpha = 1, shape = 19)+
  
  geom_sf_label_repel(data = all_stations,
                      aes(label = station, col = Season),
                      max.overlaps = 50,
                      force = 2,
                      force_pull = 0.1,
                      size = 4, 
                      alpha = 0.8,
                      min.segment.length = 0,
                      fontface = "bold",
                      show_guide = FALSE)+
  scale_color_manual(values=c(Season_Cols[1],"black",Season_Cols[2]))+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = "none") +
  ggtitle("")+
  coord_sf(xlim = c(-3500000,2800000),
           ylim=c(-500000,4200000))+
  
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5),
        axis.text=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position=c(.89,.87),
        legend.background = element_rect(fill = NA))

map_combined

png(file ="1_output/Results_MainText/Figure1.png", units = "in", width = 7.5, height = 7, res = 600)
map_combined
dev.off()

# *******************************************************************
# *******************************************************************
# Summarize and plot estimates of annual indices from each program
#  - rescale annual indices based on relative abundance from BAM and eBird
# *******************************************************************
# *******************************************************************

# ---------------------------------
# LOAD RELATIVE ABUNDANCE ESTIMATES
# ---------------------------------

Relative_Abundance_Sum <- read.csv("1_output/Results_MainText/relabund_eBird_BAM.csv")

# ************************************************
# Estimate national trend using both eBird and BAM relative abundance rasters
# ************************************************

Summary_Table <- data.frame()
analysis_summary <- data.frame()

for (source_of_estimate in c("BAM","eBird")){
  
  relabund <- subset(Relative_Abundance_Sum, Source == source_of_estimate)
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
  
  analysis_summary_tmp <- all_indices %>%
    group_by(Year,Stratum,Source) %>%
    summarize(
      Relabund_raster = source_of_estimate,
      index_q025 = quantile(indices_rescaled,0.025),
      index_q500 = quantile(indices_rescaled,0.500),
      index_q975 = quantile(indices_rescaled,0.975),
      
      
      log_change_q025 = quantile(log_change,0.025) ,
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
  
  analysis_summary <- rbind(analysis_summary, analysis_summary_tmp)
  
  # -----------------------
  # Output a table
  # -----------------------
  
  # Summary in final year (represents overall trend from 1998 to 2008)
  summary_tmp <- subset(analysis_summary_tmp, Year == 2018) %>%
    group_by(Stratum,Source, Relabund_raster = source_of_estimate) %>%
    summarize(`20 year trend` = paste0(trend_q500," (",trend_q025," to ",trend_q975,")"),
              `Prob trend is positive` = prob_positive,
              `% change since 1998` = paste0(percent_change_q500," (",percent_change_q025," to ",percent_change_q975,")"))
  
  
  Summary_Table <- rbind(Summary_Table, summary_tmp)
}

write.csv(Summary_Table,file = paste0("1_output/Results_MainText/Summary_Table_",source_of_estimate,".csv"), row.names = FALSE)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot percent change relative to baseline year within each stratum
#  (Figure 2 in manuscript)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analysis_summary_stratum <- subset(analysis_summary, Stratum %in% c("West","East"))

strata_colours <- c("#016b9b","#D18A80")

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
                           
                           Source = factor(c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey"),
                                           levels = c("Pre-breeding migration","Post-breeding migration","Breeding Bird Survey")),
                           Stratum = factor(c("West","East"),levels = c("West","East")))
text_labels$panel <- c("(a)","(b)","(c)","(d)","(e)","(f)")

# ------------
# Plot itself
# ------------


Fig2 <- ggplot()+
  
  geom_ribbon(data = analysis_summary_stratum,
              
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975, 
                  fill = Stratum, 
                  col = Stratum),alpha = 0.2, col = "transparent")+
  
  geom_line(data = analysis_summary_stratum,
            aes(x = Year, 
                y = log_change_q500, 
                col = Stratum),linewidth = 2)+
  geom_text(data = text_labels, aes(x = Year, y = log_y, label = panel),hjust=0, size = 6,vjust=1)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  
  facet_grid(Stratum~Source)+
  
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=11))
Fig2

png(file = "1_output/Results_MainText/Figure2.png", units = "in", width = 7.5, height = 6, res = 600)
Fig2
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot continental indices to baseline year within each stratum
#  (Figure 3 in manuscript)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analysis_summary_continental <- subset(analysis_summary, Stratum %in% c("Continental")) %>%
  subset(.,!(Source == "Breeding Bird Survey" & Relabund_raster == "eBird"))
analysis_summary_continental$label <- paste0(analysis_summary_continental$Source,"\n(",analysis_summary_continental$Relabund_raster,")")
analysis_summary_continental$label[analysis_summary_continental$Source == "Breeding Bird Survey"] <- "Breeding Bird Survey"

analysis_summary_continental$label <- factor(analysis_summary_continental$label,
                                             levels = c("Pre-breeding migration\n(BAM)",
                                                        "Post-breeding migration\n(BAM)",
                                                        "Breeding Bird Survey",
                                                        "Pre-breeding migration\n(eBird)",
                                                        "Post-breeding migration\n(eBird)"))

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

text_labels <- data.frame(Year = 1998,
                          log_y = max(log_breaks),
                          
                          label = c("Pre-breeding migration\n(BAM)",
                                    "Post-breeding migration\n(BAM)",
                                    "Breeding Bird Survey",
                                    "Pre-breeding migration\n(eBird)",
                                    "Post-breeding migration\n(eBird)"),
                          
                          panel = c("(a)","(b)","(e)","(c)","(d)"))

text_labels$label <- factor(text_labels$label,
                            levels = c("Pre-breeding migration\n(BAM)",
                                       "Post-breeding migration\n(BAM)",
                                       "Breeding Bird Survey",
                                       "Pre-breeding migration\n(eBird)",
                                       "Post-breeding migration\n(eBird)"))

Fig3_pt1 <- ggplot()+
  
  geom_ribbon(data = subset(analysis_summary_continental, label != "Breeding Bird Survey"),
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = subset(analysis_summary_continental, label != "Breeding Bird Survey"),
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  
  geom_text(data = subset(text_labels, label != "Breeding Bird Survey"), aes(x = Year, y = log_y, label = panel), hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")+
  theme(plot.margin=unit(c(0,0.5,0,0.1), "lines"))
Fig3_pt1

Fig3_pt2 <- ggplot()+
  
  geom_ribbon(data = subset(analysis_summary_continental, label == "Breeding Bird Survey"),
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = subset(analysis_summary_continental, label == "Breeding Bird Survey"),
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  
  geom_text(data = subset(text_labels, label == "Breeding Bird Survey"), aes(x = Year, y = log_y, label = panel), hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab("")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")+
  theme(plot.margin=unit(c(0,0.5,0,-1.5), "lines"))
Fig3_pt2

png(file = "1_output/Results_MainText/Figure3_revised.png", units = "in", width = 7.5, height = 6, res = 600)
grid.arrange(Fig3_pt1, Fig3_pt2,ncol = 3,nrow=4, 
             layout_matrix = cbind(c(1,1,1,1), 
                                   c(1,1,1,1),
                                   c(NA,2,2,NA)),
             
             heights = c(1,1.1,1.1,1),
             widths = c(1,1,1.035))
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot continental indices to baseline year within each stratum
#  (Figure 3 in manuscript)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(gridExtra)

analysis_summary_continental <- subset(analysis_summary, Stratum %in% c("Continental")) %>%
  subset(.,!(Source == "Breeding Bird Survey" & Relabund_raster == "eBird"))
analysis_summary_continental$label <- paste0(analysis_summary_continental$Source,"\n (",analysis_summary_continental$Relabund_raster," weighting)")

analysis_summary_continental$label[analysis_summary_continental$Source == "Breeding Bird Survey"] <- "Breeding Bird Survey"

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
# Plot itself
# ------------

dat_a <- subset(analysis_summary_continental, Source == "Pre-breeding migration" & Relabund_raster == "BAM")
Fig3a <- ggplot()+
  
  geom_ribbon(data = dat_a,
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = dat_a,
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  geom_text(aes(x = 1998, y = max(log_breaks), label = "(a)"),hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")
Fig3a

dat_b <- subset(analysis_summary_continental, Source == "Pre-breeding migration" & Relabund_raster == "eBird")
Fig3b <- ggplot()+
  
  geom_ribbon(data = dat_b,
              
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),
              alpha = 0.2, col = "transparent")+
  
  geom_line(data = dat_b,
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  geom_text(aes(x = 1998, y = max(log_breaks), label = "(b)"),hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab(" ")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")
Fig3b

dat_c <- subset(analysis_summary_continental, Source == "Post-breeding migration" & Relabund_raster == "BAM")
Fig3c <- ggplot()+
  
  geom_ribbon(data = dat_c,
              
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = dat_c,
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  geom_text(aes(x = 1998, y = max(log_breaks), label = "(c)"),hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab("Percent change relative to 1998")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")
Fig3c

dat_d <- subset(analysis_summary_continental, Source == "Post-breeding migration" & Relabund_raster == "eBird")
Fig3d <- ggplot()+
  
  geom_ribbon(data = dat_d,
              
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = dat_d,
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  geom_text(aes(x = 1998, y = max(log_breaks), label = "(d)"),hjust=0, size = 6,vjust=1)+
  
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab(" ")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")
Fig3d

dat_e <- subset(analysis_summary_continental, Source == "Breeding Bird Survey")
Fig3e <- ggplot()+
  
  geom_ribbon(data = dat_e,
              
              aes(x = Year, 
                  y = log_change_q500, 
                  ymin = log_change_q025, 
                  ymax = log_change_q975),alpha = 0.2, col = "transparent")+
  
  geom_line(data = dat_e,
            aes(x = Year, 
                y = log_change_q500),linewidth = 2)+
  geom_text(aes(x = 1998, y = max(log_breaks), label = "(e)"),hjust=0, size = 6,vjust=1)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  coord_cartesian(ylim=range(log_breaks))+
  facet_wrap(label~.)+
  ylab(" ")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("")
Fig3e


png(file = "1_output/Results_MainText/Figure3.png", units = "in", 
    width = 10, 
    height = 6, res = 600)
grid.arrange(Fig3a, Fig3b, Fig3c, Fig3d, Fig3e,ncol = 3,nrow=4, 
             layout_matrix = cbind(c(1,1,3,3), 
                                   c(2,2,4,4),
                                   c(NA,5,5,NA)))
dev.off()