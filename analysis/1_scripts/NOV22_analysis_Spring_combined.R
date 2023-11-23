# *******************************************************************
# *******************************************************************
# Analysis of pre-breeding (i.e., spring) migration data
# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse','viridis','RColorBrewer','cowplot','readxl','ggrepel',
              'rgeos','raster','sp','sf','rgbif','rmapshaper','spatialEco','gstat',
              'foreign',"rnaturalearth","rnaturalearthdata",'jagsUI','reshape2','parallel',"MCMCvis","scales") 
if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)
library(ggsflabel) # custom library download from devtools::install_github("yutannihilation/ggsflabel")

#------------------------------------------------
# Clear working memory
#------------------------------------------------

rm(list=ls())

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

# ***************************************************************
# ***************************************************************
# PART 1: PREPARE DATA AND FIT MODEL
# ***************************************************************
# ***************************************************************

#------------------------------------------------
# Load strata shapefile; initially created using combine_strata.R
#------------------------------------------------

load("analysis/0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
nstrata <- nrow(strata_sf)
strata_sf$stratum_number <- 1:nstrata

#------------------------------------------------
# Load / process migration and isotope data
#------------------------------------------------

focal_season <- "Fall"
start_year <- 1998
end_year <- 2018

# Relevant directories / set working directory
data_directory <- paste0("analysis/0_data/")
output_directory <- paste0("analysis/2_output/",focal_season,"/")
figure_directory <- paste0(output_directory,"/figures/")
table_directory <- paste0(output_directory,"/tables/")

# PART 1: FORMAT COUNTS
{
  
  # *****
  # Read/format data from BSC (Canadian data)
  # Data provided by Danielle Ethier (BSC).  Has been pre-processed/cleaned.  Does not include offsets or effort (i.e., hours nets were open)
  # *****
  dat_can <- rbind(read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/ACBO.BLPW.2018.csv") %>% add_column(., station = "ACBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/BPBO.BLPW.2018.csv") %>% add_column(., station = "BPBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/IPBO.BLPW.2018.csv") %>% add_column(., station = "IPBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/LMBO.BLPW.2018.csv") %>% add_column(., station = "LMBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/LPBO.BLPW.2018.csv") %>% add_column(., station = "LPBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/MGBO.BLPW.2018.csv") %>% add_column(., station = "MGBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/MNO.BLPW.2018.csv") %>% add_column(., station = "MNO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/PEPBO.BLPW.2018.csv") %>% add_column(., station = "PEPBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/PIBO.BLPW.2018.csv") %>% add_column(., station = "PIBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/RUTH.BLPW.2018.csv") %>% add_column(., station = "RUTH"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/TCBO.BLPW.2018.csv") %>% add_column(., station = "TCBO"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/TLBBS.BLPW.2018.csv") %>% add_column(., station = "TLBBS"),
                   read.csv("analysis/0_data/migration_counts/CAN/2020_01_03/TTPBRS.BLPW.2018.csv") %>% add_column(., station = "TTPBRS"))
  
  # Create a column to distinguish sub-areas at LPBO (area will be 1 for all other stations)
  dat_can$area <- 1
  dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO2")] <- 2
  dat_can$area[which(dat_can$SurveyAreaIdentifier == "LPBO3")] <- 3
  
  # Number of days with data each year at each station
  day_range_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = range)
  ndays_per_station <- aggregate(doy~YearCollected + SurveyAreaIdentifier + season, data = dat_can, FUN = length)
  years_per_station <- aggregate(YearCollected ~ SurveyAreaIdentifier + season, data = dat_can, FUN = range)
  
  # Ensure that all stations/days/years/seasons are included as data
  area_season_combinations_can <- unique(dat_can[,c("SurveyAreaIdentifier","season")])
  
  # Combine Canadian data
  dat_combined_can = data.frame()
  for (i in 1:nrow(area_season_combinations_can)){
    
    dat <- subset(dat_can, SurveyAreaIdentifier == area_season_combinations_can$SurveyAreaIdentifier[i] & season == area_season_combinations_can$season[i])
    if (nrow(dat) == 0) next
    
    min_doy <- min(dat$doy)
    max_doy <- max(dat$doy)
    min_year <- min(dat$YearCollected)
    max_year <- max(dat$YearCollected)
    
    # Create a "full" dataframe to store counts on all days (including NAs)
    dat_full <- expand.grid(YearCollected = seq(min_year,max_year),doy = seq(min_doy,max_doy))
    
    # Fill with counts
    dat_full <- merge(dat_full, dat, all.x = TRUE)
    
    # Ensure relevant data is filled in
    dat_full$SurveyAreaIdentifier = dat$SurveyAreaIdentifier[1]
    dat_full$station = dat$station[1]
    dat_full$area = dat$area[1]
    dat_full$season = dat$season[1]
    dat_full$min_doy = min_doy
    dat_full$max_doy = max_doy
    dat_full$min_year = min_year
    dat_full$max_year = max_year
    
    dat_combined_can = rbind(dat_combined_can, dat_full)
  }
  dat_combined_can$country = "CAN"
  
  # *****
  # Read/format data from Ricky Dunn (USA data)
  # *****
  dat_usa <- rbind(readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - AIMS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - BIBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - BSBO fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - BSBO spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - FBBS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - FBBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - KWRS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - MCCS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                   readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - MCCS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"))
  
  # datasets with different column names
  tmp = readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - PARC fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
  colnames(tmp) <- colnames(dat_usa)
  dat_usa <- rbind(dat_usa, tmp)
  rm(tmp)    
  
  tmp = readxl::read_xlsx("analysis/0_data/migration_counts/USA/Cleaned BLPW - CFMS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
  colnames(tmp) <- colnames(dat_usa)
  dat_usa <- rbind(dat_usa, tmp)
  rm(tmp) 
  
  # Fill net N/100 net-hr column
  dat_usa$`N/100 net-hr` = dat_usa$N/dat_usa$`Net-hrs`*100
  dat_usa = na.omit(dat_usa)
  summary(dat_usa)
  
  # Change column names
  colnames(dat_usa) = c("station","YearCollected","doy","net.hrs","ObservationCount","N.per.100.net.hrs","season")
  dat_usa$area = 1
  
  # Load analysis windows for each station and restrict data to those dates
  us_windows = read.csv("analysis/0_data/migration_counts/USA/US_station_windows.csv")
  dat_usa2 = data.frame()
  for (i in 1:nrow(us_windows)){
    x = us_windows[i,]
    station = us_windows$station[i]
    area = us_windows$area[i]
    season = us_windows$season[i]
    
    # data inside that range
    dat_usa2 = rbind(dat_usa2, subset(dat_usa, station == x$station & area == x$area & season == x$season & doy >= x$start_date & doy <= x$end_date))
    rm(x)
  }
  
  dat_usa = dat_usa2
  dat_usa = merge(dat_usa, us_windows, all.x = TRUE)
  
  area_season_combinations_usa <- unique(dat_usa[,c("station","season")])
  dat_combined_usa = data.frame()
  for (i in 1:nrow(area_season_combinations_usa)){
    dat <- subset(dat_usa, station == area_season_combinations_usa$station[i] & season == area_season_combinations_usa$season[i])
    if (nrow(dat) == 0) next
    min_doy <- min(dat$doy)
    max_doy <- max(dat$doy)
    min_year <- min(dat$YearCollected)
    max_year <- max(dat$YearCollected)
    
    # Create a "full" dataframe to store counts on all days (including NAs)
    dat_full <- expand.grid(YearCollected = seq(min_year,max_year),
                            doy = seq(min_doy,max_doy))
    
    # Fill with counts
    dat_full <- merge(dat_full, dat, all.x = TRUE)
    
    # Ensure relevant data is filled in
    dat_full$station = dat$station[1]
    dat_full$station = dat$station[1]
    dat_full$area = dat$area[1]
    dat_full$season = dat$season[1]
    dat_full$min_doy = min_doy
    dat_full$max_doy = max_doy
    dat_full$min_year = min_year
    dat_full$max_year = max_year
    
    dat_combined_usa = rbind(dat_combined_usa, dat_full)
  }
  dat_combined_usa$country = "USA"
  
  # *****
  # Combine CAN and USA data into a single dataframe
  # *****
  dat_combined = dplyr::bind_rows(dat_combined_can, dat_combined_usa)
  
  # --------------------------------------------------------------------
  # Fix station labels
  # --------------------------------------------------------------------
  dat_combined$station[dat_combined$station == "MCCS"] <- "MBO"
  # --------------------------------------------------------------------
  
  
  # Limit to data collected in year range
  dat_combined = subset(dat_combined, YearCollected >= start_year & YearCollected <= end_year)
  
  # Restrict to season of interest
  dat_combined <- subset(dat_combined, season == focal_season)
  dat_combined$site <- paste0(dat_combined$station,dat_combined$area)
  
  # Dummy variable that denotes whether a site should have a different intercept estimated for it (for LPBO)
  dat_combined$dummy_site <- (dat_combined$area > 1) %>% as.numeric()
  
  #----------
  # For any sites without recorded net hours, fill in with mean (this also fills all CMMN stations)
  dat_combined$net.hrs[which(is.na(dat_combined$net.hrs))] = mean(dat_combined$net.hrs, na.rm = TRUE)
  #----------
  
  #----------
  # Summaries of data availability at each station
  #----------
  
  # Load coordinates of cmmn stations
  station_coordinates <- read_xlsx("analysis/0_data/locations/station_locations.xlsx")
  
  # --------------------------------------------------------------------
  # Fix station labels
  # --------------------------------------------------------------------
  station_coordinates$station[station_coordinates$station == "MCCS"] <- "MBO"
  # --------------------------------------------------------------------
  
  station_data_summarized <- dat_combined %>%
    group_by(YearCollected,season,station) %>%
    summarize(mean_daily_count = mean(ObservationCount, na.rm = TRUE),
              total_count = sum(ObservationCount, na.rm = TRUE),
              ndays_observed = sum(ObservationCount > 0 & !is.na(ObservationCount))) %>%
    group_by(season,station) %>%
    summarize(mean_seasonal_total = mean(total_count, na.rm = TRUE),
              nyears_operational = length(unique(YearCollected)),
              nyears_observed = sum(ndays_observed > 0),
              nyears_observed_10count = sum(total_count >= 10),
              nyears_observed_20count = sum(total_count >= 20),
              nyears_observed_5days = sum(ndays_observed >= 5),
              nyears_20count_and_5days = sum(total_count >= 20 & ndays_observed >= 5),
              first_year = min(YearCollected),
              last_year = max(YearCollected))
  
  station_data_summarized$season <- factor(station_data_summarized$season, levels = c("Spring","Fall"))
  
  # Merge station summaries with coordinates dataframe
  station_data_summarized <- full_join(station_data_summarized, station_coordinates)
  
  # Remove stations with no season column
  station_data_summarized <- subset(station_data_summarized, !is.na(season))
  
  # Remove stations with no geographic information
  station_data_summarized <- station_data_summarized %>% subset(!is.na(lat) & !is.na(lon))
  
  # Create spatial object
  station_data_summarized_sf <-  station_data_summarized %>% st_as_sf(coords = c("lon", "lat"),crs = 4269, agr = "constant", remove = FALSE)
  
  #----------
  # Identify survey window at each station
  #----------
  
  date_ranges_start <- aggregate(doy ~ YearCollected + season + site, data = na.omit(dat_combined[,c("YearCollected","season","site","doy","ObservationCount")]), FUN = min)
  date_ranges_end <- aggregate(doy ~ YearCollected + season + site, data = na.omit(dat_combined[,c("YearCollected","season","site","doy","ObservationCount")]), FUN = max)
  colnames(date_ranges_start)[4] <- "start"
  colnames(date_ranges_end)[4] <- "end"
  date_ranges <- merge(date_ranges_start, date_ranges_end, all = TRUE)
  
  date_ranges <- date_ranges %>%
    group_by(site, season) %>%
    summarize(start = min(start, na.rm = TRUE),
              end = max(end, na.rm = TRUE))
  
  date_ranges$window_length <- date_ranges$end - date_ranges$start + 1
  
  # *****
  # FORMAT DATA FOR JAGS (note that isotope/breeding origin data is formatted later in script)
  # *****
  
  # nsite (note that 3 sub-areas at LPBO will be treated as different "sites")
  site_vec <- sort(unique(dat_combined$site))
  nsite <- length(site_vec)
  
  # nyear
  year_vec <- min(dat_combined$YearCollected):max(dat_combined$YearCollected)
  nyear <- length(year_vec)
  
  # nday
  nday <- max(date_ranges$window_length)
  
  # counts
  count_df <- dat_combined[,c("site","YearCollected","doy","ObservationCount","net.hrs","dummy_site")] %>%
    rename(site_name = site, year_abs = YearCollected, day_number = doy, count = ObservationCount, net_hrs = net.hrs) %>%
    add_column(site_number = match(.$site_name,site_vec),
               year_number = match(.$year_abs,year_vec),
               station = gsub('[[:digit:]]+', '', .$site_name)) %>%
    subset(!is.na(count))
  
  count_df$station_number <- factor(count_df$station) %>% as.numeric()
  station_names <- factor(count_df$station) %>% levels(.)
  
  jags_data <- list(nstrata = nstrata,
                    nday = max(count_df$day_number),
                    nyear = max(count_df$year_number),
                    nsite = max(count_df$site_number),
                    nstation = max(count_df$station_number),
                    nobs = nrow(count_df),
                    
                    day = count_df$day_number,
                    year = count_df$year_number,
                    site = count_df$site_number,
                    station = count_df$station_number,
                    
                    count = count_df$count,
                    net_hrs = count_df$net_hrs,
                    dummy_site = count_df$dummy_site)
  
} # End section that processes count data













# ******************************************************************************
# ******************************************************************************
# SELECT DATA FOR ANALYSIS
# ******************************************************************************
# ******************************************************************************

CAN_sites <- dat_combined %>%
  subset(country == "CAN" & season == focal_season) %>%
  dplyr::select(site) %>%
  unique()
CAN_sites <- CAN_sites$site

# ----------------------------------
# Summary of data availability
# ----------------------------------

# Only select sites with a 
count_summary <- count_df %>%
  group_by(site_name,year_abs,year_number) %>%
  summarize(max_count = max(count),
            total_count = sum(count)) %>%
  group_by(site_name) %>%
  summarize(total_count = mean(total_count),
            max_count = mean(max_count))

# ----------------------------------
# Data to include in analysis
# ----------------------------------

# Must collect an average of 10 birds per season to be included
sites_to_include <- subset(count_summary, total_count >= 10 & site_name %in% CAN_sites)

count_df_for_analysis <- subset(count_df, site_name %in% sites_to_include$site_name)
count_df_for_analysis$site_number <- as.numeric(as.factor(count_df_for_analysis$site_name))













# ******************************************************************************
# ******************************************************************************
# METHOD 1: SUM UP ANNUAL INDICES EACH YEAR AFTER CORRECTING FOR EFFORT
# ******************************************************************************
# ******************************************************************************

# ------------------
# FIT THE MODEL IN JAGS
# ------------------

# The jags script to fit the model
sink("./analysis/1_scripts/migration_model_site.jags")
cat("
    model {

  #---------------------------------------------
  # Population process
  #---------------------------------------------
  
  for (s in 1:nsite){
  
    # Model dynamics through time
    process_sd[s] ~ dunif(0,2)
    intercept[s] ~ dnorm(0,0.1)
    slope[s] ~ dnorm(0,0.1)
    for (t in 1:nyear){
      logN[s,t] ~ dnorm(intercept[s] + slope[s]*(t-1),pow(process_sd[s],-2) )
      N[s,t] <- exp(logN[s,t])
    }
    
    #---------------------------------------------
    # Parameters for within-season dynamics at each site
    #---------------------------------------------
  
    # Migration phenology at each site
    migration_phenology_mean[s] ~ dunif(1,360) # Date of peak migration
    migration_phenology_sd[s] ~ dunif(0,20)    # Width of migration period
    
    # Daily overdispersion in counts at each site (e.g., due to daily weather)
    siteday_sd[s] ~ dunif(0,2)

    
  }
  
  for (i in 1:nobs){

    mu[i] <- log(f[i]) + log(N[site[i],year[i]]) + log(net_hrs[i])
    expected_count[i] <- exp(mu[i])

    # Likelihood for counts
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[site[i]], pow(migration_phenology_sd[site[i]],-2)    ))

    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i], pow(siteday_sd[site[i]],-2))
    count[i] ~ dpois(exp(mu[i]))
    
  }

  # *********************************************************
  # Derived quantities
  # *********************************************************

  for (t in 1:nyear){
    Nsum[t] <- sum(N[1:nsite,t])
  }

}
    ",fill = TRUE)
sink()

# ----------------------------------
# Package data for jags
# ----------------------------------
year_vec <- min(count_df_for_analysis$year_abs):max(count_df_for_analysis$year_abs)
year_num <- 1:length(year_vec)

jags_data <- list(count = count_df_for_analysis$count,
                  net_hrs = count_df_for_analysis$net_hrs,
                  year = match(count_df_for_analysis$year_abs,year_vec),
                  site = count_df_for_analysis$site_number,
                  day = count_df_for_analysis$day,
                  
                  pi = pi,
                  nyear = length(year_vec),
                  nsite = length(unique(count_df_for_analysis$site_name)),
                  nobs = nrow(count_df_for_analysis))


# ----------------------------------
# Fit Model
# ----------------------------------

parameters.to.save = c("Nsum",
                       "intercept",
                       "slope",
                       "N",
                       "process_sd",
                       "siteday_sd",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       "expected_count"
)



nsamp <- 1000
nb <- 1000
nt <- 5
ni <- nb + nsamp*nt

out <- jags(data = jags_data,
            model.file = "analysis/1_scripts/migration_model_site.jags",
            parameters.to.save = parameters.to.save,
            inits = NULL,
            n.chains = 3,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = TRUE)

out$mcmc.info$elapsed.mins # 2 minutes
out

#-------------------------------------------------------------------
# Observed counts each year (overlaid with model-derived expected counts)
#-------------------------------------------------------------------

expected_vs_observed <- count_df_for_analysis %>%
  group_by(site_name,year_abs) %>%
  summarize(sum_obs = sum(count),
            sum_N = sum(count/net_hrs))

# Calculate uncertainty in expected count for each site-year combination
for (i in 1:nrow(expected_vs_observed)){
  j <- which(count_df_for_analysis$year_abs == expected_vs_observed$year_abs[i] & 
               count_df_for_analysis$site_name == expected_vs_observed$site_name[i])
  
  expected <- out$sims.list$expected_count[,j] %>% apply(.,1,sum) %>% quantile(c(0.025,0.5,0.975))
  
  expected_vs_observed$expected_q025[i] <-  expected[1]
  expected_vs_observed$expected_q50[i]  <-  expected[2]
  expected_vs_observed$expected_q975[i] <-  expected[3]
  
}

# Correlation between observed and expected counts
cor_obs_expected <- expected_vs_observed %>%
  group_by(site_name) %>%
  summarize(cor = round(cor(sum_obs,expected_q50),2),
            max = max(c(sum_obs,expected_q975)))

cor_obs_expected$cor_label <- paste0("cor = ",cor_obs_expected$cor)

plot_obs_vs_expected <- ggplot(data = expected_vs_observed)+
  geom_errorbar(aes(x = year_abs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = year_abs, y = expected_q50, col = "Expected", shape = "Expected"))+
  geom_point(aes(x = year_abs, y = sum_obs, col = "Observed", shape = "Observed"))+
  geom_text(data = cor_obs_expected, aes(x = 2025, y = 0, label = cor_label), 
            hjust = 1, size = 3)+
  facet_wrap(site_name~.)+
  scale_color_manual(values=c("gray75","black"), name = "", guide = FALSE)+
  scale_shape_manual(values=c(19,4), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Total Seasonal Count")+
  xlim(c(1990,2025))+
  ggtitle("Observed vs Expected Seasonal Total Counts")

plot_obs_vs_expected

#-------------------------------------------------------------------
# Annual indices at each station
#-------------------------------------------------------------------

site_names <- count_df_for_analysis[,c("site_name","site_number")] %>%
  unique() %>%
  arrange(site_number)

N_df <- out$sims.list$N %>%
  reshape2::melt(., varnames = list("samp","site_number","year_number"), value.name = "N") %>%
  group_by(site_number,year_number) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))

N_df$site_name <- site_names$site_name[N_df$site_number]
N_df$Year <- year_vec[N_df$year_number]

plot_fit <- ggplot()+
  geom_ribbon(data = N_df, aes(x = Year, ymin = N_q025,ymax = N_q975, col = "Expected", fill = "Expected"))+
  geom_line(data = N_df, aes(x = Year, y = N_q500), col = "gray50")+
  geom_point(data = expected_vs_observed, aes(x = year_abs, y = sum_N), shape = 4)+
  facet_wrap(site_name~.)+
  scale_color_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  scale_fill_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Sum (Birds/Net-hr) across season")+
  ggtitle("Observed vs Expected")+
  xlim(c(1990,2025))

plot_fit

#-------------------------------------------------------------------
# Total abundance
#-------------------------------------------------------------------

Nsum_df <- data.frame(Year = year_vec,
                      N_q025 = out$q2.5$Nsum,
                      N_q500 = out$q50$Nsum,
                      N_q975 = out$q97.5$Nsum)

plot_Nsum <- ggplot(data = subset(Nsum_df, Year >=1998))+
  geom_ribbon(aes(x = Year, ymin = N_q025, ymax = N_q975), alpha = 0.2)+
  geom_line(aes(x = Year, y = N_q500))+
  theme_bw()+
  scale_y_continuous(trans = "log10")

plot_obs_vs_expected
plot_fit
plot_Nsum

#-------------------------------------------------------------------
# Trend estimate
#-------------------------------------------------------------------

nsamps <- out$mcmc.info$n.samples
trend1 <- c()

for (i in 1:nsamps){
  trend1[i] <- mean(diff(log(out$sims.list$Nsum[i,])))
}

hist(trend1, breaks = 100)



# ******************************************************************************
# ******************************************************************************
# METHOD 2: FIT A MODEL THAT SHARES INFORMATION BETWEEN STATIONS
# ******************************************************************************
# ******************************************************************************

# ------------------
# FIT THE MODEL IN JAGS
# ------------------

# The jags script to fit the model
sink("./analysis/1_scripts/migration_model_site_shared.jags")
cat("
    model {

  #---------------------------------------------
  # NATIONAL POPULATION PROCESS
  #---------------------------------------------
  
  process_sd ~ dunif(0,2)
  # intercept <- 0
  # slope ~ dnorm(0,0.1)
  # for (t in 1:nyear){
  #   logNsum[t] <- slope * (t-1)
  #   Nsum[t] <- exp(logNsum[t])
  # }
  
  logNsum[1] <- 0
  for (t in 2:nyear){
    logNsum[t] ~ dnorm(logNsum[t-1],pow(process_sd,-2))
  }
  for (t in 1:nyear){Nsum[t] <- exp(logNsum[t])}
  #---------------------------------------------
  # SITE-LEVEL POPULATION PROCESS
  #---------------------------------------------
  
  for (s in 1:nsite){
  
    log_rho[s] ~ dnorm(0,0.01)
    annual_sd[s] ~ dunif(0,2)
    
    for (t in 1:nyear){
      logN[s,t] ~ dnorm(log(Nsum[t]*exp(log_rho[s])),pow(annual_sd[s],-2) )
      N[s,t] <- exp(logN[s,t])
    }
    
    
    #---------------------------------------------
    # Parameters for within-season dynamics at each site
    #---------------------------------------------
  
    # Migration phenology at each site
    migration_phenology_mean[s] ~ dunif(1,360) # Date of peak migration
    migration_phenology_sd[s] ~ dunif(0,20)    # Width of migration period
    
    # Daily overdispersion in counts at each site (e.g., due to daily weather)
    siteday_sd[s] ~ dunif(0,2)

    
  }
  
  for (i in 1:nobs){

    mu[i] <- log(f[i]) + log(N[site[i],year[i]]) + log(net_hrs[i])
    expected_count[i] <- exp(mu[i])

    # Likelihood for counts
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[site[i]], pow(migration_phenology_sd[site[i]],-2)    ))

    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i], pow(siteday_sd[site[i]],-2))
    count[i] ~ dpois(exp(mu[i]))
    
  }


}
    ",fill = TRUE)
sink()

# ----------------------------------
# Package data for jags
# ----------------------------------
year_vec <- min(count_df_for_analysis$year_abs):max(count_df_for_analysis$year_abs)
year_num <- 1:length(year_vec)

jags_data <- list(count = count_df_for_analysis$count,
                  net_hrs = count_df_for_analysis$net_hrs,
                  year = match(count_df_for_analysis$year_abs,year_vec),
                  site = count_df_for_analysis$site_number,
                  day = count_df_for_analysis$day,
                  
                  pi = pi,
                  nyear = length(year_vec),
                  nsite = length(unique(count_df_for_analysis$site_name)),
                  nobs = nrow(count_df_for_analysis))


# ----------------------------------
# Fit Model
# ----------------------------------

parameters.to.save = c("Nsum",
                       "intercept",
                       "slope",
                       "N",
                       "log_rho",
                       "annual_sd",
                       "process_sd",
                       "siteday_sd",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       "expected_count"
)



nsamp <- 1000
nb <- 1000
nt <- 5
ni <- nb + nsamp*nt

out2 <- jags(data = jags_data,
            model.file = "analysis/1_scripts/migration_model_site_shared.jags",
            parameters.to.save = parameters.to.save,
            inits = NULL,
            n.chains = 3,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = TRUE)

out2$mcmc.info$elapsed.mins
out2

#-------------------------------------------------------------------
# Expected counts/net hour each day of the year (overlaid with observed counts)
#-------------------------------------------------------------------

expected_df <- out$sims.list$E %>%
  reshape2::melt(., varnames = list("samp","day_number","site_number","year_number"), value.name = "N") %>%
  group_by(day_number,site_number,year_number) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))

expected_df$day_number <- day_seq[expected_df$day_number]

expected_df$site_name <- site_names$site_name[expected_df$site_number]
expected_df$Year <- year_vec[expected_df$year_number]

count_df_for_analysis$Year <- count_df_for_analysis$year_abs
plot_fit <- ggplot()+
  geom_point(data = count_df_for_analysis, aes(x = day_number, y = count/net_hrs), shape = 19, size =0.1)+
  geom_ribbon(data = expected_df, aes(x = day_number, ymin = N_q025,ymax = N_q975, col = "Expected", fill = "Expected"), alpha = 0.5)+
  geom_line(data = expected_df, aes(x = day_number, y = N_q500), col = "gray50")+
  facet_grid(site_name~Year, scales = "free")+
  scale_color_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  scale_fill_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Sum (Birds/Net-hr) across season")+
  ggtitle("Observed vs Expected")

plot_fit

#-------------------------------------------------------------------
# Observed counts each year (overlaid with model-derived expected counts)
#-------------------------------------------------------------------

expected_vs_observed <- count_df_for_analysis %>%
  group_by(site_name,year_abs) %>%
  summarize(sum_obs = sum(count),
            sum_N = sum(count/net_hrs))

# Calculate uncertainty in expected count for each site-year combination
for (i in 1:nrow(expected_vs_observed)){
  j <- which(count_df_for_analysis$year_abs == expected_vs_observed$year_abs[i] & 
               count_df_for_analysis$site_name == expected_vs_observed$site_name[i])
  
  expected <- out2$sims.list$expected_count[,j] %>% apply(.,1,sum) %>% quantile(c(0.025,0.5,0.975))
  
  expected_vs_observed$expected_q025[i] <-  expected[1]
  expected_vs_observed$expected_q50[i]  <-  expected[2]
  expected_vs_observed$expected_q975[i] <-  expected[3]
  
}

# Correlation between observed and expected counts
cor_obs_expected <- expected_vs_observed %>%
  group_by(site_name) %>%
  summarize(cor = round(cor(sum_obs,expected_q50),2),
            max = max(c(sum_obs,expected_q975)))

cor_obs_expected$cor_label <- paste0("cor = ",cor_obs_expected$cor)

plot_obs_vs_expected <- ggplot(data = expected_vs_observed)+
  geom_errorbar(aes(x = year_abs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = year_abs, y = expected_q50, col = "Expected", shape = "Expected"))+
  geom_point(aes(x = year_abs, y = sum_obs, col = "Observed", shape = "Observed"))+
  geom_text(data = cor_obs_expected, aes(x = 2025, y = 0, label = cor_label), 
            hjust = 1, size = 3)+
  facet_wrap(site_name~.)+
  scale_color_manual(values=c("gray75","black"), name = "", guide = FALSE)+
  scale_shape_manual(values=c(19,4), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Total Seasonal Count")+
  xlim(c(1990,2025))+
  ggtitle("Observed vs Expected Seasonal Total Counts")

plot_obs_vs_expected

#-------------------------------------------------------------------
# Annual indices at each station
#-------------------------------------------------------------------

site_names <- count_df_for_analysis[,c("site_name","site_number")] %>%
  unique() %>%
  arrange(site_number)

N_df <- out2$sims.list$N %>%
  reshape2::melt(., varnames = list("samp","site_number","year_number"), value.name = "N") %>%
  group_by(site_number,year_number) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))

N_df$site_name <- site_names$site_name[N_df$site_number]
N_df$Year <- year_vec[N_df$year_number]

plot_fit <- ggplot()+
  geom_ribbon(data = N_df, aes(x = Year, ymin = N_q025,ymax = N_q975, col = "Expected", fill = "Expected"))+
  geom_line(data = N_df, aes(x = Year, y = N_q500), col = "gray50")+
  geom_point(data = expected_vs_observed, aes(x = year_abs, y = sum_N), shape = 4)+
  facet_wrap(site_name~.)+
  scale_color_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  scale_fill_manual(values=c("gray85","black"), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Sum (Birds/Net-hr) across season")+
  ggtitle("Observed vs Expected")+
  xlim(c(1990,2025))

plot_fit

#-------------------------------------------------------------------
# Total abundance
#-------------------------------------------------------------------

Nsum_df <- data.frame(Year = year_vec,
                      N_q025 = out2$q2.5$Nsum,
                      N_q500 = out2$q50$Nsum,
                      N_q975 = out2$q97.5$Nsum)

plot_Nsum <- ggplot(data = subset(Nsum_df, Year >=1998))+
  geom_ribbon(aes(x = Year, ymin = N_q025, ymax = N_q975), alpha = 0.2)+
  geom_line(aes(x = Year, y = N_q500))+
  theme_bw()+
  scale_y_continuous(trans = "log10")

plot_obs_vs_expected
plot_fit
plot_Nsum

#-------------------------------------------------------------------
# Trend estimate
#-------------------------------------------------------------------

nsamps <- out2$mcmc.info$n.samples
trend2 <- c()

for (i in 1:nsamps){
  trend2[i] <- mean(diff(log(out2$sims.list$Nsum[i,])))
}

hist(trend2, breaks = 100)

