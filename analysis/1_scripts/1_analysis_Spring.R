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

focal_season <- "Spring"
start_year <- 2000
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


# PART 2: FORMAT BREEDING ORIGINS
{
  
  # *****
  # Breeding origin assignments
  # *****
  assignment_file = "analysis/0_data/Isotopes/isotope_assignments_2strata.xlsx"
  
  # First 7 columns contain the relevant data
  assignments <- read_xlsx(assignment_file)[,c("location","season","year","lat","lon","assigned_to_West","assigned_to_East")] %>% 
    subset(season == focal_season & year %in% year_vec)
  
  # Spatial datasets used to assign breeding origin data to each migration monitoring station
  assignments_sf <- assignments %>% st_as_sf(coords = c("lon", "lat"),crs = 4269, agr = "constant", remove = FALSE)
  station_data_summarized_sf <- station_data_summarized_sf %>% st_transform(crs(assignments_sf)) # convert station locations to same crs
  
  # *****
  # For each source of breeding origin data, determine the nearest migration station (within 300 km) and assign those
  # data to be used at that station
  # *****
  
  # An array to store assignment information at each migration monitoring station
  N_origin <- array(NA, dim = c(jags_data$nstrata,jags_data$nstation, jags_data$nyear))
  dimnames(N_origin)[[1]] <- c("East","West")
  dimnames(N_origin)[[2]] <- station_names
  dimnames(N_origin)[[3]] <- year_vec
  
  # Loop through isotope assignments, and assign them to nearest station (up to max distance of 250 km)
  
  for (i in 1:nrow(assignments_sf)){
    
    # Determine nearest station
    dists <- st_distance(assignments_sf[i,],station_data_summarized_sf) %>% as.numeric()
    if (min(dists) >= 250000) next
    
    nearest_station <- station_data_summarized_sf[which(dists == min(dists)),]
    
    assignment_year <- assignments_sf$year[i]
    
    if (sum(is.na(N_origin[,nearest_station$station,as.character(assignment_year)]))>0) N_origin[1:(dim(N_origin)[1]),nearest_station$station,as.character(assignment_year)] <- rep(0,length(N_origin[,nearest_station$station,as.character(assignment_year)]))
    
    
    N_origin["West",nearest_station$station,as.character(assignment_year)] <- N_origin["West",nearest_station$station,as.character(assignment_year)] + as.numeric(assignments[i,"assigned_to_West"])
    N_origin["East",nearest_station$station,as.character(assignment_year)] <- N_origin["East",nearest_station$station,as.character(assignment_year)] + as.numeric(assignments[i,"assigned_to_East"])
    
  }
  
  # Sample sizes
  N_station_sampled <- apply(N_origin,c(2,3), sum)
  N_station_sampled[is.na(N_station_sampled)] <- 999  # Placeholder for sample size in years with no stable isotope information
  
  # Append to jags data package
  jags_data$N_origin <- N_origin
  jags_data$N_station_sampled <- N_station_sampled
  
} # End breeding origin processing script

site_vec <- unique(count_df$site_name) %>% sort()
year_vec <- min(count_df$year_abs):max(count_df$year_abs)

#------------------------------------------------
# Fix migration parameters for several stations
#------------------------------------------------

# Fix certain migration rates to zero; stations with 0 for "rho_fix" will not have migration parameters estimated
rho_fix <- matrix(1,nrow = jags_data$nstrata, 
                  ncol = jags_data$nstation, 
                  dimnames = list(c("East","West"),station_names)) 

# Western stations cannot receive eastern birds
rho_fix["East", c("ACBO","LMBO")] <- 0 

# Eastern coastal sites cannot receive western birds
rho_fix["West", c("AIMS","MCCS")] <- 0 

#------------------------------------------------
# Add several extra quantities that are needed for JAGS script
#------------------------------------------------

jags_data$rho_fix <- rho_fix
jags_data$pi <- pi
jags_data$logN0 <- rep(0,jags_data$nstrata)

# "Station" at each site
station_site <- count_df %>% group_by(station,site_number) %>% summarize(station_number = mean(station_number))

jags_data$station_site <- station_site$station_number

# ------------------
# FIT THE MODEL IN JAGS
# ------------------

# The jags script to fit the model
sink("./analysis/1_scripts/migration_model.jags")
cat("
    model {
  
  # *********************************************************
  # Priors and likelihood
  # *********************************************************

  #---------------------------------------------
  # Model for population dynamics in each region
  #---------------------------------------------
  
  sigma_proc ~ dunif(0,2)
  tau_proc <- pow(sigma_proc,-2)
    
  for (j in 1:nstrata){
    
    logX[j,1] <- 0
    X[j,1] <- exp(logX[j,1])
    
    for (y in 2:nyear){

      logX[j,y] ~ dnorm(logX[j,y-1],tau_proc)
      X[j,y] <- exp(logX[j,y])
    
    }
  
  } # j
  
  #---------------------------------------------
  # Model for breeding origins of migrants arriving at each station
  #---------------------------------------------

  for (s in 1:nstation){

    for (j in 1:nstrata){
      rho[j,s] ~ dlnorm(0,0.25)
    }
    
    for (y in 1:nyear){
      for (j in 1:nstrata){
        M[j,s,y] <- X[j,y] * rho[j,s] * rho_fix[j,s]
      }
    }
    
  }
  
  sigma_rho ~ dunif(0,2)
  tau_rho <- pow(sigma_rho,-2)

  for (s in 1:nstation){
    for (y in 1:nyear){
      
      # Total number of migrants in year [y] at station [s]
      T[s,y] <- sum(M[1:nstrata,s,y])
      
      # Proportion of birds from each stratum
      for (j in 1:nstrata){ 
        p[j,s,y] <- M[j,s,y]/T[s,y] 
      }
      
      # Multinomial likelihood for observed breeding origins in sample of birds
      N_origin[1:nstrata,s,y] ~ dmulti(p[1:nstrata,s,y],N_station_sampled[s,y])

      T_star[s,y] ~ dlnorm(log(T[s,y]) - 1/(2*tau_rho),tau_rho)
      
    } # y
  } # s 
  
  #---------------------------------------------
  # Within-season model for migration counts
  #---------------------------------------------

  # Daily overdispersion in counts at each site (e.g., due to daily weather)
  sigma_stationday ~ dunif(0,2)
  tau_stationday <- pow(sigma_stationday,-2)
  
  migration_phenology_sd ~ dunif(0,20)
  migration_phenology_tau <- pow(migration_phenology_sd,-2)
  
  for (s in 1:nstation){
    migration_phenology_mean[s] ~ dunif(1,360)
  }
  
  # Site-level fixed effects (for LPBO that contains multiple sub-stations)
  for (k in 1:nsite){

    site_effect[k] ~ dnorm(0,0.25)
    
  }
  
  for (i in 1:nobs){
    
    mu[i] <- log(f[i]) + log(T_star[station[i],year[i]]) + log(net_hrs[i]) + site_effect[site[i]]*dummy_site[i]
    expected_count[i] <- exp(mu[i])

    # Likelihood for counts
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[station[i]], migration_phenology_tau))
    
    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday), tau_stationday) 
    count[i] ~ dpois(exp(log_lambda[i]))
    
    # *********************************************************
    # Simulate counts for posterior predictive checking
    # *********************************************************

    sim_log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday), tau_stationday) 
    sim_count[i] ~ dpois(exp(sim_log_lambda[i] ))
   
    # chi-square statistics
    X2_sim[i] <- pow(sim_count[i] - expected_count[i],2)/expected_count[i]
    X2_obs[i] <- pow(count[i] - expected_count[i],2)/expected_count[i]

  }
  
  # *********************************************************
  # Derived quantities
  # *********************************************************

  #---------------------------------------------
  # Regional composition of each annual index
  #---------------------------------------------
  
  for (y in 1:nyear){
    for (s in 1:nstation){
      for (j in 1:nstrata){
        station_composition[j,s,y] <- T_star[s,y] * p[j,s,y]
      }
    }
  }
  
}
    ",fill = TRUE)
sink()

# ------------------
# Fit Model
# ------------------

parameters.to.save = c("sigma_proc",
                       "sigma_rho",
                       "sigma_stationday",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       
                       "rho",
                       "eta",
                       "T",
                       "T_star",
                       "station_composition",
                       "expected_count",
                       "site_effect",
                       "sim_count",
                       "X",
                       "X2_sim",
                       "X2_obs"
)


inits <- NULL
nsamp <- 2000
nb <- 20000
nt <- 100
ni <- nb + nsamp*nt

# # Need to rerun with time-varying model
# out <- jags(data = jags_data,
#             model.file = "analysis/1_scripts/migration_model.jags",
#             parameters.to.save = parameters.to.save,
#             inits = inits,
#             n.chains = 3,
#             n.thin = nt,
#             n.iter = ni,
#             n.burnin = nb,
#             parallel = TRUE)
# 
# out$mcmc.info$elapsed.mins # 153 mins
# save.image(paste0(output_directory,"analysis_",focal_season,".RData"))

# ------------------
# Save/load workspace
# ------------------

# Load fitted model
load(paste0(output_directory,"analysis_",focal_season,".RData"))

# ***************************************************************
# ***************************************************************
# PART 2: ASSESS MODEL CONVERGENCE AND EFFECTIVE SAMPLE SIZE
# ***************************************************************
# ***************************************************************

#-------------------------------------------------------------------
# Assess model convergence
#-------------------------------------------------------------------

latent_parameters <- c("sigma_proc",
                       "sigma_rho",
                       "sigma_stationday",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       "log_rho_mu")

latent_states <- names(out$Rhat)[which(!(names(out$Rhat) %in% latent_parameters | names(out$Rhat) %in% c("sim_count","X2_sim","X2_obs")))]

max(unlist(out$Rhat[latent_parameters]), na.rm = TRUE)
length(unlist(out$Rhat[latent_parameters]))
mean(unlist(out$Rhat[latent_parameters]) > 1.1, na.rm = TRUE)

max(unlist(out$Rhat[latent_states]), na.rm = TRUE)
length(unlist(out$Rhat[latent_states]))
mean(unlist(out$Rhat[latent_states]) > 1.1, na.rm = TRUE)
sum(unlist(out$Rhat[latent_states]) > 1.1, na.rm = TRUE)

# Effective sample sizes
n.eff <- unlist(out$n.eff[!(names(out$Rhat) %in% c("sim_count","X2_sim","X2_obs"))])
n.eff[n.eff > 1 & n.eff <= 1000] # Parameters with fewer than 1000 samples
min(n.eff[n.eff>1])
mean(n.eff[n.eff>1]<1000)

# ***************************************************************
# ***************************************************************
# PART 3: SUMMARIZE DATA AVAILABILITY ACROSS THE STUDY REGION
# ***************************************************************
# ***************************************************************
station_coordinates <- read_xlsx("analysis/0_data/locations/station_locations.xlsx")

# Summary of migration count availability at each station
station_summary <- count_df %>%
  rename(Year = year_abs, Station = station) %>%
  group_by(Station,Year) %>%
  summarize(Total_Count = sum(count),
            Min_Day = min(day_number),
            Max_Day = max(day_number),
            n_days = length(unique(day_number))) %>%
  group_by(Station) %>%
  summarize(n_Years = length(unique(Year)),
            First_Year = min(Year),
            Last_Year = max(Year),
            Year_Range = paste0(min(Year)," - ", max(Year)),
            min_Count = min(Total_Count),
            max_Count = max(Total_Count),
            mean_Count = round(mean(Total_Count),1)) %>%
  left_join(station_coordinates, by = c("Station" = "station")) %>%
  dplyr::select(Station,name,country,lat,lon,Year_Range,mean_Count,min_Count,max_Count) %>%
  rename("Station code" = Station,
         "Station name" = name,
         "Country" = country,
         "Lat" = lat,
         "Lon" = lon,
         "Year range" = Year_Range,
         "Mean annual count" = mean_Count,
         "Min annual count" = min_Count,
         "Max annual count" = max_Count) %>%
  mutate(Lat = round(Lat,1),
         Lon = round(Lon,1))

write.csv(station_summary, file = paste0(output_directory,"/tables/",focal_season,"_station_summary.csv"),row.names = FALSE)

# ***************************************************************
# ***************************************************************
# PART 4: PLOT RAW DATA
# ***************************************************************
# ***************************************************************
dsn = paste0(system.file("maps", package = "bbsBayes"),"/BBS_USGS_strata.shp")
BBS_strata_boundaries <- sf::read_sf(dsn = dsn)

# -----------------------------------------
# Map of study sites and which ones have isotopes available
# -----------------------------------------

station_assignments <- jags_data$N_origin %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, Year = Var3, n = value)
station_assignments$Stratum = factor(station_assignments$Stratum, levels = c("West","East"))


BCR_poly <- strata$BCR_poly
xlim <- c(-162, -56)
ylim <- c(35, 70)

# Add indicator for assignment information available
station_data_summarized_sf$assignments <- "No"
station_data_summarized_sf$assignments[station_data_summarized_sf$station %in% na.omit(station_assignments)$Station] <- "Yes"

origins_known <- colSums(rho_fix == 0)>0
origins_known <- origins_known[which(origins_known)]
station_data_summarized_sf$assignments[which(station_data_summarized_sf$station %in% names(origins_known))] <- "Yes"

# Convert to AEA 
strata_sf <- strata_sf %>% st_transform(st_crs(BBS_strata_boundaries))

strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))
strata_fig <- ggplot(data = strata_sf) +
  geom_sf(data = BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7)+
  xlab("")+
  ylab("")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", size = 0.5))+
  geom_sf(data = station_data_summarized_sf ,size = 1, alpha = 1, shape = 19,col = "black")+
  geom_sf_label_repel(data = station_data_summarized_sf,
                      aes(label = station),
                      max.overlaps = 50,
                      force = 1.2,
                      size = 1.2, 
                      alpha = 0.8,
                      min.segment.length = 0)+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = FALSE) +
  ggtitle("Pre-breeding migration")#+
  #coord_sf(xlim = xlim, ylim = ylim)
#strata_fig

png(file = paste0(output_directory,"figures/Manuscript_Fig1_Station_Map_",focal_season,".png"), units = "in", width = 5, height = 4, res = 600)
strata_fig 
dev.off()

save(strata_fig, file = paste0("analysis/2_output/Figures_Main_Text/",focal_season,"_map.RData"))

#-------------------------------------------------------------------
# Observed counts each year (overlaid with model-derived expected counts)
#-------------------------------------------------------------------

expected_vs_observed <- count_df %>%
  group_by(station,year_abs) %>%
  summarize(sum_obs = sum(count))

# Calculate uncertainty in expected count for each station-year combination
for (i in 1:nrow(expected_vs_observed)){
  j <- which(count_df$year_abs == expected_vs_observed$year_abs[i] & 
               count_df$station == expected_vs_observed$station[i])
  
  expected <- out$sims.list$sim_count[,j] %>% apply(.,1,sum) %>% quantile(c(0.025,0.5,0.975))
  
  expected_vs_observed$expected_q025[i] <-  expected[1]
  expected_vs_observed$expected_q50[i]  <-  expected[2]
  expected_vs_observed$expected_q975[i] <-  expected[3]
  
}

# Correlation between observed and expected counts
cor_obs_expected <- expected_vs_observed %>%
  group_by(station) %>%
  summarize(cor = round(cor(sum_obs,expected_q50),2),
            max = max(c(sum_obs,expected_q975)))

cor_obs_expected$cor_label <- paste0("cor = ",cor_obs_expected$cor)

plot_obs_vs_expected <- ggplot(data = expected_vs_observed)+
  geom_errorbar(aes(x = year_abs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = year_abs, y = expected_q50, col = "Expected", shape = "Expected"))+
  geom_point(aes(x = year_abs, y = sum_obs, col = "Observed", shape = "Observed"))+
  geom_text(data = cor_obs_expected, aes(x = 2000, y = max*1.1, label = cor_label), hjust = 0, size = 3)+
  facet_wrap(station~., scales = "free_y")+
  scale_color_manual(values=c("gray75","black"), name = "", guide = FALSE)+
  scale_shape_manual(values=c(19,4), name = "", guide = FALSE)+
  xlab("Year")+
  ylab("Total Seasonal Count")+
  ggtitle("Observed vs Expected Seasonal Total Counts\n\nPre-breeding migration")
plot_obs_vs_expected

png(file = paste0(output_directory,"figures/Appendix_S1_migration_counts.png"), units = "in", width = 8, height = 8, res = 600)
plot_obs_vs_expected
dev.off()

plot_obs_vs_expected2 <- ggplot(data = expected_vs_observed)+
  geom_abline(slope=1,intercept=0)+
  geom_errorbar(aes(x = sum_obs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = sum_obs, y = expected_q50, col = "Expected", shape = "Expected"))+
  geom_text(data = cor_obs_expected, aes(x = 2000, y = max*1.1, label = cor_label), hjust = 0, size = 3)+
  facet_wrap(station~., scales = "free")+
  scale_color_manual(values=c("gray75","black"), name = "", guide = FALSE)+
  scale_shape_manual(values=c(19,4), name = "", guide = FALSE)+
  xlab("Observed Total Count")+
  ylab("Expected Total Count")+
  ggtitle("Observed vs Expected Seasonal Total Counts\n\nPre-breeding migration")
plot_obs_vs_expected2

png(file = paste0(output_directory,"figures/Appendix_S1_obs_vs_expected.png"), units = "in", width = 8, height = 8, res = 600)
plot_obs_vs_expected2
dev.off()

#-------------------------------------------------------------------
# Breeding origin assignments at each station
#-------------------------------------------------------------------

station_assignment_plot <- ggplot(station_assignments, 
                                  aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~., scales = "free")+
  xlim(c(range(station_assignments$Year)))+
  ggtitle(focal_season)+
  ylab("Number of bird samples analyzed")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_assignment_plot

png(file = paste0(output_directory,"figures/Appendix_S2_breeding_origins.png"), units = "in", width = 8, height = 6, res = 600)
station_assignment_plot
dev.off()

# ***************************************************************
# ***************************************************************
# PART 5: PLOT MODEL ESTIMATES
# ***************************************************************
# ***************************************************************

# -------------------------------------------------------------------
# Time series of regional population indices
# -------------------------------------------------------------------

X_summary <- out$sims.list$X %>%
  reshape2::melt() %>%
  rename(sample = Var1, stratum_number = Var2, year_number = Var3, X = value) %>%
  group_by(stratum_number,year_number) %>%
  summarize(X_mean = mean(X),
            X_lcl = quantile(X,0.025),
            X_ucl = quantile(X,0.975))
X_summary$Year = year_vec[X_summary$year_number]
X_summary$Stratum = c("East","West")[X_summary$stratum_number] %>% factor(levels = c("West","East"))

regional_trajectory_plot <- ggplot(X_summary,
                                   aes(x = Year, 
                                       y = X_mean, 
                                       ymin = X_lcl, 
                                       ymax = X_ucl, 
                                       fill = Stratum, 
                                       col = Stratum))+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 2)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  facet_grid(Stratum~.)+
  ylab("Regional index")+
  xlab("Year")+
  ggtitle("Regional trajectory\n\nPre-breeding migration")

regional_trajectory_plot

png(file = paste0(output_directory,"figures/Appendix_S3_regional_index.png"), units = "in", width = 6, height = 6, res = 600)
regional_trajectory_plot
dev.off()

# -------------------------------------------------------------------
# Estimates of change since 2000
# -------------------------------------------------------------------

change_since_2000 <- data.frame()

for (j in 1:jags_data$nstrata){
  for (t in 1:jags_data$nyear){
    
    log_change <- log(out$sims.list$X[,j,t]/out$sims.list$X[,j,1])
    percent_change <- 100 * (out$sims.list$X[,j,t]-out$sims.list$X[,j,1])/out$sims.list$X[,j,1]
    
    change_since_2000 <- rbind(change_since_2000,
                               data.frame(stratum_number = j,
                                          year_number = t,
                                          log_change_q025 = quantile(log_change,0.025),
                                          log_change_q50 = quantile(log_change,0.5),
                                          log_change_q975 = quantile(log_change,0.975),
                                          
                                          percent_change_q025 = quantile(percent_change,0.025),
                                          percent_change_q50 = quantile(percent_change,0.5),
                                          percent_change_q975 = quantile(percent_change,0.975),
                                          
                                          prob_decline = mean(log_change<0)
                               ))
  }
}


change_since_2000$Year = year_vec[change_since_2000$year_number]
change_since_2000$Stratum = c("East","West")[change_since_2000$stratum_number] %>% factor(levels = c("West","East"))

head(change_since_2000)

# Convert log-scale change to percent change using: 100*(exp(log_change)-1)
y_axis_breaks <- c(0,100,300,900)
log_breaks <- log(y_axis_breaks/100+1)
log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
y_axis_breaks <- 100*(exp(log_breaks)-1)

y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
y_axis_labels <- paste0(y_axis_labels," %")

plot_change_since_2000 <- ggplot(change_since_2000,
                                 aes(x = Year, 
                                     y = log_change_q50, 
                                     ymin = log_change_q025, 
                                     ymax = log_change_q975, 
                                     fill = Stratum, 
                                     col = Stratum))+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 2)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  
  facet_grid(Stratum~.)+
  ylab("Percent change since 2000")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("Estimated percent change relative to 2000\n\nPre-breeding migration")

plot_change_since_2000

png(file = paste0(output_directory,"figures/Appendix_S4_change_since_2000.png"), units = "in", width = 6, height = 6, res = 600)
plot_change_since_2000
dev.off()

# -------------------------------------------------------------------
# Estimates of change since 2008
# -------------------------------------------------------------------

change_since_2008 <- data.frame()

for (j in 1:jags_data$nstrata){
  for (t in (jags_data$nyear-10):jags_data$nyear){
    
    log_change <- log(out$sims.list$X[,j,t]/out$sims.list$X[,j,jags_data$nyear-10])
    percent_change <- 100 * (out$sims.list$X[,j,t]-out$sims.list$X[,j,jags_data$nyear-10])/out$sims.list$X[,j,jags_data$nyear-10]
    
    change_since_2008 <- rbind(change_since_2008,
                               data.frame(stratum_number = j,
                                          year_number = t,
                                          log_change_q025 = quantile(log_change,0.025),
                                          log_change_q50 = quantile(log_change,0.5),
                                          log_change_q975 = quantile(log_change,0.975),
                                          
                                          percent_change_q025 = quantile(percent_change,0.025),
                                          percent_change_q50 = quantile(percent_change,0.5),
                                          percent_change_q975 = quantile(percent_change,0.975),
                                          
                                          prob_decline = mean(log_change<0)
                               ))
  }
}


change_since_2008$Year = year_vec[change_since_2008$year_number]
change_since_2008$Stratum = c("East","West")[change_since_2008$stratum_number] %>% factor(levels = c("West","East"))

head(change_since_2008)

# Convert log-scale change to percent change using: 100*(exp(log_change)-1)

y_axis_breaks <- c(0,100,300,900)
log_breaks <- log(y_axis_breaks/100+1)
log_breaks <- c(-log_breaks,log_breaks) %>% unique() %>% sort()
y_axis_breaks <- 100*(exp(log_breaks)-1)

y_axis_labels <- y_axis_breaks %>% round() %>% paste0()
y_axis_labels[which(y_axis_breaks>0)] <- paste0("+",y_axis_labels[which(y_axis_breaks>0)])
y_axis_labels <- paste0(y_axis_labels," %")

plot_change_since_2008 <- ggplot(change_since_2008,
                                 aes(x = Year, 
                                     y = log_change_q50, 
                                     ymin = log_change_q025, 
                                     ymax = log_change_q975, 
                                     fill = Stratum, 
                                     col = Stratum))+
  geom_ribbon(alpha = 0.2, col = "transparent")+
  geom_line(linewidth = 2)+
  scale_color_manual(values = strata_colours, guide = FALSE)+
  scale_fill_manual(values = strata_colours, guide = FALSE)+
  
  scale_y_continuous(breaks = log_breaks, labels = y_axis_labels)+
  
  facet_grid(Stratum~.)+
  ylab("Percent change since 2008")+
  xlab("Year")+
  geom_hline(yintercept = 0, linetype = 2)+
  ggtitle("Estimated percent change relative to 2008\n\nPre-breeding migration")

plot_change_since_2008

png(file = paste0(output_directory,"figures/Appendix_S5_change_since_2008.png"), units = "in", width = 6, height = 6, res = 600)
plot_change_since_2008
dev.off()

#-------------------------------------------------------------------
# Station-level composition each year
#-------------------------------------------------------------------

station_composition_full <- out$sims.list$station_composition %>% 
  reshape2::melt() %>%
  rename(samp = Var1, stratum_number = Var2, station_number = Var3, year_number = Var4, index = value)
station_composition_full$Station = station_names[station_composition_full$station_number]
station_composition_full$Year = year_vec[station_composition_full$year_number]
station_composition_full$Stratum = c("East","West")[station_composition_full$stratum_number]

station_composition <- station_composition_full %>%
  group_by(Station,Year,Stratum) %>%
  summarize(index_mean = mean(index),
            index_q50 = quantile(index,0.5),
            index_q025 = quantile(index,0.025),
            index_q975 = quantile(index,0.975))
station_composition$Stratum = factor(station_composition$Stratum, levels = c("West","East"))

# Determine the years in which counts were available
station_years_with_counts <- count_df %>%
  group_by(station,year_abs) %>%
  summarize(counts_available = sum(count)>0,
            sum_count = sum(count),
            sum_net_hrs = sum(net_hrs),
            sum_obs_index = sum(count)/sum(net_hrs),
            index2 = sum(count/net_hrs)) %>%
  ungroup() %>%
  rename(Station = station, Year = year_abs)

station_composition <- full_join(station_composition, station_years_with_counts)
station_composition$counts_available[!is.na(station_composition$counts_available)] <- "Yes"
station_composition$counts_available[is.na(station_composition$counts_available)] <- "No"

station_composition_plot <- ggplot(station_composition, aes(x = Year, y = index_q50, ymin = index_q025, ymax = index_q975, fill = Stratum, alpha = counts_available)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  
  scale_alpha_manual(values=c(0.2,1), guide = "none")+
  facet_wrap(Station~., scales = "free")+
  ggtitle("Esimates of annual station composition\n\nPre-breeding migration")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_composition_plot

png(file = paste0(output_directory,"figures/Appendix_S6_station_composition.png"), units = "in", width = 6, height = 6, res = 600)
station_composition_plot
dev.off()

#-------------------------------------------------------------------
# Trend estimates within each stratum, and nationally
#-------------------------------------------------------------------

trend_fn <- function(N2,N1,year_interval){
  percent_change = (N2-N1)/N1 * 100
  trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
  return(trend)
}

trend_df <- data.frame()
for (s in 1:jags_data$nstrata){
  for (start_year in 1:jags_data$nyear){
    for (end_year in 1:jags_data$nyear){
      if (end_year <= start_year) next
      
      trend_est <- trend_fn(N2 = out$sims.list$X[,s,end_year], N1 = out$sims.list$X[,s,start_year], year_interval = end_year - start_year)
      
      trend_df <- rbind(trend_df, data.frame(stratum_number = s,
                                             start_year = start_year,
                                             end_year = end_year,
                                             trend_est_q50 = quantile(trend_est,0.5),
                                             trend_est_q025 = quantile(trend_est,0.025),
                                             trend_est_q975 = quantile(trend_est,0.975),
                                             prob_decline = round(mean(trend_est<0),2)))
    }
  }
  
}
trend_df <- trend_df %>% arrange(start_year)
trend_df$start_year_abs <- year_vec[trend_df$start_year]
trend_df$end_year_abs <- year_vec[trend_df$end_year]
trend_df$Stratum <- c("East","West")[trend_df$stratum_number]


trend_df$'Trend Period' <- paste0(trend_df$start_year_abs," to ", trend_df$end_year_abs)
trend_df$'Trend Length (Yrs)' <- trend_df$end_year_abs - trend_df$start_year_abs
trend_df$'Trend Estimate' <- paste0(round(trend_df$trend_est_q50,1)," (",round(trend_df$trend_est_q025,1)," to ",round(trend_df$trend_est_q975,1),")")
trend_df$'P(Decline)' <- trend_df$prob_decline
trend_df$'Width of 95% CRI' <- round(trend_df$trend_est_q975 - trend_df$trend_est_q025,2)

trend_summary <- trend_df %>%
  subset(end_year_abs == max(trend_df$end_year_abs) & trend_df$'Trend Length' %in% c(10,max(trend_df$'Trend Length'))) %>%
  dplyr::select(Stratum,'Trend Period','Trend Length (Yrs)',"Trend Estimate",'Width of 95% CRI','P(Decline)') %>%
  arrange('Trend Length (Yrs)')
trend_summary

write.csv(trend_summary, file = paste0(output_directory,"/tables/",focal_season,"_trend_summary.csv"),row.names = FALSE)

#-------------------------------------------------------------------
# Estimate relative abundance within each stratum
#-------------------------------------------------------------------
require(ebirdst)
require(terra)

# Download eBird relative abundance map

# Need to set access key (only once) before downloading ranges
# ebirdst::set_ebirdst_access_key()

# Ensure path is correctly set
usethis::edit_r_environ()

# This should read:
# EBIRDST_KEY='ntqm1ha68fov'
# EBIRDST_DATA_DIR='D:/Working_Files/1_Projects/Landbirds/BLPW-migration-trends/0_data/eBird/'

ebirdst_download("Blackpoll Warbler") 

# #-------------------------------------------------------------------
# # Station-level indices each year
# #-------------------------------------------------------------------
# 
# station_indices_full <- out$sims.list$T %>% 
#   reshape2::melt() %>%
#   rename(samp = Var1, 
#          station_number = Var2, 
#          year_number = Var3, 
#          index = value)
# 
# station_indices_full$Station = station_names[station_indices_full$station_number]
# station_indices_full$Year = year_vec[station_indices_full$year_number]
# 
# station_indices <- station_indices_full %>%
#   group_by(Station,Year) %>%
#   summarize(index_mean = mean(index),
#             index_q50 = quantile(index,0.5),
#             index_q025 = quantile(index,0.025),
#             index_q975 = quantile(index,0.975))
# 
# station_indices <- full_join(station_indices, station_years_with_counts)
# 
# station_index_plot <- ggplot(data = na.omit(station_indices), 
#                              aes(x = Year, 
#                                  y = index_mean, 
#                                  ymin = index_q025, 
#                                  ymax = index_q975)) +
#   geom_point()+
#   geom_errorbar(width=0)+
#   facet_wrap(Station~., scales = "free_y")
# 
# station_index_plot
# 
# #write.csv(station_indices,paste0(output_directory,"tables/station_indices_Spring.csv"),row.names = FALSE)
# 
# # png(file = paste0(output_directory,"figures/Results_Station_Indices.png"), units = "in", width = 8, height = 5, res = 600)
# # station_index_plot
# # dev.off()
# 
# # -------------------------------------------------------------------
# # Calculate study-wide and 10-year trends within each stratum
# # -------------------------------------------------------------------
# 
# trend_fn <- function(N2,N1,year_interval){
#   percent_change = (N2-N1)/N1 * 100
#   trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
#   return(trend)
# }
# 
# trend_East_full <- trend_fn(N2 = out$sims.list$X[,1,19],
#                             N1 = out$sims.list$X[,1,1],
#                             year_interval = 19-1)
# quantile(trend_East_full,c(0.025,0.5,0.975))
# mean(trend_East_full<0)
# 
# trend_West_full <- trend_fn(N2 = out$sims.list$X[,2,19],
#                             N1 = out$sims.list$X[,2,1],
#                             year_interval = 19-1)
# quantile(trend_West_full,c(0.025,0.5,0.975))
# mean(trend_West_full<0)
# 
# trend_East_10 <- trend_fn(N2 = out$sims.list$X[,1,19],
#                             N1 = out$sims.list$X[,1,9],
#                             year_interval = 19-9)
# quantile(trend_East_10,c(0.025,0.5,0.975))
# mean(trend_East_10<0)
# 
# trend_West_10 <- trend_fn(N2 = out$sims.list$X[,2,19],
#                             N1 = out$sims.list$X[,2,9],
#                             year_interval = 19-9)
# quantile(trend_West_10,c(0.025,0.5,0.975))
# mean(trend_West_10<0)
# 
# 
# 
# # -------------------------------------------------------------------
# # Calculate regional and national abundance in each year using BAM's bootstrap replicates
# # -------------------------------------------------------------------
# 
# # Specify year represented by BAM density raster
# bam_year <- 2011
# bam_year_index <- which(year_vec == bam_year)
# bam_bootstrap_files <- list.files("analysis/0_data/bam_density_raster/bam_bootstrap/")
# 
# N_strata_df_total <- data.frame()
# for (i in 1:length(bam_bootstrap_files)){
#   
#   # load raster
#   bam <- raster(paste0("analysis/0_data/bam_density_raster/bam_bootstrap/",bam_bootstrap_files[i]))
#   
#   # Crop to stratum boundaries
#   strata_sp <- strata_sf %>% st_transform(crs = projection(bam)) %>% as('Spatial')
#   bam <- crop(bam,strata_sp) %>% mask(strata_sp)
#   
#   # Extract cumulative abundance in each stratum
#   stratum_N <- rep(NA,jags_data$nstrata)
#   for (r in 1:nstrata){
#     bam_stratum <- crop(bam,strata_sp[r,]) %>% mask(strata_sp[r,]) 
#     stratum_N[r] <- sum(values(bam_stratum),na.rm = TRUE) * 100 # Multiply by 100 because each pixel is 100 ha (and bam is density in males/ha)
#   }
#   
#   # Empty array to store total abundances
#   logN_array <- array(NA,c(out$mcmc.info$n.samples,jags_data$nyear,jags_data$nstrata))
#   
#   # Loop through strata
#   for (j in 1:nstrata){
#     
#     # BAM estimate for this stratum in the year of the BAM estimate
#     logN_array[,bam_year_index,j] <- log(stratum_N[j])
#     
#     # Work forwards from the year of the BAM estimate
#     if (bam_year_index < jags_data$nyear){
#       for (y in bam_year_index:(jags_data$nyear-1)){
#         
#         # Posterior samples of annual growth rate for this year
#         annual_growth <- log(out$sims.list$X[,j,y+1]/out$sims.list$X[,j,y])
#         logN_array[,y+1,j] <- logN_array[,y,j] + annual_growth
#       }
#     }
#     
#     # Work backwards from the year of the BAM estimate
#     if (bam_year_index > 1 ){
#       for (y in bam_year_index:2){
#         
#         # Posterior samples of annual growth rate for this year
#         annual_growth <- log(out$sims.list$X[,j,y]/out$sims.list$X[,j,y-1])
#         logN_array[,y-1,j] <- logN_array[,y,j] - annual_growth
#       }
#     }
#   }
#   
#   N_strata_df <- melt(logN_array) %>% 
#     rename(mcmc = Var1, year = Var2, stratum = Var3, logN = value) %>%
#     add_column(bootstrap_rep = bam_bootstrap_files[i])
#   N_strata_df$N <- exp(N_strata_df$logN)
#   
#   N_strata_df_total <- rbind(N_strata_df_total,N_strata_df)
#   
#   print(i)
#   
# }
# 
# N_strata_df_total$year_abs <- year_vec[N_strata_df_total$year]
# 
# # Summarize estimates
# N_strata_estimates <- N_strata_df_total %>% 
#   group_by(stratum,year) %>%
#   summarize(N_q025 = quantile(N,0.025),
#             N_q500 = quantile(N,0.500),
#             N_q975 = quantile(N,0.975)
#   )
# 
# N_strata_estimates$year_abs <- year_vec[N_strata_estimates$year]
# N_strata_estimates$stratum_name <- dimnames(jags_data$N_origin)[[1]][N_strata_estimates$stratum]
# # save(N_strata_estimates, file = paste0(output_directory,"/R_objects/",focal_season,"_N_strata_estimates.RData"))
# 
# # National abundance
# N_national_estimates_boot <- N_strata_df_total %>% 
#   group_by(year,mcmc,bootstrap_rep) %>%
#   
#   # Sum across all strata
#   summarize(N = sum(N))
# 
# N_national_estimates <- N_national_estimates_boot %>%
#   group_by(year) %>%                    
#   
#   # Summarize across BAM bootstrap replicates
#   summarize(N_q025 = quantile(N,0.025),
#             N_q500 = quantile(N,0.500),
#             N_q975 = quantile(N,0.975))
# N_national_estimates$year_abs <- year_vec[N_national_estimates$year]
# # save(N_national_estimates, file = paste0(output_directory,"/R_objects/",focal_season,"_N_national_estimates.RData"))
# 
# #-------------------------------------------------------------------
# # Calculate trends and percent change across specific intervals
# #-------------------------------------------------------------------
# national_trend_df <- data.frame()
# 
# # 10-year trends and study-wide trends
# for (year_interval in c(10,length(unique(N_strata_df_total$year))-1)){
#   for (year in 1:(max(N_national_estimates_boot$year)-1)){
#     for (bootrep in unique(N_national_estimates_boot$bootstrap_rep)){
#       
#       start_year_pchange <- year
#       end_year_pchange <- year+year_interval
#       
#       if (!(end_year_pchange %in% N_national_estimates_boot$year)) next
#       
#       N1 <- subset(N_national_estimates_boot,year == start_year_pchange & bootstrap_rep == bootrep)$N
#       N2 <- subset(N_national_estimates_boot, year == end_year_pchange & bootstrap_rep == bootrep)$N
#       
#       percent_change = (N2-N1)/N1 * 100
#       trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
#       
#       tmp <- data.frame(mcmc = 1:length(trend),
#                         end_year_pchange = end_year_pchange,
#                         year_interval = year_interval,
#                         percent_change = percent_change,
#                         trend = trend)
#       
#       national_trend_df <- rbind(national_trend_df, tmp)
#     } # bootrep
#     
#     print(year)
#   } # Start year
# } # Year interval
# 
# national_trend_df$stratum_name <- "National"
# national_trend_df$end_year <- year_vec[national_trend_df$end_year_pchange]
# national_trend_df$start_year <- national_trend_df$end_year - national_trend_df$year_interval
# 
# 
# # Summarize credible intervals across bootstrap replicates
# national_trend_df_summary <- national_trend_df %>%
#   group_by(stratum_name,start_year,end_year, year_interval) %>%
#   summarize(percent_change_mean = mean(percent_change),
#             percent_change_q025 = quantile(percent_change,0.025),
#             percent_change_q975 = quantile(percent_change,0.975),
#             
#             trend_mean = mean(trend),
#             trend_q025 = quantile(trend,0.025),
#             trend_q975 = quantile(trend,0.975),
#             
#             prob_increase = mean(percent_change > 0),
#             prob_30_decline = mean(percent_change <= -30),
#             prob_50_decline = mean(percent_change <= -50))
# 
# 
# # Stratum-level analysis
# stratum_trend_df <- data.frame()
# for (year_interval in c(10,length(unique(N_strata_df_total$year))-1)){
#   for (year in 1:(max(N_strata_df_total$year)-1)){
#     for (current_stratum in 1:nstrata){
#       
#       start_year_pchange <- year
#       end_year_pchange <- year+year_interval
#       
#       if (!(end_year_pchange %in% N_strata_df_total$year)) next
#       
#       N1 <- subset(N_strata_df_total,year == start_year_pchange & bootstrap_rep == N_national_estimates_boot$bootstrap_rep[1] & stratum == current_stratum)$N
#       N2 <- subset(N_strata_df_total, year == end_year_pchange & bootstrap_rep == N_national_estimates_boot$bootstrap_rep[1] & stratum == current_stratum)$N
#       
#       percent_change = (N2-N1)/N1 * 100
#       trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
#       
#       tmp <- data.frame(mcmc = 1:length(r),
#                         stratum = current_stratum,
#                         end_year_pchange = end_year_pchange,
#                         year_interval = year_interval,
#                         percent_change = percent_change,
#                         trend = trend)
#       
#       stratum_trend_df <- rbind(stratum_trend_df, tmp)
#       
#     } # stratum
#     
#     print(year)
#   } # Start year
# } # Year interval
# 
# stratum_trend_df$stratum_name <- dimnames(jags_data$N_origin)[[1]][stratum_trend_df$stratum]
# stratum_trend_df$end_year <- year_vec[stratum_trend_df$end_year_pchange]
# stratum_trend_df$start_year <- stratum_trend_df$end_year - stratum_trend_df$year_interval
# 
# # Summarize credible intervals across bootstrap replicates
# stratum_trend_df_summary <- stratum_trend_df %>%
#   group_by(stratum_name,start_year,end_year, year_interval) %>%
#   summarize(percent_change_mean = mean(percent_change),
#             percent_change_q025 = quantile(percent_change,0.025),
#             percent_change_q975 = quantile(percent_change,0.975),
#             
#             trend_mean = mean(trend),
#             trend_q025 = quantile(trend,0.025),
#             trend_q975 = quantile(trend,0.975),
#             
#             prob_increase = mean(percent_change > 0),
#             prob_30_decline = mean(percent_change <= -30),
#             prob_50_decline = mean(percent_change <= -50))
# 
# 
# trend_df_summary <- bind_rows(national_trend_df_summary,stratum_trend_df_summary) %>%
#   add_column(season = focal_season) %>%
#   relocate(stratum_name,season)
# 
# # Rearrange columns
# trend_df_summary <- trend_df_summary %>% 
#   ungroup() %>%
#   dplyr::select(stratum_name,season,start_year, end_year,
#                 trend_mean,trend_q025,trend_q975,
#                 percent_change_mean,percent_change_q025,percent_change_q975,
#                 prob_increase,prob_30_decline,prob_50_decline)
# 
# # SAVE RELEVANT OUTPUT / SUMMARIES IN TABLE FORMAT
# # write.csv(trend_df_summary,paste0(output_directory,"tables/trend_and_change_estimates.csv"),row.names = FALSE)
# 
# # *****************************************************
# # *****************************************************
# # PART 3: GENERATE SOME PLOTS OF DATA AND RESULT SUMMARIES
# # *****************************************************
# # *****************************************************
# 
# # Trajectory over last 10 years
# N_strata_estimates$stratum_name = factor(N_strata_estimates$stratum_name, levels = c("West","East"))
# 
# N_strata_10yr_plot <- ggplot(subset(N_strata_estimates,year_abs %in% ((end_year-10):end_year))) +
#   geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975, fill = stratum_name), alpha = 0.3)+
#   geom_line(aes(x = year_abs, y = N_q500, col = stratum_name))+
#   scale_fill_manual(values = strata_colours, name = "Stratum")+
#   scale_color_manual(values = strata_colours, name = "Stratum")+
#   xlab("Year")+
#   ylab("Abundance")+
#   facet_grid(.~stratum_name)+
#   scale_y_continuous(label=comma)+
#   scale_x_continuous(breaks = seq(2008,2018,length.out = 3))+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     panel.spacing = unit(1, "lines"),
#     legend.position=c(.93,.75))+
#   coord_cartesian(ylim = c(0,75000000))
# N_strata_10yr_plot
# 
# # png(file = paste0(output_directory,"figures/Trajectory_Strata_10yr.png"), units = "in", width = 6, height = 5, res = 600)
# # N_strata_10yr_plot
# # dev.off()
# 
# # Trajectory over entire time series
# N_strata_18yr_plot <- ggplot(N_strata_estimates) +
#   geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975, fill = stratum_name), alpha = 0.3)+
#   geom_line(aes(x = year_abs, y = N_q500, col = stratum_name))+
#   scale_fill_manual(values = strata_colours, name = "Stratum", guide = FALSE)+
#   scale_color_manual(values = strata_colours, name = "Stratum", guide = FALSE)+
#   xlab("Year")+
#   ylab("Abundance")+
#   facet_grid(stratum_name~.)+
#   scale_y_continuous(label=comma)+
#   theme(
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     panel.spacing = unit(1, "lines"),
#     legend.position=c(.93,.75))+
#   coord_cartesian(ylim = c(0,75000000))
# 
# # png(file = paste0(output_directory,"figures/Trajectory_Strata_18yr.png"), units = "in", width = 6, height = 5, res = 600)
# # N_strata_18yr_plot
# # dev.off()
# 
# N_national_10yr_plot <- ggplot(subset(N_national_estimates,year_abs %in% ((end_year-10):end_year))) +
#   geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975), alpha = 0.3)+
#   geom_line(aes(x = year_abs, y = N_q500))+
#   xlab("Year")+
#   ylab("National Abundance")+
#   ggtitle("National trajectory")+
#   scale_x_continuous(breaks = seq(2008,2018,length.out = 3))+
#   scale_y_continuous(label=comma)+
#   coord_cartesian(ylim = c(0,75000000))
# 
# # png(file = paste0(output_directory,"figures/Trajectory_10yr.png"), units = "in", width = 6, height = 5, res = 600)
# # N_national_10yr_plot
# # dev.off()
# 
# N_national_18yr_plot <- ggplot(N_national_estimates) +
#   geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975), alpha = 0.3)+
#   geom_line(aes(x = year_abs, y = N_q500))+
#   xlab("Year")+
#   ylab("National Abundance")+
#   ggtitle("National trajectory")+
#   scale_y_continuous(label=comma)+
#   coord_cartesian(ylim = c(0,75000000))
# # 
# # png(file = paste0(output_directory,"figures/Trajectory_18yr.png"), units = "in", width = 6, height = 5, res = 600)
# # N_national_18yr_plot
# # dev.off()
# 
# # ---------------------------------------------
# # Violin plot of posterior trend estimates
# # ---------------------------------------------
# 
# # Individual mcmc samples
# trend_df = bind_rows(stratum_trend_df,national_trend_df)
# trend_df$stratum_name = factor(trend_df$stratum_name, levels = c("National","West","East"))
# 
# decline_text = subset(trend_df_summary, end_year == 2018 & start_year == 2008)
# 
# # 10 year trend
# trend_violin_10yr = ggplot(subset(trend_df, end_year == 2018 & year_interval == 10), 
#                            aes(x = stratum_name, 
#                                y = trend, 
#                                fill = stratum_name)) +
#   geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
#   geom_violin(draw_quantiles = c(0.025,0.5,0.975), 
#               alpha = 0.7,
#               col = "gray35", size = 1) +
#   scale_fill_manual(values = c("gray80",strata_colours), guide = "none")+
#   geom_text(data = decline_text, aes(x = stratum_name, y = -25, label = 1-round(prob_increase,2)), fontface = "bold")+
#   xlab("Stratum")+
#   ylab("Trend\n(% change per year)")+
#   coord_cartesian(ylim = c(-25,25))+
#   ggtitle("10 year trend")
# print(trend_violin_10yr)
# 
# png(file = paste0(output_directory,"figures/Trend_Violin_10yr.png"), units = "in", width = 6, height = 5, res = 600)
# trend_violin_10yr
# dev.off()
# 
# # png(file = paste0(figure_directory,"/trend_violin_10yr.png"), units = "in", width = 6, height = 3, res = 1000)
# # trend_violin_10yr
# # dev.off()
# 
# # 18 year trend
# decline_text = subset(trend_df_summary, end_year == 2018 & start_year == 2000)
# trend_violin_18yr = ggplot(subset(trend_df, end_year == 2018 & year_interval == 18), 
#                            aes(x = stratum_name, 
#                                y = trend, 
#                                fill = stratum_name)) +
#   geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
#   geom_violin(draw_quantiles = c(0.025,0.5,0.975), 
#               alpha = 0.7,
#               col = "gray35", size = 1) +
#   scale_fill_manual(values = c("gray80",strata_colours), guide = "none")+
#   geom_text(data = decline_text, aes(x = stratum_name, y = -25, label = 1-round(prob_increase,2)), fontface = "bold")+
#   xlab("Stratum")+
#   ylab("Trend\n(% change per year)")+
#   coord_cartesian(ylim = c(-25,25))+
#   ggtitle("18 year trend")
# print(trend_violin_18yr)
# 
# png(file = paste0(output_directory,"figures/Trend_Violin_18yr.png"), units = "in", width = 6, height = 5, res = 600)
# trend_violin_18yr
# dev.off()
# 
# # 10 year trend with comparison to BBS
# load("./analysis/2_BBS_analysis/trend_10yr_samples.RData")
# bbs_trend_10yr_samples = trend_10yr_samples
# trend_df_BBS = data.frame(mcmc = 1:length(bbs_trend_10yr_samples),
#                           trend = bbs_trend_10yr_samples,
#                           stratum_name = "BBS")
# 
# trend_df_comparison = bind_rows(subset(trend_df,stratum_name == "National" & end_year == 2018 & year_interval == 10),trend_df_BBS)
# trend_df_comparison$stratum_name[which(trend_df_comparison$stratum_name == "National")] = "Migration"
# 
# # Label for plot
# text = trend_df_comparison %>%
#   group_by(stratum_name) %>%
#   summarize(prob_decline = round(mean(trend<0),2))
# 
# trend_comparison_plot = ggplot(trend_df_comparison,
#                                aes(x = stratum_name, 
#                                    y = trend)) +
#   geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
#   geom_violin(draw_quantiles = c(0.025,0.5,0.975), alpha = 0.7,
#               fill = "gray80",
#               col = "gray35", size = 1) +
#   geom_text(data = text, aes(x = stratum_name, y = -25, label = prob_decline), fontface = "bold")+
#   xlab("")+
#   ylab("Trend\n(% change per year)")+
#   coord_cartesian(ylim = c(-25,25))+
#   ggtitle(focal_season)
# trend_comparison_plot
# 
# png(file = paste0(output_directory,"figures/National_Trend_Comparison_BBS.png"), units = "in", width = 6, height = 5, res = 600)
# trend_comparison_plot 
# dev.off()
