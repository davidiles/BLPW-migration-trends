# *******************************************************************
# Fits a standard, stratified analysis to spring migration data
# *******************************************************************

# ***************************************************************
# ***************************************************************
# PART 1: PREPARE DATA AND FIT MODEL
# ***************************************************************
# ***************************************************************

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

# Location of script on machine
proj_loc <- "D:/Working_Files/1_Projects/Landbirds/BLPW-migration-trends-2022-2strata/"
setwd(proj_loc)

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
start_year <- 2000
end_year <- 2018

# Relevant directories / set working directory
data_directory <- paste0(proj_loc,"analysis/0_data/")
output_directory <- paste0(proj_loc,"analysis/2_output/",focal_season,"/")
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
  station_coordinates <- read.csv("analysis/0_data/locations/station_locations.csv")
  
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
  
  # # For each station, apply breeding origin assignments within buffer
  # for (i in 1:nrow(station_data_summarized_sf)){
  #   
  #   # Calculate 
  #   station_name = station_data_summarized_sf$station[i]
  #   dists <- st_distance(station_data_summarized_sf[i,],assignments_sf) %>% as.numeric()
  #   within_buffer <- which(dists <= 250000)
  #   
  #   # If no assignments are available, move to next station
  #   if (length(within_buffer) < 1) next
  #   
  #   breeding_origins_at_station <- assignments[within_buffer,] %>%
  #     group_by(year) %>%
  #     summarize(assigned_to_West = sum(assigned_to_West),
  #               assigned_to_East = sum(assigned_to_East))
  #   
  #   # Fill in data
  #   for (j in 1:nrow(breeding_origins_at_station)){
  #     N_origin["East",station_name,as.character(breeding_origins_at_station$year[j])] = breeding_origins_at_station$assigned_to_East[j]
  #     N_origin["West",station_name,as.character(breeding_origins_at_station$year[j])] = breeding_origins_at_station$assigned_to_West[j]
  #     
  #   }
  #   
  # }
  
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

# Fix certain migration rates to zero; stations with 0 for "rho.fix" will not have migration parameters estimated
rho.fix <- matrix(1,nrow = jags_data$nstrata, 
                  ncol = jags_data$nstation, 
                  dimnames = list(c("East","West"),station_names)) 

# Western stations cannot receive eastern birds
rho.fix["East", c("CFMS","TLBBS","MNO","LMBO")] <- 0 

#------------------------------------------------
# Add several extra quantities that are needed for JAGS script
#------------------------------------------------

jags_data$rho.fix <- rho.fix
jags_data$pi <- pi
jags_data$logN0 <- rep(0,jags_data$nstrata)

# "Station" at each site
station_site <- count_df %>% group_by(station,site_number) %>% summarize(station_number = mean(station_number))

jags_data$station_site <- station_site$station_number

# ------------------
# FIT THE MODEL IN JAGS
# ------------------

parameters.to.save = c("trend",
                       "sigma_process",
                       "sigma_stationyear",
                       "sigma_stationday",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       
                       "log_rho_mu",
                       "rho",
                       "T",
                       "station_composition",
                       "expected_count",
                       "sim_count",
                       
                       "X",
                       "stratum_annual_growth",
                       "percent_change",
                       "X2_sim",
                       "X2_obs"
)


inits <- NULL #function()list(sigma_stationyear = 0.3)

# ------------------
# The jags script to fit the model
# ------------------

sink("./analysis/1_scripts/migration_model_0.jags")
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
    
    trend[j] ~ dnorm(0,0.25)
    
    # Work forwards from baseline year
    for (y in 1:nyear){

      logX[j,y] ~ dnorm( trend[j] * (y-1) , tau_proc)
      X[j,y] <- exp(logX[j,y])

    }
  
  } # j
  
  #---------------------------------------------
  # Model for breeding origins of migrants arriving at each station
  #---------------------------------------------
  
  sigma_stationyear ~ dunif(0,2)
  tau_stationyear <- pow(sigma_stationyear,-2)

  for (s in 1:nstation){

    for (j in 1:nstrata){
      
      rho_mu[j,s] ~ dunif(0,10)
      log_rho_mu[j,s] <- log(rho_mu[j,s])

      for (y in 1:nyear){

	      rho[j,s,y] ~ dlnorm(log_rho_mu[j,s] - 1/(2*tau_stationyear),tau_stationyear)
        M[j,s,y] <- X[j,y] * rho[j,s,y] * rho.fix[j,s]

      }
    }
  }
  
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

    } # y
  } # s 
  
  #---------------------------------------------
  # Within-season model for migration counts
  #---------------------------------------------

  migration_phenology_sd ~ dunif(0,20)
  migration_phenology_tau <- pow(migration_phenology_sd,-2)
  
  for (s in 1:nstation){
    migration_phenology_mean[s] ~ dunif(1,360)
    
  }
  
  # Site-level fixed effects (for LPBO that contains multiple sub-stations)
  for (k in 1:nsite){

    site_effect[k] ~ dnorm(0,0.25)
    
    # Daily overdispersion in counts at each site (e.g., due to daily weather)
    sigma_stationday[k] ~ dunif(0,2)
    tau_stationday[k] <- pow(sigma_stationday[k],-2)
  }
  
  for (i in 1:nobs){
    
    mu[i] <- log(f[i]) + log(T[station[i],year[i]]) + log(net_hrs[i]) + site_effect[site[i]]*dummy_site[i]
    expected_count[i] <- exp(mu[i])

    # Likelihood for counts
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[station[i]], migration_phenology_tau))
    
    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday[site[i]]), tau_stationday[site[i]]) 
    count[i] ~ dpois(exp(log_lambda[i]))
    
    # *********************************************************
    # Simulate counts for posterior predictive checking
    # *********************************************************

    sim_log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday[site[i]]), tau_stationday[site[i]]) 
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
        station_composition[j,s,y] <- T[s,y] * p[j,s,y]
      }
    }
  }
  
  #---------------------------------------------
  # Percent change from first to last year
  #---------------------------------------------
  
  for (j in 1:nstrata){
    percent_change[j] <- 100 * (X[j,nyear] - X[j,1]) / X[j,1]
  }
  
  #---------------------------------------------
  # Annual percent change
  #---------------------------------------------
  
  for (y in 2:nyear){
    for (j in 1:nstrata){
      stratum_annual_growth[j,y-1] <- 100*(X[j,y]-X[j,y-1])/X[j,y-1]
    }
  }
  
}
    ",fill = TRUE)
sink()

# ------------------
# Fit Model
# ------------------

nsamp <- 5000
nb <- 20000
nt <- 50
ni <- nb + nsamp*nt

# Need to rerun with time-varying model
out <- jags(data = jags_data,
            model.file = "analysis/1_scripts/migration_model_0.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 3,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = TRUE)

out$mcmc.info$elapsed.mins


# # ------------------
# # Save/load workspace
# # ------------------
# 
save.image(paste0(output_directory,"analysis_",focal_season,"_0.RData"))

# Load fitted model
load(paste0(output_directory,"analysis_",focal_season,"_0.RData"))


# ***************************************************************
# ***************************************************************
# PART 2: GOODNESS-OF-FIT
# ***************************************************************
# ***************************************************************

#-------------------------------------------------------------------
# Assess model convergence
#-------------------------------------------------------------------

MCMCtrace(out,c("trend",
                "sigma_process",
                "sigma_migration",
                "sigma_daily",
                "migration_phenology_sd",
                "migration_phenology_mean",
                "log_rho_mu"), Rhat = TRUE,
          pdf = FALSE)

latent_parameters <- c("trend",
                       "sigma_process",
                       "sigma_migration",
                       "sigma_daily",
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
n.eff <- unlist(out$n.eff)
n.eff[n.eff > 1 & n.eff <= 1000] # Parameters with fewer than 1000 samples



# ***************************************************************
# ***************************************************************
# PART 3: INTERPRET RESULTS
# ***************************************************************
# ***************************************************************

#-------------------------------------------------------------------
# Breeding origin assignments at each station
#-------------------------------------------------------------------

station_assignments <- jags_data$N_origin %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, Year = Var3, n = value)
station_assignments$Stratum = factor(station_assignments$Stratum, levels = c("West","East"))

station_assignment_plot <- ggplot(station_assignments, 
                                  aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~., scales = "free")+
  xlim(c(range(station_assignments$Year)))+
  ggtitle(focal_season)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_assignment_plot

# png(file = paste0(output_directory,"figures/Results_Origins_Observed.png"), units = "in", width = 8, height = 6, res = 600)
# station_assignment_plot
# dev.off()

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
  ggtitle(focal_season)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_composition_plot

# png(file = paste0(output_directory,"figures/Results_Origins_Estimated.png.png"), units = "in", width = 8, height = 6, res = 600)
# station_composition_plot
# dev.off()

# -----------------------------------------
# Map of study sites and which ones have isotopes available
# -----------------------------------------

BCR_poly <- strata$BCR_poly
xlim <- c(-162, -56)
ylim <- c(35, 70)

# Add indicator for assignment information available
station_data_summarized_sf$assignments <- "No"
station_data_summarized_sf$assignments[station_data_summarized_sf$station %in% na.omit(station_assignments)$Station] <- "Yes"

origins_known <- colSums(rho.fix == 0)>0
origins_known <- origins_known[which(origins_known)]
station_data_summarized_sf$assignments[which(station_data_summarized_sf$station %in% names(origins_known))] <- "Yes"

strata_sf$Strata <- factor(strata_sf$Strata, levels = c("West","East"))
strata_fig <- ggplot(data = strata_sf) +
  geom_sf(data = BCR_poly, fill = "gray96", col = "gray85")+
  geom_sf(data = strata_sf, aes(fill = Strata), alpha = 0.7)+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_fill_manual(values = strata_colours, name = "Stratum") +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed", size = 0.5))+
  geom_sf(data = station_data_summarized_sf ,size = 1, alpha = 0.5, shape = 7,col = "black")+
  geom_sf_label_repel(data = station_data_summarized_sf,
                      aes(label = station,
                          col = assignments),
                      max.overlaps = 50,
                      size = 1.5, alpha = 0.8)+
  scale_color_manual(values = c("gray50","black"), name = "Breeding Origins Estimated") +
  ggtitle(focal_season)+
  coord_sf(xlim = xlim, ylim = ylim)
#strata_fig

# png(file = paste0(output_directory,"figures/Results_Origins_Assignment_Map.png"), units = "in", width = 8, height = 5, res = 600)
# strata_fig 
# dev.off()

#-------------------------------------------------------------------
# Station-level indices each year
#-------------------------------------------------------------------

station_indices_full <- out$sims.list$T %>% 
  reshape2::melt() %>%
  rename(samp = Var1, 
         station_number = Var2, 
         year_number = Var3, 
         index = value)

station_indices_full$Station = station_names[station_indices_full$station_number]
station_indices_full$Year = year_vec[station_indices_full$year_number]

station_indices <- station_indices_full %>%
  group_by(Station,Year) %>%
  summarize(index_mean = mean(index),
            index_q50 = quantile(index,0.5),
            index_q025 = quantile(index,0.025),
            index_q975 = quantile(index,0.975))

station_indices <- full_join(station_indices, station_years_with_counts)

station_index_plot <- ggplot(data = na.omit(station_indices), 
                             aes(x = Year, 
                                 y = index_mean, 
                                 ymin = index_q025, 
                                 ymax = index_q975)) +
  geom_point()+
  geom_errorbar(width=0)+
  facet_wrap(Station~., scales = "free_y")

station_index_plot

#write.csv(station_indices,paste0(output_directory,"tables/station_indices_Spring.csv"),row.names = FALSE)

# png(file = paste0(output_directory,"figures/Results_Station_Indices.png"), units = "in", width = 8, height = 5, res = 600)
# station_index_plot
# dev.off()

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

regional_trajectory_plot <- ggplot(X_summary,aes(x = Year, y = X_mean, ymin = X_lcl, ymax = X_ucl, fill = Stratum, col = Stratum))+
  geom_ribbon(alpha = 0.2)+
  geom_line()+
  facet_grid(Stratum~.)+
  ylab("Regional index")+
  xlab("Year")+
  ggtitle("Regional trajectory")

regional_trajectory_plot
# 
# png(file = paste0(output_directory,"figures/Results_Regional_Trajectory.png"), units = "in", width = 8, height = 5, res = 600)
# regional_trajectory_plot
# dev.off()



# -------------------------------------------------------------------
# Calculate regional and national abundance in each year using BAM's bootstrap replicates
# -------------------------------------------------------------------

# Specify year represented by BAM density raster
bam_year <- 2011
bam_year_index <- which(year_vec == bam_year)
bam_bootstrap_files <- list.files("analysis/0_data/bam_density_raster/bam_bootstrap/")

N_strata_df_total <- data.frame()
for (i in 1:length(bam_bootstrap_files)){
  
  # load raster
  bam <- raster(paste0("analysis/0_data/bam_density_raster/bam_bootstrap/",bam_bootstrap_files[i]))
  
  # Crop to stratum boundaries
  strata_sp <- strata_sf %>% st_transform(crs = projection(bam)) %>% as('Spatial')
  bam <- crop(bam,strata_sp) %>% mask(strata_sp)
  
  # Extract cumulative abundance in each stratum
  stratum_N <- rep(NA,jags_data$nstrata)
  for (r in 1:nstrata){
    bam_stratum <- crop(bam,strata_sp[r,]) %>% mask(strata_sp[r,]) 
    stratum_N[r] <- sum(values(bam_stratum),na.rm = TRUE) * 100 # Multiply by 100 because each pixel is 100 ha (and bam is density in males/ha)
  }
  
  # Empty array to store total abundances
  logN_array <- array(NA,c(out$mcmc.info$n.samples,jags_data$nyear,jags_data$nstrata))
  
  # Loop through strata
  for (j in 1:nstrata){
    
    # BAM estimate for this stratum in the year of the BAM estimate
    logN_array[,bam_year_index,j] <- log(stratum_N[j])
    
    # Work forwards from the year of the BAM estimate
    if (bam_year_index < jags_data$nyear){
      for (y in bam_year_index:(jags_data$nyear-1)){
        
        # Posterior samples of annual growth rate for this year
        annual_growth <- log(out$sims.list$X[,j,y+1]/out$sims.list$X[,j,y])
        logN_array[,y+1,j] <- logN_array[,y,j] + annual_growth
      }
    }
    
    # Work backwards from the year of the BAM estimate
    if (bam_year_index > 1 ){
      for (y in bam_year_index:2){
        
        # Posterior samples of annual growth rate for this year
        annual_growth <- log(out$sims.list$X[,j,y]/out$sims.list$X[,j,y-1])
        logN_array[,y-1,j] <- logN_array[,y,j] - annual_growth
      }
    }
  }
  
  N_strata_df <- melt(logN_array) %>% 
    rename(mcmc = Var1, year = Var2, stratum = Var3, logN = value) %>%
    add_column(bootstrap_rep = bam_bootstrap_files[i])
  N_strata_df$N <- exp(N_strata_df$logN)
  
  N_strata_df_total <- rbind(N_strata_df_total,N_strata_df)
  
  print(i)
  
}

N_strata_df_total$year_abs <- year_vec[N_strata_df_total$year]

# Summarize estimates
N_strata_estimates <- N_strata_df_total %>% 
  group_by(stratum,year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975)
  )

N_strata_estimates$year_abs <- year_vec[N_strata_estimates$year]
N_strata_estimates$stratum_name <- dimnames(jags_data$N_origin)[[1]][N_strata_estimates$stratum]
# save(N_strata_estimates, file = paste0(output_directory,"/R_objects/",focal_season,"_N_strata_estimates.RData"))

# National abundance
N_national_estimates_boot <- N_strata_df_total %>% 
  group_by(year,mcmc,bootstrap_rep) %>%
  
  # Sum across all strata
  summarize(N = sum(N))

N_national_estimates <- N_national_estimates_boot %>%
  group_by(year) %>%                    
  
  # Summarize across BAM bootstrap replicates
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))
N_national_estimates$year_abs <- year_vec[N_national_estimates$year]
# save(N_national_estimates, file = paste0(output_directory,"/R_objects/",focal_season,"_N_national_estimates.RData"))

#-------------------------------------------------------------------
# Calculate trends and percent change across specific intervals
#-------------------------------------------------------------------
national_trend_df <- data.frame()

# 10-year trends and study-wide trends
for (year_interval in c(10,length(unique(N_strata_df_total$year))-1)){
  for (year in 1:(max(N_national_estimates_boot$year)-1)){
    for (bootrep in unique(N_national_estimates_boot$bootstrap_rep)){
      
      start_year_pchange <- year
      end_year_pchange <- year+year_interval
      
      if (!(end_year_pchange %in% N_national_estimates_boot$year)) next
      
      N1 <- subset(N_national_estimates_boot,year == start_year_pchange & bootstrap_rep == bootrep)$N
      N2 <- subset(N_national_estimates_boot, year == end_year_pchange & bootstrap_rep == bootrep)$N
      
      percent_change = (N2-N1)/N1 * 100
      trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
      
      tmp <- data.frame(mcmc = 1:length(trend),
                        end_year_pchange = end_year_pchange,
                        year_interval = year_interval,
                        percent_change = percent_change,
                        trend = trend)
      
      national_trend_df <- rbind(national_trend_df, tmp)
    } # bootrep
    
    print(year)
  } # Start year
} # Year interval

national_trend_df$stratum_name <- "National"
national_trend_df$end_year <- year_vec[national_trend_df$end_year_pchange]
national_trend_df$start_year <- national_trend_df$end_year - national_trend_df$year_interval


# Summarize credible intervals across bootstrap replicates
national_trend_df_summary <- national_trend_df %>%
  group_by(stratum_name,start_year,end_year, year_interval) %>%
  summarize(percent_change_mean = mean(percent_change),
            percent_change_q025 = quantile(percent_change,0.025),
            percent_change_q975 = quantile(percent_change,0.975),
            
            trend_mean = mean(trend),
            trend_q025 = quantile(trend,0.025),
            trend_q975 = quantile(trend,0.975),
            
            prob_increase = mean(percent_change > 0),
            prob_30_decline = mean(percent_change <= -30),
            prob_50_decline = mean(percent_change <= -50))


# Stratum-level analysis
stratum_trend_df <- data.frame()
for (year_interval in c(10,length(unique(N_strata_df_total$year))-1)){
  for (year in 1:(max(N_strata_df_total$year)-1)){
    for (current_stratum in 1:nstrata){
      
      start_year_pchange <- year
      end_year_pchange <- year+year_interval
      
      if (!(end_year_pchange %in% N_strata_df_total$year)) next
      
      N1 <- subset(N_strata_df_total,year == start_year_pchange & bootstrap_rep == N_national_estimates_boot$bootstrap_rep[1] & stratum == current_stratum)$N
      N2 <- subset(N_strata_df_total, year == end_year_pchange & bootstrap_rep == N_national_estimates_boot$bootstrap_rep[1] & stratum == current_stratum)$N
      
      percent_change = (N2-N1)/N1 * 100
      trend <- 100*((N2/N1)^(1/year_interval)-1) # Equation from Smith et al. 2014
      
      tmp <- data.frame(mcmc = 1:length(r),
                        stratum = current_stratum,
                        end_year_pchange = end_year_pchange,
                        year_interval = year_interval,
                        percent_change = percent_change,
                        trend = trend)
      
      stratum_trend_df <- rbind(stratum_trend_df, tmp)
      
    } # stratum
    
    print(year)
  } # Start year
} # Year interval

stratum_trend_df$stratum_name <- dimnames(jags_data$N_origin)[[1]][stratum_trend_df$stratum]
stratum_trend_df$end_year <- year_vec[stratum_trend_df$end_year_pchange]
stratum_trend_df$start_year <- stratum_trend_df$end_year - stratum_trend_df$year_interval

# Summarize credible intervals across bootstrap replicates
stratum_trend_df_summary <- stratum_trend_df %>%
  group_by(stratum_name,start_year,end_year, year_interval) %>%
  summarize(percent_change_mean = mean(percent_change),
            percent_change_q025 = quantile(percent_change,0.025),
            percent_change_q975 = quantile(percent_change,0.975),
            
            trend_mean = mean(trend),
            trend_q025 = quantile(trend,0.025),
            trend_q975 = quantile(trend,0.975),
            
            prob_increase = mean(percent_change > 0),
            prob_30_decline = mean(percent_change <= -30),
            prob_50_decline = mean(percent_change <= -50))


trend_df_summary <- bind_rows(national_trend_df_summary,stratum_trend_df_summary) %>%
  add_column(season = focal_season) %>%
  relocate(stratum_name,season)

# Rearrange columns
trend_df_summary <- trend_df_summary %>% 
  ungroup() %>%
  dplyr::select(stratum_name,season,start_year, end_year,
                trend_mean,trend_q025,trend_q975,
                percent_change_mean,percent_change_q025,percent_change_q975,
                prob_increase,prob_30_decline,prob_50_decline)

# SAVE RELEVANT OUTPUT / SUMMARIES IN TABLE FORMAT
# write.csv(trend_df_summary,paste0(output_directory,"tables/trend_and_change_estimates.csv"),row.names = FALSE)

# *****************************************************
# *****************************************************
# PART 3: GENERATE SOME PLOTS OF DATA AND RESULT SUMMARIES
# *****************************************************
# *****************************************************

# Trajectory over last 10 years
N_strata_estimates$stratum_name = factor(N_strata_estimates$stratum_name, levels = c("West","East"))

N_strata_10yr_plot <- ggplot(subset(N_strata_estimates,year_abs %in% ((end_year-10):end_year))) +
  geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975, fill = stratum_name), alpha = 0.3)+
  geom_line(aes(x = year_abs, y = N_q500, col = stratum_name))+
  scale_fill_manual(values = strata_colours, name = "Stratum")+
  scale_color_manual(values = strata_colours, name = "Stratum")+
  xlab("Year")+
  ylab("Abundance")+
  facet_grid(.~stratum_name)+
  scale_y_continuous(label=comma)+
  scale_x_continuous(breaks = seq(2008,2018,length.out = 3))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position=c(.93,.75))+
  coord_cartesian(ylim = c(0,75000000))
N_strata_10yr_plot

# png(file = paste0(output_directory,"figures/Trajectory_Strata_10yr.png"), units = "in", width = 6, height = 5, res = 600)
# N_strata_10yr_plot
# dev.off()

# Trajectory over entire time series
N_strata_18yr_plot <- ggplot(N_strata_estimates) +
  geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975, fill = stratum_name), alpha = 0.3)+
  geom_line(aes(x = year_abs, y = N_q500, col = stratum_name))+
  scale_fill_manual(values = strata_colours, name = "Stratum", guide = FALSE)+
  scale_color_manual(values = strata_colours, name = "Stratum", guide = FALSE)+
  xlab("Year")+
  ylab("Abundance")+
  facet_grid(stratum_name~.)+
  scale_y_continuous(label=comma)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position=c(.93,.75))+
  coord_cartesian(ylim = c(0,75000000))

# png(file = paste0(output_directory,"figures/Trajectory_Strata_18yr.png"), units = "in", width = 6, height = 5, res = 600)
# N_strata_18yr_plot
# dev.off()

N_national_10yr_plot <- ggplot(subset(N_national_estimates,year_abs %in% ((end_year-10):end_year))) +
  geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975), alpha = 0.3)+
  geom_line(aes(x = year_abs, y = N_q500))+
  xlab("Year")+
  ylab("National Abundance")+
  ggtitle("National trajectory")+
  scale_x_continuous(breaks = seq(2008,2018,length.out = 3))+
  scale_y_continuous(label=comma)+
  coord_cartesian(ylim = c(0,75000000))

# png(file = paste0(output_directory,"figures/Trajectory_10yr.png"), units = "in", width = 6, height = 5, res = 600)
# N_national_10yr_plot
# dev.off()

N_national_18yr_plot <- ggplot(N_national_estimates) +
  geom_ribbon(aes(x = year_abs, ymin = N_q025, ymax = N_q975), alpha = 0.3)+
  geom_line(aes(x = year_abs, y = N_q500))+
  xlab("Year")+
  ylab("National Abundance")+
  ggtitle("National trajectory")+
  scale_y_continuous(label=comma)+
  coord_cartesian(ylim = c(0,75000000))
# 
# png(file = paste0(output_directory,"figures/Trajectory_18yr.png"), units = "in", width = 6, height = 5, res = 600)
# N_national_18yr_plot
# dev.off()

# ---------------------------------------------
# Violin plot of posterior trend estimates
# ---------------------------------------------

# Individual mcmc samples
trend_df = bind_rows(stratum_trend_df,national_trend_df)
trend_df$stratum_name = factor(trend_df$stratum_name, levels = c("National","West","East"))

decline_text = subset(trend_df_summary, end_year == 2018 & start_year == 2008)

# 10 year trend
trend_violin_10yr = ggplot(subset(trend_df, end_year == 2018 & year_interval == 10), 
                           aes(x = stratum_name, 
                               y = trend, 
                               fill = stratum_name)) +
  geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975), 
              alpha = 0.7,
              col = "gray35", size = 1) +
  scale_fill_manual(values = c("gray80",strata_colours), guide = "none")+
  geom_text(data = decline_text, aes(x = stratum_name, y = -25, label = 1-round(prob_increase,2)), fontface = "bold")+
  xlab("Stratum")+
  ylab("Trend\n(% change per year)")+
  coord_cartesian(ylim = c(-25,25))+
  ggtitle("10 year trend")
print(trend_violin_10yr)

png(file = paste0(output_directory,"figures/Trend_Violin_10yr.png"), units = "in", width = 6, height = 5, res = 600)
trend_violin_10yr
dev.off()

# png(file = paste0(figure_directory,"/trend_violin_10yr.png"), units = "in", width = 6, height = 3, res = 1000)
# trend_violin_10yr
# dev.off()

# 18 year trend
decline_text = subset(trend_df_summary, end_year == 2018 & start_year == 2000)
trend_violin_18yr = ggplot(subset(trend_df, end_year == 2018 & year_interval == 18), 
                           aes(x = stratum_name, 
                               y = trend, 
                               fill = stratum_name)) +
  geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975), 
              alpha = 0.7,
              col = "gray35", size = 1) +
  scale_fill_manual(values = c("gray80",strata_colours), guide = "none")+
  geom_text(data = decline_text, aes(x = stratum_name, y = -25, label = 1-round(prob_increase,2)), fontface = "bold")+
  xlab("Stratum")+
  ylab("Trend\n(% change per year)")+
  coord_cartesian(ylim = c(-25,25))+
  ggtitle("18 year trend")
print(trend_violin_18yr)

png(file = paste0(output_directory,"figures/Trend_Violin_18yr.png"), units = "in", width = 6, height = 5, res = 600)
trend_violin_18yr
dev.off()

# 10 year trend with comparison to BBS
load("./analysis/2_BBS_analysis/trend_10yr_samples.RData")
bbs_trend_10yr_samples = trend_10yr_samples
trend_df_BBS = data.frame(mcmc = 1:length(bbs_trend_10yr_samples),
                          trend = bbs_trend_10yr_samples,
                          stratum_name = "BBS")

trend_df_comparison = bind_rows(subset(trend_df,stratum_name == "National" & end_year == 2018 & year_interval == 10),trend_df_BBS)
trend_df_comparison$stratum_name[which(trend_df_comparison$stratum_name == "National")] = "Migration"

# Label for plot
text = trend_df_comparison %>%
  group_by(stratum_name) %>%
  summarize(prob_decline = round(mean(trend<0),2))

trend_comparison_plot = ggplot(trend_df_comparison,
                               aes(x = stratum_name, 
                                   y = trend)) +
  geom_hline(yintercept = 0, col = "gray80", size = 1.5)+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975), alpha = 0.7,
              fill = "gray80",
              col = "gray35", size = 1) +
  geom_text(data = text, aes(x = stratum_name, y = -25, label = prob_decline), fontface = "bold")+
  xlab("")+
  ylab("Trend\n(% change per year)")+
  coord_cartesian(ylim = c(-25,25))+
  ggtitle(focal_season)
trend_comparison_plot

png(file = paste0(output_directory,"figures/National_Trend_Comparison_BBS.png"), units = "in", width = 6, height = 5, res = 600)
trend_comparison_plot 
dev.off()
