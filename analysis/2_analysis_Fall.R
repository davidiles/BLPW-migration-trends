# *******************************************************************
# *******************************************************************

# Analysis of Postbreeding ("Fall") migration data

# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse',
              'jagsUI',
              'readxl',
              'reshape2',
              'sf') 
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

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/BLPW-migration-trends/analysis")

# ***************************************************************
# ***************************************************************
# PART 1: PREPARE DATA AND FIT MODEL
# ***************************************************************
# ***************************************************************

#------------------------------------------------
# Load strata shapefile; initially created using combine_strata.R
#------------------------------------------------

load("0_data/strata/strata_East_West.RData")
strata_sf <- strata$strata_sf
nstrata <- nrow(strata_sf)
strata_sf$stratum_number <- 1:nstrata

#------------------------------------------------
# Prepare to format data
#------------------------------------------------

focal_season <- "Fall"
start_year <- 1998
end_year <- 2018

# Relevant directories / set working directory
data_directory <- paste0("0_data/")
output_directory <- paste0("1_output/",focal_season,"/")
figure_directory <- paste0(output_directory,"/figures/")
table_directory <- paste0(output_directory,"/tables/")

#------------------------------------------------
# Format Canadian migration monitoring data
#------------------------------------------------

dat_can <- rbind(read.csv("0_data/migration_counts/CAN/2020_01_03/ACBO.BLPW.2018.csv") %>% add_column(., station = "ACBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/BPBO.BLPW.2018.csv") %>% add_column(., station = "BPBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/IPBO.BLPW.2018.csv") %>% add_column(., station = "IPBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/LMBO.BLPW.2018.csv") %>% add_column(., station = "LMBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/LPBO.BLPW.2018.csv") %>% add_column(., station = "LPBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/MGBO.BLPW.2018.csv") %>% add_column(., station = "MGBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/MNO.BLPW.2018.csv") %>% add_column(., station = "MNO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/PEPBO.BLPW.2018.csv") %>% add_column(., station = "PEPBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/PIBO.BLPW.2018.csv") %>% add_column(., station = "PIBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/RUTH.BLPW.2018.csv") %>% add_column(., station = "RUTH"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/TCBO.BLPW.2018.csv") %>% add_column(., station = "TCBO"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/TLBBS.BLPW.2018.csv") %>% add_column(., station = "TLBBS"),
                 read.csv("0_data/migration_counts/CAN/2020_01_03/TTPBRS.BLPW.2018.csv") %>% add_column(., station = "TTPBRS"))

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

#------------------------------------------------
# Read/format data from US stations
#------------------------------------------------

dat_usa <- rbind(readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - AIMS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - BIBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - BSBO fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - BSBO spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - FBBS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - FBBS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - KWRS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - MCCS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall"),
                 readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - MCCS spring.xlsx") %>% as.data.frame() %>% add_column(., season = "Spring"))

# datasets with different column names
tmp = readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - PARC fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
colnames(tmp) <- colnames(dat_usa)
dat_usa <- rbind(dat_usa, tmp)
rm(tmp)    

tmp = readxl::read_xlsx("0_data/migration_counts/USA/Cleaned BLPW - CFMS fall.xlsx") %>% as.data.frame() %>% add_column(., season = "Fall")
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
us_windows = read.csv("0_data/migration_counts/USA/US_station_windows.csv")
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

#------------------------------------------------
# Combine CAN and USA data into a single dataframe
#------------------------------------------------

dat_combined = dplyr::bind_rows(dat_combined_can, dat_combined_usa)
dat_combined$station[dat_combined$station == "MCCS"] <- "MBO" # Fix station label for manomet

# Limit to data collected in year range
dat_combined = subset(dat_combined, YearCollected >= start_year & YearCollected <= end_year)

# Restrict to season of interest
dat_combined <- subset(dat_combined, season == focal_season)
dat_combined$site <- paste0(dat_combined$station,dat_combined$area)

# Dummy variable for sub-sites within a station
dat_combined$dummy_site <- (dat_combined$area > 1) %>% as.numeric()

# For any sites without recorded net hours, fill in with mean for that station
# also fill all CMMN stations
mean_daily_effort <- dat_combined %>%
  group_by(station) %>%
  summarize(mean_net.hrs = mean(net.hrs, na.rm = TRUE))
overall_mean <- mean(mean_daily_effort$mean_net.hrs,na.rm = TRUE)
mean_daily_effort$mean_net.hrs[is.na(mean_daily_effort$mean_net.hrs)] <- overall_mean

for (s in unique(dat_combined$station)){
  rows_to_fill <- which(dat_combined$station == s & is.na(dat_combined$net.hrs))
  dat_combined$net.hrs[rows_to_fill] <- mean_daily_effort$mean_net.hrs[which(mean_daily_effort$station == s)]
}


#----------
# Summaries of data availability at each station
#----------

# Load coordinates of cmmn stations
station_coordinates <- read.csv("0_data/locations/station_locations.csv") 

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

# nyear
year_vec <- min(dat_combined$YearCollected):max(dat_combined$YearCollected)
nyear <- length(year_vec)

# nday
nday <- max(date_ranges$window_length)

# counts
count_df <- dat_combined[,c("site","YearCollected","doy","ObservationCount","net.hrs","dummy_site")] %>%
  rename(site_name = site, year_abs = YearCollected, day_number = doy, count = ObservationCount, net_hrs = net.hrs) %>%
  add_column(year_number = match(.$year_abs,year_vec),
             station = gsub('[[:digit:]]+', '', .$site_name)) %>%
  subset(!is.na(count))

count_df$station_number <- factor(count_df$station) %>% as.numeric()
station_names <- factor(count_df$station) %>% levels(.)

count_df <- count_df %>%
  group_by(station,station_number,year_abs,year_number,day_number)%>%
  summarize(count = sum(count),
            net_hrs = sum(net_hrs))

jags_data <- list(nstrata = nstrata,
                  nyear = max(count_df$year_number),
                  nstation = max(count_df$station_number),
                  nobs = nrow(count_df),
                  
                  day = count_df$day_number,
                  year = count_df$year_number,
                  station = count_df$station_number,
                  
                  count = count_df$count,
                  net_hrs = count_df$net_hrs)

# End section that processes count data
# ------------------------------------------------------



# ------------------------------------------------------
# PART 2: FORMAT BREEDING ORIGINS
# ------------------------------------------------------

# Breeding origin assignments
assignment_file = "0_data/Isotopes/isotope_assignments.xlsx"

# First 7 columns contain the relevant data
assignments <- read_xlsx(assignment_file)[,c("location","season","year","lat","lon","stratum_West","stratum_Central","stratum_East")] %>% 
  subset(season == focal_season & year %in% year_vec) %>%
  mutate(assigned_to_West = stratum_West + stratum_Central,
         assigned_to_East = stratum_East)

# Spatial datasets used to assign breeding origin data to each migration monitoring station
assignments_sf <- assignments %>% st_as_sf(coords = c("lon", "lat"),crs = 4269, agr = "constant", remove = FALSE)
station_data_summarized_sf <- station_data_summarized_sf %>% st_transform(st_crs(assignments_sf)) # convert station locations to same crs

# *****
# For each source of breeding origin data, assign it to all stations within 250 km to estimate their catchment
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
  
  inrange <- station_data_summarized_sf[which(dists < 250000),]
  
  assignment_year <- assignments_sf$year[i]
  
  for (station in inrange$station){
    if (sum(is.na(N_origin[,station,as.character(assignment_year)]))>0) N_origin[1:(dim(N_origin)[1]),station,as.character(assignment_year)] <- rep(0,length(N_origin[,station,as.character(assignment_year)]))
    N_origin["West",station,as.character(assignment_year)] <- N_origin["West",station,as.character(assignment_year)] + as.numeric(assignments[i,"assigned_to_West"])
    N_origin["East",station,as.character(assignment_year)] <- N_origin["East",station,as.character(assignment_year)] + as.numeric(assignments[i,"assigned_to_East"])
    
  }
  
}

# Sample sizes
N_station_sampled <- apply(N_origin,c(2,3), sum)
N_station_sampled[is.na(N_station_sampled)] <- 999  # Placeholder for sample size in years with no stable isotope information

# Append to jags data package
jags_data$N_origin <- N_origin
jags_data$N_station_sampled <- N_station_sampled

# End section processing breeding origins
# ------------------------------------------------------

# ******************************************************************
# PREPARE FOR ANALYSIS WITH JAGS
# ******************************************************************

#------------------------------------------------
# Fix migration parameters for several stations
#------------------------------------------------

# Fix certain migration rates to zero; stations with 0 for "rho_fix" will not have migration parameters estimated
rho_fix <- matrix(1,nrow = jags_data$nstrata, 
                  ncol = jags_data$nstation, 
                  dimnames = list(c("East","West"),station_names)) 

western_stations <- c("CFMS","TLBBS","MNO","LMBO")

# Western stations cannot receive eastern birds
rho_fix["East", western_stations] <- 0 

# Do not use isotope data in model for those stations
jags_data$N_origin[,western_stations,] <- NA
jags_data$rho_fix <- rho_fix

# ------------------------------------------------------
# Fit Model
# ------------------------------------------------------

parameters.to.save = c("slope",
                       "sigma_rho",
                       "sigma_X",
                       "sigma_stationday",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       "rho",
                       "M",
                       "T",
                       "expected_count",
                       "sim_count",
                       "X")


inits <- NULL
nsamp <- 2000
nb <- 10000
nt <- 50
ni <- nb + nsamp*nt

# out <- jags(data = jags_data,
#             model.file = "migration_model.jags",
#             parameters.to.save = parameters.to.save,
#             inits = inits,
#             n.chains = 3,
#             n.thin = nt,
#             n.iter = ni,
#             n.burnin = nb,
#             parallel = TRUE)
# 
# out$mcmc.info$elapsed.mins 

# ***************************************************************
# ***************************************************************
# PART 2: ASSESS MODEL CONVERGENCE AND EFFECTIVE SAMPLE SIZE
# ***************************************************************
# ***************************************************************

#-------------------------------------------------------------------
# Assess model convergence
#-------------------------------------------------------------------

latent_parameters <- c("slope",
                       "sigma_proc",
                       "sigma_rho",
                       "sigma_stationday",
                       "migration_phenology_sd",
                       "migration_phenology_mean",
                       "log_rho_mu")

latent_states <- names(out$Rhat)[which(!(names(out$Rhat) %in% latent_parameters | names(out$Rhat) %in% c("sim_count","X2_sim","X2_obs")))]

# Rhat
max(unlist(out$Rhat[latent_parameters]), na.rm = TRUE)
length(unlist(out$Rhat[latent_parameters]))
mean(unlist(out$Rhat[latent_parameters]) > 1.1, na.rm = TRUE)

max(unlist(out$Rhat[latent_states]), na.rm = TRUE)
length(unlist(out$Rhat[latent_states]))
mean(unlist(out$Rhat[latent_states]) > 1.1, na.rm = TRUE)
sum(unlist(out$Rhat[latent_states]) > 1.1, na.rm = TRUE)
unlist(out$Rhat[latent_states])[which(unlist(out$Rhat[latent_states]) > 1.1)]

# Effective sample sizes
n.eff <- unlist(out$n.eff[!(names(out$Rhat) %in% c("sim_count","X2_sim","X2_obs","daily_index"))])

# Parameters with fewer than 1000 samples
n.eff[n.eff > 1 & n.eff <= 1000] 
mean(n.eff[n.eff>1]<1000)

# Traceplot
MCMCvis::MCMCtrace(out, params = c("slope",
                                   "sigma_rho",
                                   "sigma_stationday",
                                   "migration_phenology_sd",
                                   "migration_phenology_mean",
                                   "X"),
                   Rhat = TRUE,n.eff = TRUE, ind = TRUE,
                   pdf = TRUE, filename = paste0("Traceplot_",focal_season,".pdf"), 
                   wd = paste0("1_output/",focal_season))

# ***************************************************************
# ***************************************************************
# PART 3: SUMMARIZE DATA AVAILABILITY FOR APPENDIX 1
# ***************************************************************
# ***************************************************************

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
  arrange(Lon) %>%
  mutate(Lat = round(Lat,1),
         Lon = round(Lon,1)) 

write.csv(station_summary, file = paste0("1_output/",focal_season,"/tables/station_summary.csv"),row.names = FALSE)

# ***************************************************************
# ***************************************************************
# PART 4: GOODNESS-OF-FIT
# ***************************************************************
# ***************************************************************

#-------------------------------------------------------------------
# Plot model predictions overlaid with observed counts (black x)
#-------------------------------------------------------------------

expected_vs_observed <- count_df %>%
  group_by(station,year_abs) %>%
  summarize(sum_obs = sum(count))

# Calculate uncertainty in expected count for each station-year combination
for (i in 1:nrow(expected_vs_observed)){
  
  # Which observations correspond to this year/station combination?
  j <- which(count_df$year_abs == expected_vs_observed$year_abs[i] & 
               count_df$station == expected_vs_observed$station[i])
  
  expected <- out$sims.list$expected_count[,j] %>% apply(.,1,sum) %>% quantile(c(0.025,0.5,0.975))
  expected_vs_observed$expected_q025[i] <-  expected[1]
  expected_vs_observed$expected_q50[i]  <-  expected[2]
  expected_vs_observed$expected_q975[i] <-  expected[3]
  expected_vs_observed$expected_mean[i] <-  out$sims.list$expected_count[,j] %>% apply(.,1,sum) %>% mean()
  
  
}

# Correlation between observed and expected counts
cor_obs_expected <- expected_vs_observed %>%
  group_by(station) %>%
  summarize(cor = round(cor(sum_obs,expected_mean),2),
            max = max(c(sum_obs,expected_q975)))

# Plot
cor_obs_expected$cor_label <- paste0("cor = ",cor_obs_expected$cor)

plot_obs_vs_expected <- ggplot(data = expected_vs_observed)+
  geom_errorbar(aes(x = year_abs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = year_abs, y = expected_mean, col = "Expected", shape = "Expected"))+
  geom_point(aes(x = year_abs, y = sum_obs, col = "Observed", shape = "Observed"))+
  geom_text(data = cor_obs_expected, aes(x = 1998, y = max*1.1, label = cor_label), hjust = 0, size = 3)+
  facet_wrap(station~., scales = "free_y")+
  scale_color_manual(values=c("gray75","black"), name = "", guide = "none")+
  scale_shape_manual(values=c(19,4), name = "", guide = "none")+
  xlab("Year")+
  ylab("Total Seasonal Count")+
  ggtitle("Observed vs Expected Seasonal Total Counts\n\nPostbreeding migration")
plot_obs_vs_expected

png(file = paste0(output_directory,"figures/Appendix_GOF_mig_counts_at_each_station.png"), units = "in", width = 8, height = 6, res = 600)
plot_obs_vs_expected
dev.off()

#-------------------------------------------------------------------
# Breeding origin assignments at each station
#-------------------------------------------------------------------

station_assignments <- jags_data$N_origin %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, Year = Var3, n = value)
station_assignments$Stratum = factor(station_assignments$Stratum, levels = c("West","East"))

station_fixed_rho <- jags_data$rho_fix %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, fix = value)
station_fixed_rho$Stratum = factor(station_fixed_rho$Stratum, levels = c("West","East"))
station_fixed_rho$Fixed <- NA
station_fixed_rho$Fixed[station_fixed_rho$fix == 0 & station_fixed_rho$Stratum == "East"] <- "Model will assume this station\nonly captures birds from West\nbased on its location"
station_fixed_rho$Fixed[station_fixed_rho$fix == 0 & station_fixed_rho$Stratum == "West"] <- "Model will assume this station\nonly captures birds from East\nbased on its location"

station_assignments$Station <- factor(station_assignments$Station,levels = station_summary$`Station code`)
station_assignment_plot <- ggplot(data = station_assignments, 
                                  aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  geom_text(data = station_fixed_rho, aes(x = 1998, y = max(station_assignments$n,na.rm = TRUE)*1.5, label = Fixed),
            hjust=0,vjust=1, fontface = "italic", size = 2)+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~.)+
  xlim(c(range(station_assignments$Year)))+
  ggtitle("Postbreeding Migration")+
  ylab("Number of bird samples analyzed")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_assignment_plot

png(file = paste0(output_directory,"figures/Appendix_GOF_breeding_origin_data.png"), units = "in", width = 8, height = 6, res = 600)
station_assignment_plot
dev.off()

#-------------------------------------------------------------------
# Posterior predictive checks
#-------------------------------------------------------------------

# For each station, in each year, calculate Bayesian p-value using chi-squared statistic 
# and plot simulated vs empirical discrepancy measures

Bayesian_pvals <- count_df %>%
  group_by(station,year_abs) %>%
  summarize(Bayesian_pval = NA)

# Calculate uncertainty in expected count for each station-year combination
for (i in 1:nrow(Bayesian_pvals)){
  
  # Which observations correspond to this year/station combination?
  j <- which(count_df$year_abs == Bayesian_pvals$year_abs[i] & 
               count_df$station == Bayesian_pvals$station[i])
  
  observed_count_sum  <- sum(count_df$count[j])
  expected_count_sum  <- out$sims.list$expected_count[,j] %>% apply(.,1,sum)
  simulated_count_sum <- out$sims.list$sim_count[,j] %>% apply(.,1,sum)
  
  X2_obs <- (observed_count_sum - expected_count_sum)^2/expected_count_sum
  X2_sim <- (simulated_count_sum - expected_count_sum)^2/expected_count_sum
  BPval <- mean(X2_obs>X2_sim)
  Bayesian_pvals$Bayesian_pval[i] <- BPval
  
}

BPval_plot <- ggplot(data = Bayesian_pvals,
                     aes(x = year_abs, y = Bayesian_pval)) +
  geom_point() +
  geom_hline(yintercept = c(0.2,0.8), col = "red", linewidth = 0.5,linetype = 2)+
  coord_cartesian(ylim=c(0,1))+
  ggtitle("Postbreeding Migration\n\nPosterior predictive checks\n\n(Proportion of observed datasets with\nlarger X2-statistic than simulated datasets)")+
  xlab("Year")+
  ylab("Bayesian p-value")+
  facet_wrap(station~.)+
  theme_bw()
BPval_plot

png(file = paste0("1_output/",focal_season,"/figures/Appendix_GOF_Posterior_Predictive_Checks.png"), units = "in", width = 8, height = 6, res = 600)
BPval_plot
dev.off()

# -------------------------------------------------------------------------
# Create a plot of daily observed counts for each station
# -------------------------------------------------------------------------

daily_count_plot <- ggplot(count_df)+
  geom_point(aes(x = day_number, y = count/net_hrs))+
  scale_y_continuous(trans="log10")+
  ylab("log(count / net hour)")+
  xlab("Day of the year")+
  facet_grid(station~year_abs, scales = "free_y")+
  ggtitle("Postbreeding migration")

png(file = paste0("1_output/",focal_season,"/figures/Appendix_daily_counts.png"), units = "in", width = 20, height = 12, res = 600)
daily_count_plot
dev.off()

# *******************************************************************
# *******************************************************************
# Calculate regional trends based on stratum-level indices (X)
# - Also use relative abundance estimates to calculate continental trends
# - Create figure illustrating trajectories
# - NOTE THAT RELATIVE ABUNDANCES IN EACH STRATUM NEED TO BE CALCULATED FIRST,
#   DONE IN SCRIPT 0
# *******************************************************************
# *******************************************************************

# ---------------------------------
# LOAD RELATIVE ABUNDANCE ESTIMATES
# ---------------------------------

year_vec <- seq(1998,2018)
relabund <- read.csv("1_output/Results_MainText/relabund_eBird_BAM.csv")

source_of_estimate <- "BAM"
relabund <- subset(relabund, Source == source_of_estimate)
relabund$Sum <- relabund$Sum/relabund$Sum[1]

# Year to which relative abundance will be "pinned"
if (source_of_estimate == "eBird") relabund_year_number <- which(year_vec == 2018)
if (source_of_estimate == "BAM") relabund_year_number <- which(year_vec == 2011)

# ---------------------------------
# Postbreeding migration
# ---------------------------------

# ~~~~~~~~~~
# Extract and rescale annual indices
# ~~~~~~~~~~

# Extract annual indices
X <- out$sims.list$X

# For each sample from posterior, re-scale estimates based on relative abundance (for continental trend estimate)
stratum_indices_df <- X %>% reshape2::melt() %>% 
  rename(sample = Var1, stratum_number = Var2, year_number = Var3, indices = value) %>%
  mutate(Stratum = c("East","West")[stratum_number],
         Year = year_vec[year_number],
         indices_rescaled = NA)

annual_means <- stratum_indices_df %>%
  group_by(Year,year_number,Stratum,stratum_number) %>%
  summarize(index_median = median(indices))

for (mcmc in 1:dim(X)[1]){
  
  # Estimates of log-linear slope
  for (stratum_number in 1:jags_data$nstrata){
    
    indices <- X[mcmc,stratum_number,]
    
    # Mean index in year of relative abundance estimate
    median_in_relabund_year <- annual_means$index_median[annual_means$year_number == relabund_year_number & annual_means$stratum_number == stratum_number]
    indices_rescaled <- indices/median_in_relabund_year * relabund$Sum[stratum_number]
    
    # Fill in dataframe
    row <- which(stratum_indices_df$sample == mcmc & stratum_indices_df$stratum_number == stratum_number)
    stratum_indices_df$indices_rescaled[row] <- indices_rescaled
    
  }
  
  print(mcmc)
  
}

# ~~~~~~~~~~
# Summarize stratum-level indices
# ~~~~~~~~~~

# Annual indices
stratum_indices_summarized <- stratum_indices_df %>%
  group_by(Stratum,Year) %>%
  summarize(
    
    # On original scale
    indices_q025 = quantile(indices,0.025),
    indices_q500 = quantile(indices,0.500),
    indices_q975 = quantile(indices,0.975),
    
    # Rescaled
    indices_rescaled_q025 = quantile(indices_rescaled,0.025),
    indices_rescaled_q500 = quantile(indices_rescaled,0.500),
    indices_rescaled_q975 = quantile(indices_rescaled,0.975)
  )

# ~~~~~~~~~~
# Continental estimates based on sum of rescaled X estimates
# ~~~~~~~~~~

continental_indices_df <- stratum_indices_df %>%
  group_by(sample,Year) %>%
  summarize(indices_rescaled = sum(indices_rescaled)) %>%
  mutate(Stratum = "Continental")

# Summarize
continental_indices_summarized <- continental_indices_df %>%
  group_by(Year,Stratum) %>%
  summarize(
    
    indices_rescaled_q025 = quantile(indices_rescaled,0.025),
    indices_rescaled_q500 = quantile(indices_rescaled,0.500),
    indices_rescaled_q975 = quantile(indices_rescaled,0.975)
  )

# ~~~~~~~~~~
# Combine stratum-level and continental estimates into a single dataframe (for plotting)
# ~~~~~~~~~~

indices_df <- bind_rows(stratum_indices_df,continental_indices_df)
indices_summarized <- bind_rows(stratum_indices_summarized,continental_indices_summarized)

# ~~~~~~~~~~
# Calculate trends
# ~~~~~~~~~~

yend = 2018 # End year
y0 = 1998   # Start year

trend_df <- indices_df %>%
  group_by(Stratum,sample) %>%
  mutate(trend = 100*((indices_rescaled[Year == yend]/indices_rescaled[Year == y0])^(1/(yend-y0))-1),
         perc_change_1998 = 100*(indices_rescaled[Year == 2018] - indices_rescaled[Year == 1998])/indices_rescaled[Year == 1998],
         perc_change_2008 = 100*(indices_rescaled[Year == 2018] - indices_rescaled[Year == 2008])/indices_rescaled[Year == 2008]
  )

trend_summarized <- trend_df %>%
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

Fall_results <- list(indices_df = indices_df,
                       indices_summarized = indices_summarized,
                       trend_df = trend_df,
                       trend_summarized = trend_summarized)

Fall_results$trend_summarized %>% as.data.frame()

saveRDS(Fall_results,paste0(output_directory,"results_Fall.rds"))

# Save/load entire workspace
#save.image(paste0(output_directory,"wksp_Fall.RData"))
