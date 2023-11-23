# *******************************************************************
# *******************************************************************

# Analysis of post-breeding (fall) migration data

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

western_stations <- c("CFMS","TLBBS","MNO","LMBO")

# Western stations cannot receive eastern birds
rho_fix["East", western_stations] <- 0 

# Do not use isotope data in model for those stations
jags_data$N_origin[,western_stations,] <- NA

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

# Prepared by previous script

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
nsamp <- 1000
nb <- 10000
nt <- 100
ni <- nb + nsamp*nt

# Need to rerun with time-varying model
out <- jags(data = jags_data,
            model.file = "analysis/1_scripts/migration_model.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 3,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = TRUE)

out$mcmc.info$elapsed.mins # 408 mins
save.image(paste0(output_directory,"analysis_",focal_season,".RData"))

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

rho_fix <- jags_data$rho_fix
rho_fix[rho_fix == 0] <- NA
hist(out$q50$rho*rho_fix, breaks = seq(0,10,0.1), xlim = c(0,2.5))
hist(out$mean$rho*rho_fix, breaks = seq(0,10,0.1), xlim = c(0,2.5))

mean(log(out$q50$rho*rho_fix),na.rm = TRUE)
median(log(out$q50$rho*rho_fix),na.rm = TRUE)
sd(log(out$q50$rho*rho_fix),na.rm = TRUE)

a = rlnorm(100000,-3,1.5)
#-------------------------------------------------------------------
# Assess model convergence
#-------------------------------------------------------------------

# MCMCtrace(out,c("trend",
#                 "sigma_process",
#                 "sigma_migration",
#                 "sigma_daily",
#                 "migration_phenology_sd",
#                 "migration_phenology_mean",
#                 "log_rho_mu"), Rhat = TRUE,
#           pdf = FALSE)

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
# PART 4: GOODNESS OF FIT
# ***************************************************************
# ***************************************************************

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
  ggtitle("Observed vs Expected Seasonal Total Counts\n\nPost-breeding migration")
plot_obs_vs_expected

png(file = paste0(output_directory,"figures/Appendix_S1_migration_counts.png"), units = "in", width = 8, height = 8, res = 600)
plot_obs_vs_expected
dev.off()

plot_obs_vs_expected2 <- ggplot(data = expected_vs_observed)+
  geom_abline(slope=1,intercept=0)+
  geom_errorbar(aes(x = sum_obs, ymin = expected_q025,ymax = expected_q975, col = "Expected"), width = 0)+
  geom_point(aes(x = sum_obs, y = expected_q50, col = "Expected", shape = "Expected"))+
  geom_text(data = cor_obs_expected, aes(x = 0, y = max*1.1, label = cor_label), hjust = 0, size = 3)+
  geom_point(data = cor_obs_expected, aes(x = max*1.1, y = max*1.1), col = "transparent")+
  facet_wrap(station~., scales = "free")+
  scale_color_manual(values=c("gray50"), name = "", guide = FALSE)+
  scale_shape_manual(values=c(19,4), name = "", guide = FALSE)+
  xlab("Observed Total Count")+
  ylab("Expected Total Count")+
  ggtitle("Observed vs Expected Seasonal Total Counts\n\nPost-breeding migration")
plot_obs_vs_expected2

png(file = paste0(output_directory,"figures/Appendix_S1_obs_vs_expected.png"), units = "in", width = 8, height = 8, res = 600)
plot_obs_vs_expected2
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

station_assignment_plot <- ggplot(data = station_assignments, 
                                  aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  geom_text(data = station_fixed_rho, aes(x = 2000, y = max(station_assignments$n,na.rm = TRUE)*1.5, label = Fixed),
            hjust=0,vjust=1, fontface = "italic", size = 2)+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station~.)+
  xlim(c(range(station_assignments$Year)))+
  ggtitle(focal_season)+
  ylab("Number of bird samples analyzed")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

station_assignment_plot

png(file = paste0(output_directory,"figures/Appendix_S2_breeding_origins.png"), units = "in", width = 8, height = 6, res = 600)
station_assignment_plot
dev.off()


regional_change_since_2000 <- data.frame()

for (j in 1:jags_data$nstrata){
  for (t in 1:jags_data$nyear){
    
    # ------------
    # Spring estimates
    # ------------
    
    log_change <- log(out$sims.list$X[,j,t]/out$sims.list$X[,j,1])
    
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
    
  }
}

regional_change_since_2000$Year = (2000:2018)[regional_change_since_2000$year_number]
regional_change_since_2000$Stratum = c("East","West")[regional_change_since_2000$stratum_number] %>% factor(levels = c("West","East"))

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
