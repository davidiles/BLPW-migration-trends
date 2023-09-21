# *******************************************************************
# *******************************************************************

# Analysis of spring migration data

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
# LOAD FITTED MODEL
# ***************************************************************
# ***************************************************************


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

# ------------------
# Load workspace
# ------------------

# Load fitted model
load(paste0(output_directory,"analysis_",focal_season,"_0.RData"))

# ***************************************************************
# ***************************************************************
# SIMULATE NEW DATASET
# ***************************************************************
# ***************************************************************

# JAGS script to simulate new data
sink("./analysis/1_scripts/migration_model_sim.jags")
cat("
    model {
  
  # *********************************************************
  # Priors and likelihood
  # *********************************************************

  #---------------------------------------------
  # Model for population dynamics in each region
  #---------------------------------------------
  
  tau_proc <- pow(sigma_proc,-2)
    
  for (j in 1:nstrata){
    
    # Work forwards from baseline year
    for (y in 1:nyear){

      logX[j,y] ~ dnorm( trend[j] * (y-1) , tau_proc)
      X[j,y] <- exp(logX[j,y])

    }
  
  } # j
  
  #---------------------------------------------
  # Model for breeding origins of migrants arriving at each station
  #---------------------------------------------
  
  tau_stationyear <- pow(sigma_stationyear,-2)

  for (s in 1:nstation){

    for (j in 1:nstrata){
      
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

  migration_phenology_tau <- pow(migration_phenology_sd,-2)
  
  for (s in 1:nstation){
    migration_phenology_mean[s] ~ dunif(1,360)
    
  }
  
  # Site-level fixed effects (for LPBO that contains multiple sub-stations)
  for (k in 1:nsite){

    # Daily overdispersion in counts at each site (e.g., due to daily weather)
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
    
  }
  
}
    ",fill = TRUE)
sink()


# ------------------
# Use JAGS machinery to simulate new data
# ------------------

jags_data_sim <- jags_data

# Remove existing data
jags_data_sim$N_origin <- jags_data$N_origin * NA
jags_data_sim$count <- jags_data$count * NA

# Assume up to 10 isotope samples will be collected per year at each site
jags_data_sim$N_station_sampled <- array(10,dim = dim(jags_data$N_station_sampled), dimnames = dimnames(jags_data$N_station_sampled))

# Specify biologically plausible population parameters
jags_data_sim$trend <- out$mean$trend
jags_data_sim$sigma_proc <- out$mean$sigma_proc
jags_data_sim$sigma_stationyear <- out$mean$sigma_stationyear
jags_data_sim$rho_mu <- exp(out$mean$log_rho_mu)
jags_data_sim$migration_phenology_sd <- out$mean$migration_phenology_sd
jags_data_sim$site_effect <- rep(0,jags_data$nsite)
jags_data_sim$sigma_stationday <- out$mean$sigma_stationday

parameters.to.save = c("N_origin","count","X")

inits <- NULL
nsamp <- 1
nb <- 1
nt <- 1
ni <- nb + nsamp*nt

out_sim <- jags(data = jags_data_sim,
            model.file = "analysis/1_scripts/migration_model_sim.jags",
            parameters.to.save = parameters.to.save,
            inits = inits,
            n.chains = 1,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = FALSE)

out_sim$mcmc.info$elapsed.mins

# **************************************************************************************
# **************************************************************************************
# Loop through "isotope collection scenarios" and estimate trend based on those data
# **************************************************************************************
# **************************************************************************************

isotope_scenario <- 3

if (isotope_scenario == 1) years_with_isotopes = 1
if (isotope_scenario == 2) years_with_isotopes = seq(1,jags_data_sim$nyear,5)
if (isotope_scenario == 3) years_with_isotopes = seq(1,jags_data_sim$nyear,1)

jags_data_refit <- jags_data
jags_data_refit$N_station_sampled <- jags_data_sim$N_station_sampled
jags_data_refit$N_origin <- jags_data$N_origin * NA
jags_data_refit$N_origin[,,years_with_isotopes] <- out_sim$sims.list$N_origin[1,,,years_with_isotopes]

jags_data_refit$count <- out_sim$sims.list$count[1,]

# ------------------
# Fit model to simulated data; evaluate trends
# ------------------
parameters.to.save = c("X")

inits <- NULL
nsamp <- 1000
nb <- 1000
nt <- 5
ni <- nb + nsamp*nt

out_refit <- jags(data = jags_data_refit,
                model.file = "analysis/1_scripts/migration_model_0.jags",
                parameters.to.save = parameters.to.save,
                inits = inits,
                n.chains = 3,
                n.thin = nt,
                n.iter = ni,
                n.burnin = nb,
                parallel = TRUE)
