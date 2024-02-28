# *******************************************************************
# *******************************************************************
# Simulate trends within each of two regions, and simulate migration monitoring data
# - reanalyze data to confirm that trends can be recovered
# *******************************************************************
# *******************************************************************

#------------------------------------------------
# Load/install packages
#------------------------------------------------

my_packs <- c('tidyverse','jagsUI') 
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

# ***************************************************************
# ***************************************************************
# LOAD FITTED MODEL
# ***************************************************************
# ***************************************************************

#------------------------------------------------
# Load / process migration and isotope data
#------------------------------------------------

focal_season <- "Spring"
start_year <- 1998
end_year <- 2018

# Relevant directories / set working directory
data_directory <- paste0("0_data/")
output_directory <- paste0("1_output/",focal_season,"/")
figure_directory <- paste0(output_directory,"/figures/")
table_directory <- paste0(output_directory,"/tables/")

# ------------------
# Load workspace
# ------------------

# Load fitted model
load(paste0(output_directory,"wksp_",focal_season,".RData"))

# ***************************************************************
# ***************************************************************
# SIMULATE NEW DATASET
# ***************************************************************
# ***************************************************************

# Script used to simulate new data, based on values of X, rho, and variance components supplied as data
sink("migration_model_simulation.jags")
cat("
    model {
  
  #---------------------------------------------
  # Model for breeding origins of migrants arriving at each station
  #---------------------------------------------
  
  tau_rho <- pow(sigma_rho,-2)
  
  # Model fluctuations in numbers of migrants from each stratum in each year
  for (s in 1:nstation){

    for (j in 1:nstrata){
    
      for (y in 1:nyear){
         M[j,s,y] <- X[j,y] * rho[j,s] * rho_fix[j,s]
      }
      
    }

  }

  # Model for known origins of migrating birds
  for (s in 1:nstation){
  
    for (y in 1:nyear){

      # Total number of migrants in year [y] at station [s]
      sumM[s,y] <- sum(M[1:nstrata,s,y])
      T[s,y] ~ dlnorm(log(sumM[s,y]),tau_rho)
      
      # Proportion of birds from each stratum
      for (j in 1:nstrata){
        p[j,s,y] <- M[j,s,y]/T[s,y]
      }

      # Multinomial likelihood for observed breeding origins in sample of birds
      N_origin[1:nstrata,s,y] ~ dmulti(p[1:nstrata,s,y],N_station_sampled[s,y])

    } # y
  } # s

  # ---------------------------------------------
  # Model counts of migrants arriving at each station, within each season
  # ---------------------------------------------

  # Daily overdispersion in counts at each site (e.g., due to daily weather, stopover behavior, etc.)
  tau_stationday <- pow(sigma_stationday,-2)
  
  # Seasonal migration window
  for (s in 1:nstation){
    migration_phenology_tau[s] <- pow(migration_phenology_sd[s],-2)
  }

  for (i in 1:nobs){

    # Distribute migrants throughout the season
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[station[i]], migration_phenology_tau[station[i]]))
    mu[i] <- log(f[i]) + log(T[station[i],year[i]])
    
    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i] + log(net_hrs[i]), tau_stationday)
    count[i] ~ dpois(exp(log_lambda[i]))

  }
  
}
    ",fill = TRUE)
sink()

sim_results <- data.frame()
sim_file <- "1_output/simulation/sim_results.Rdata"
if (file.exists(sim_file)) load(file = sim_file)

# 181 onwards
for (simulation_rep in (1:250)){
  
  print(simulation_rep)
  
  set.seed(simulation_rep)
  
  if (file.exists(sim_file)) load(file = sim_file)
  
  if (nrow(sim_results) > 0){
    if (simulation_rep %in% sim_results$Simulation_Rep) next
  }
  
  # ------------------
  # Simulate true trends for each of two strata (drawing uniformly from -0.1 to 0.1 for each stratum)
  # ------------------
  
  slopes_true <- runif(2,-0.1,0.1)
  X_sim <- matrix(0,nrow=jags_data$nstrata,ncol=jags_data$nyear)
  for (j in 1:jags_data$nstrata){
    logX_sim <- slopes_true[j] * ((1:jags_data$nyear)-1)
    X_sim[j,] <- exp(logX_sim)
  }
  
  # ------------------
  # List that will store parameters for this simulation (use jags_data as template)
  # ------------------
  
  jags_data_sim <- jags_data
  jags_data_sim$X <- X_sim
  
  # ------------------
  # Specify biologically plausible values for various model parameters
  # ------------------
  
  # Magnitude of year-to-year variation in summed seasonal indices
  jags_data_sim$sigma_rho <- out$mean$sigma_rho
  
  # Magnitude of extra-Poisson day-to-day variation in counts at each site
  jags_data_sim$sigma_stationday <- out$mean$sigma_stationday
  
  # Migration phenology at each station
  jags_data_sim$migration_phenology_mean <- sample(out$mean$migration_phenology_mean, size = jags_data_sim$nstation,replace=TRUE)
  jags_data_sim$migration_phenology_sd <- sample(out$mean$migration_phenology_sd, size = jags_data_sim$nstation,replace=TRUE)
  
  # Assume no stations have fixed (known) catchments
  jags_data_sim$rho_fix[,] <- 1 # No stations are fixed to have 0 catchment from any stratum
  
  # Randomly select from a wide range of rho parameters among stations
  #  - Most sites capture relatively few birds (per unit effort), but some capture a lot
  jags_data_sim$rho <- array(
    exp(runif(length(jags_data_sim$rho),log(0.01),log(2))),
    dim=dim(jags_data_sim$rho))
  
  # ------------------
  # Use JAGS machinery to simulate counts and isotope data
  # ------------------
  
  # Remove empirical isotope and count data
  jags_data_sim$N_origin <- jags_data$N_origin * NA
  jags_data_sim$count <- jags_data$count * NA
  
  # Assume up to 20 isotope samples will be collected at each site
  jags_data_sim$N_station_sampled <- array(20,dim = dim(jags_data$N_station_sampled), dimnames = dimnames(jags_data$N_station_sampled))
  
  parameters.to.save = c("N_origin","count")
  
  inits <- NULL
  nsamp <- 1
  nb <- 1
  nt <- 1
  ni <- nb + nsamp*nt
  
  out_sim <- jags(data = jags_data_sim,
                  model.file = "migration_model_simulation.jags",
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
  # Summarize resultant seasonal totals (confirm they are reasonable)
  # **************************************************************************************
  # **************************************************************************************
  
  count_df_sim <- count_df
  count_df_sim$count <- out_sim$sims.list$count[1,]
  
  T_est <- count_df_sim %>%
    group_by(station_number,year_abs) %>%
    summarize(T = sum(count/(exp(0.5*0.5^2)*net_hrs))) %>%
    group_by(station_number) %>%
    summarize(T1 = median(T))
  
  # **************************************************************************************
  # **************************************************************************************
  # Use simulated data to see if "true slopes" can be recovered
  # **************************************************************************************
  # **************************************************************************************
  
  # Create new jags data object, using simulated data only
  jags_data_refit <- jags_data
  
  jags_data_refit$N_station_sampled <- jags_data_sim$N_station_sampled
  
  # Replace empirical N_origin data with simulated N_origin data
  jags_data_refit$N_origin <- jags_data$N_origin * NA
  jags_data_refit$N_origin[!is.na(jags_data$N_origin)] <- out_sim$sims.list$N_origin[!is.na(jags_data$N_origin)]
  
  # Replace empirical count data with simulated count data
  jags_data_refit$count <- out_sim$sims.list$count[1,]
  
  # Assume catchment is unknown for every station (none are fixed to zero)
  jags_data_refit$rho_fix[,] <- 1
  
  rho_init <- jags_data$rho_fix * NA
  for (s in 1:ncol(rho_fix)) rho_init[,s] <- T_est$T1[T_est$station_number == s]
  
  inits <- function(){list(rho = rho_init)}
  
  # ------------------
  # Fit model to simulated data; see if trends can be recovered
  # ------------------
  
  parameters.to.save = c("slope",
                         "sigma_X",
                         "X",
                         "sigma_rho",
                         "sigma_stationday",
                         "migration_phenology_sd",
                         "migration_phenology_mean",
                         "rho",
                         "M",
                         "T")
  
  inits <- NULL
  nsamp <- 1000
  nb <- 1000
  nt <- 5
  ni <- nb + nsamp*nt
  
  out_refit <- jags(data = jags_data_refit,
                    model.file = "migration_model.jags",
                    parameters.to.save = parameters.to.save,
                    inits = inits,
                    n.chains = 3,
                    n.thin = nt,
                    n.iter = ni,
                    n.burnin = nb,
                    parallel = TRUE)
  
  print(out_refit$mcmc.info$elapsed.mins) # 10 min
  
  # Assess convergence
  latent_parameters <- c("slope",
                         "sigma_proc",
                         "sigma_rho",
                         "sigma_stationday",
                         "sigma_X",
                         "X")
  
  # Rhat
  Rhats <- unlist(out_refit$Rhat[latent_parameters])
  which(Rhats>1.1)
  max_Rhat <- max(Rhats, na.rm = TRUE)
  
  # ----------------------------------------------------
  # Save results from this run
  # ----------------------------------------------------
  
  if (file.exists(sim_file)) load(file = sim_file)
  
  # Summarize information for this run
  sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                              Stratum = "East",
                                              slope_true = slopes_true[1],
                                              slope_est_q025 = out_refit$q2.5$slope[1],
                                              slope_est_q50 = out_refit$q50$slope[1],
                                              slope_est_q975 = out_refit$q97.5$slope[1],
                                              max_Rhat = max_Rhat))
  
  sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                              Stratum = "West",
                                              
                                              slope_true = slopes_true[2],
                                              slope_est_q025 = out_refit$q2.5$slope[2],
                                              slope_est_q50 = out_refit$q50$slope[2],
                                              slope_est_q975 = out_refit$q97.5$slope[2],
                                              max_Rhat = max_Rhat))
  
  save(sim_results, file = sim_file)
  
  # ----------------------------------------------------
  # Plot results (removing runs that did not converge)
  # ----------------------------------------------------
  
  sim_results <- subset(sim_results, max_Rhat <= 1.2)
  sim_results$Stratum <- factor(sim_results$Stratum, levels = c("West","East"))
  sim_results$cov <- sim_results$slope_est_q025 < sim_results$slope_true & sim_results$slope_est_q975 > sim_results$slope_true
  
  limits <- range(sim_results[,c("slope_est_q025","slope_est_q975","slope_true")])
  
  result_plot <- ggplot(data = sim_results, 
                        aes(x = slope_true, 
                            y = slope_est_q50, ymin = slope_est_q025, ymax = slope_est_q975,
                            col = factor(cov)))+
    geom_errorbar(width=0)+
    geom_point()+
    geom_abline(slope=1,intercept=0)+
    scale_color_manual(values = c("orangered","dodgerblue"), guide = FALSE)+
    coord_cartesian(ylim=limits,xlim=limits)+
    ylab("Estimated log-linear trend")+
    xlab("True log-linear trend")+
    facet_grid(.~Stratum)+
    theme_bw()
  
  print(result_plot)
  
  mean(sim_results$cov)
  mean(sim_results$slope_est_q50 - sim_results$slope_true)
  
}

# ------------------------------------------------------
# Summarize results of repeated simulations
# ------------------------------------------------------
if (file.exists(sim_file)) load(file = sim_file)


sim_results <- subset(sim_results, max_Rhat <= 1.1)
sim_results$Stratum <- factor(sim_results$Stratum, levels = c("West","East"))
sim_results$cov <- sim_results$slope_est_q025 < sim_results$slope_true & sim_results$slope_est_q975 > sim_results$slope_true

limits <- range(sim_results[,c("slope_est_q025","slope_est_q975","slope_true")])

result_plot <- ggplot(data = sim_results, 
                      aes(x = slope_true, 
                          y = slope_est_q50, ymin = slope_est_q025, ymax = slope_est_q975,
                          col = factor(cov)))+
  geom_errorbar(width=0)+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  scale_color_manual(values = c("orangered", "black"), guide = FALSE)+
  coord_cartesian(ylim=limits,xlim=limits)+
  ylab("Estimated log-linear trend")+
  xlab("True log-linear trend")+
  facet_grid(.~Stratum)+
  theme_bw()

print(result_plot)

mean(sim_results$cov)
mean(sim_results$slope_est_q50 - sim_results$slope_true)

png(file ="1_output/Results_Appendix/Simulation.png", units = "in", width = 8, height = 4, res = 600)
print(result_plot)
dev.off()