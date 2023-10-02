# *******************************************************************
# *******************************************************************

# Simulate trends within each of two regions, and simulate migration monitoring data
# - reanalyze data to confirm that trends can be recovered
# - evaluate benefits of 

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
load(paste0(output_directory,"analysis_",focal_season,".RData"))

# ***************************************************************
# ***************************************************************
# SIMULATE NEW DATASET
# ***************************************************************
# ***************************************************************

# The jags script to fit the model
sink("./analysis/1_scripts/migration_model_sim.jags")
cat("
    model {
  
  # *********************************************************
  # Priors and likelihood
  # *********************************************************

  #---------------------------------------------
  # Model for breeding origins of migrants arriving at each station
  #---------------------------------------------

  for (s in 1:nstation){
    for (y in 1:nyear){
      for (j in 1:nstrata){
        M[j,s,y] <- X[j,y] * rho[j,s] * rho_fix[j,s]
      }
    }
    
  }
  
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
  tau_stationday <- pow(sigma_stationday,-2)
  migration_phenology_tau <- pow(migration_phenology_sd,-2)
  
  for (i in 1:nobs){
    
    mu[i] <- log(f[i]) + log(T_star[station[i],year[i]]) + log(net_hrs[i]) + site_effect[site[i]]*dummy_site[i]
    expected_count[i] <- exp(mu[i])

    # Likelihood for counts
    f[i] <-  exp(logdensity.norm(day[i], migration_phenology_mean[station[i]], migration_phenology_tau))
    
    # Add daily overdispersion
    log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday), tau_stationday) 
    count[i] ~ dpois(exp(log_lambda[i]))
    
  }
  
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


sim_results <- data.frame()

if (file.exists("analysis/2_output/simulation/sim_results.RData")) load(file = "analysis/2_output/simulation/sim_results.RData")

for (simulation_rep in 1:125){
  
  print(simulation_rep)
  
  set.seed(simulation_rep)
  if (simulation_rep %in% sim_results$Simulation_Rep) next
  if (file.exists("analysis/2_output/simulation/sim_results.RData")) load(file = "analysis/2_output/simulation/sim_results.RData")
  
  # ------------------
  # Simulate trajectories for each of two strata
  # ------------------
  X_sim <- matrix(0,nrow=jags_data$nstrata,ncol=jags_data$nyear)
  X_sim[,1] <- 1
  for (t in 2:jags_data$nyear) X_sim[,t] <- exp(log(X_sim[,t-1]) + rnorm(2,0,out$mean$sigma_proc))
  #matplot(t(X_sim), type = "l")
  
  # ------------------
  # Use JAGS machinery to simulate counts and isotope data
  # ------------------
  
  jags_data_sim <- jags_data
  
  # Remove existing data
  jags_data_sim$N_origin <- jags_data$N_origin * NA
  jags_data_sim$count <- jags_data$count * NA
  
  # Assume up to 20 isotope samples will be collected per year at each site
  jags_data_sim$N_station_sampled <- array(20,dim = dim(jags_data$N_station_sampled), dimnames = dimnames(jags_data$N_station_sampled))
  jags_data_sim$X <- X_sim
  
  # Use biologically plausible parameters
  jags_data_sim$rho <- array(sample(out$mean$rho,length(out$mean$rho),replace=TRUE),
                              dim=dim(out$mean$rho))
                              
                              
  jags_data_sim$sigma_rho <- out$mean$sigma_rho
  jags_data_sim$sigma_stationday <- out$mean$sigma_stationday
  jags_data_sim$migration_phenology_mean <- out$mean$migration_phenology_mean
  jags_data_sim$migration_phenology_sd <- out$mean$migration_phenology_sd
  jags_data_sim$site_effect <- out$mean$site_effect
  jags_data_sim$rho_fix[,] <- 1
  
  parameters.to.save = c("N_origin","count","rho","station_composition")
  
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
  
  # True trend in each stratum
  trends_true <- 100*((X_sim[,19]/X_sim[,1])^(1/(19-1))-1)
  trends_true
  
  station_composition_full <- out_sim$sims.list$station_composition %>% 
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
  
  #station_composition_plot
  
  # **************************************************************************************
  # **************************************************************************************
  # Loop through "isotope collection scenarios" and estimate trend based on those data
  # **************************************************************************************
  # **************************************************************************************
  
  for (isotope_scenario in c(1,2,3)){
    
    if(sum(sim_results$Isotope_Scenario == isotope_scenario & sim_results$Simulation_Rep == simulation_rep)>0) next

    sim_results[which(sim_results$isotope_scenario == isotope_scenario & sim_results$simulation_rep == simulation_rep)]
    
    if (isotope_scenario == 1) years_with_isotopes = 1
    if (isotope_scenario == 2) years_with_isotopes = seq(1,jags_data_sim$nyear,5)
    if (isotope_scenario == 3) years_with_isotopes = seq(1,jags_data_sim$nyear,1)
    
    jags_data_refit <- jags_data
    jags_data_refit$N_station_sampled <- jags_data_sim$N_station_sampled
    jags_data_refit$N_origin <- jags_data$N_origin * NA
    jags_data_refit$N_origin[,,years_with_isotopes] <- out_sim$sims.list$N_origin[1,,,years_with_isotopes]
    
    if (isotope_scenario == 1){
      jags_data_refit$N_origin <- jags_data$N_origin * NA
      jags_data_refit$N_origin[!is.na(jags_data$N_origin)] <- out_sim$sims.list$N_origin[!is.na(jags_data$N_origin)]
      
    }
    
    jags_data_refit$count <- out_sim$sims.list$count[1,]
    jags_data_refit$rho_fix[,] <- 1
    
    jags_data_refit$N_origin
    
    # ------------------
    # Fit model to simulated data; evaluate trends
    # ------------------
    parameters.to.save = c("X")
    
    inits <- NULL
    nsamp <- 2000
    nb <- 1000
    nt <- 1
    ni <- nb + nsamp*nt
    
    out_refit <- jags(data = jags_data_refit,
                      model.file = "analysis/1_scripts/migration_model.jags",
                      parameters.to.save = parameters.to.save,
                      inits = inits,
                      n.chains = 3,
                      n.thin = nt,
                      n.iter = ni,
                      n.burnin = nb,
                      parallel = TRUE)
    
    out_refit$mcmc.info$elapsed.mins # about 15 min
    
    # Assess convergence of X estimates
    max_Rhat <- max(unlist(out_refit$Rhat["X"]), na.rm = TRUE)
    
    # Estimated trend in each stratum
    trends_estimated <- 100*((out_refit$sims.list$X[,,19]/out_refit$sims.list$X[,,1])^(1/(19-1))-1)
    
    # Summarize information for this run
    sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                                Isotope_Scenario = isotope_scenario,
                                                Stratum = "East",
                                                trend_true = trends_true[1],
                                                trend_est_q025 = quantile(trends_estimated[,1],0.025),
                                                trend_est_q50 = quantile(trends_estimated[,1],0.50),
                                                trend_est_q975 = quantile(trends_estimated[,1],0.975),
                                                max_Rhat = max_Rhat))
 
    sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                                Isotope_Scenario = isotope_scenario,
                                                Stratum = "West",
                                                trend_true = trends_true[2],
                                                trend_est_q025 = quantile(trends_estimated[,2],0.025),
                                                trend_est_q50 = quantile(trends_estimated[,2],0.50),
                                                trend_est_q975 = quantile(trends_estimated[,2],0.975),
                                                max_Rhat = max_Rhat))
  }
  
 
  limits <- range(sim_results[,c("trend_est_q025","trend_est_q975")])
  result_plot <- ggplot(data = sim_results, aes(x = trend_true, y = trend_est_q50, ymin = trend_est_q025, ymax = trend_est_q975))+
    geom_errorbar(width=0)+
    geom_point()+
    geom_abline(slope=1,intercept=0)+
    facet_grid(Stratum~Isotope_Scenario)+
    coord_cartesian(ylim=limits,xlim=limits)+
    theme_bw()
  
  print(result_plot)
  
  save(sim_results, file = "analysis/2_output/simulation/sim_results.RData")
  
}


sim_results$Isotope_Description <- c("Current Sampling","Collect Feathers Every 5 Years","Collect Feathers Every Year")[sim_results$Isotope_Scenario] %>% factor(levels = c("Current Sampling","Collect Feathers Every 5 Years","Collect Feathers Every Year"))

limits <- range(sim_results[,c("trend_est_q025","trend_est_q975")])
result_plot <- ggplot(data = sim_results, aes(x = trend_true, y = trend_est_q50, ymin = trend_est_q025, ymax = trend_est_q975))+
  geom_errorbar(width=0)+
  geom_point()+
  geom_abline(slope=1,intercept=0)+
  facet_grid(Stratum~Isotope_Description)+
  coord_cartesian(ylim=limits,xlim=limits)+
  theme_bw()+
  xlab("Simulated ('True') Trend")+
  ylab("Estimated Trend")

print(result_plot)

# Compare confidence interval widths between scenarios
sim_results$ciw <- sim_results$trend_est_q975 - sim_results$trend_est_q025

scen_1 <- subset(sim_results, Isotope_Scenario == 1)
scen_2 <- subset(sim_results, Isotope_Scenario == 2)
scen_3 <- subset(sim_results, Isotope_Scenario == 3)

# Sampling isotopes every year
a = 100*(scen_3$ciw - scen_1$ciw)/scen_1$ciw
hist(a)

median(a)

# Sampling isotopes every 5 years
b = 100*(scen_2$ciw - scen_1$ciw)/scen_1$ciw
hist(b)

mean(b)
median(b)
