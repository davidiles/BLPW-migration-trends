# *******************************************************************
# *******************************************************************

# Simulate trends within each of two regions, and simulate migration monitoring data
# - reanalyze data to confirm that trends can be recovered

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
start_year <- 1998
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
load(paste0(output_directory,"analysis_",focal_season,"_revised2.RData"))

# ***************************************************************
# ***************************************************************
# SIMULATE NEW DATASET
# ***************************************************************
# ***************************************************************

sink("./analysis/1_scripts/migration_model_sim.jags")
cat("
    model {
  
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

      T_star[s,y] ~ dlnorm(log(T[s,y]),tau_rho)
      
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
    # log_lambda[i] ~ dnorm(mu[i] - 1/(2*tau_stationday), tau_stationday) 
    log_lambda[i] ~ dnorm(mu[i], tau_stationday) 
    
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
sim_file <- "analysis/2_output/simulation/sim_results_nov17_randomwalk.Rdata"
if (file.exists(sim_file)) load(file = sim_file)

for (simulation_rep in 1:250){
  
  print(simulation_rep)
  
  set.seed(simulation_rep)
  
  # ------------------
  # Simulate trajectories for each of two strata
  # ------------------
  
  X_sim <- matrix(1,nrow=jags_data$nstrata,ncol=jags_data$nyear)
  for (j in 1:jags_data$nstrata){
    for (t in 2:jags_data$nyear){
      X_sim[j,t] <- exp(log(X_sim[j,t-1]) + rnorm(1,0,0.3))
    }
  }
  
  # ------------------
  # Calculate "true" slopes
  # ------------------
  
  slopes_true <- c()
  yr_vec <- 1:jags_data$nyear
  for (j in 1:jags_data$nstrata){
    slopes_true[j] <- lm(log(X_sim[j,])~yr_vec)$coef[2]
  }
  
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
  jags_data_sim$rho <- array(
    runif(length(out$mean$rho),0,2),
    dim=dim(out$mean$rho))
  
  jags_data_sim$sigma_rho <- out$mean$sigma_rho
  jags_data_sim$sigma_stationday <- out$mean$sigma_stationday
  jags_data_sim$migration_phenology_mean <- out$mean$migration_phenology_mean
  jags_data_sim$migration_phenology_sd <- out$mean$migration_phenology_sd
  jags_data_sim$site_effect <- out$mean$site_effect
  jags_data_sim$rho_fix[,] <- 1
  
  parameters.to.save = c("N_origin","count","rho","T_star","station_composition")
  
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
  trends_true <- 100*((X_sim[,jags_data$nyear]/X_sim[,1])^(1/(jags_data$nyear-1))-1)
  trends_true
  
  T_star <- out_sim$sims.list$T_star[1,,] %>% 
    reshape2::melt() %>%
    rename(station_number = Var1, year_number = Var2, T_star_true = value)
  T_star$Station = station_names[T_star$station_number]
  T_star$Year = year_vec[T_star$year_number]
  
  # Set sensible priors for migration parameters
  count_df$count <- out_sim$sims.list$count[1,]
  count_summary <- count_df %>%
    group_by(station,site_number,year_abs) %>%
    summarize(T_star_obs = sum(count/net_hrs)) %>%
    group_by(station,year_abs) %>%
    summarize(T_star_obs = mean(T_star_obs)) %>%
    ungroup() %>%
    rename(Station = station, Year = year_abs) %>%
    left_join(.,T_star)
  
  # **************************************************************************************
  # **************************************************************************************
  # Loop through "isotope collection scenarios" and estimate trend based on those data
  # **************************************************************************************
  # **************************************************************************************
  
  for (isotope_scenario in c(1,2)){
    
    # Skip if already complete
    if (file.exists(sim_file)) load(file = sim_file)
    tmp <- subset(sim_results,Simulation_Rep == simulation_rep & Isotope_Scenario == isotope_scenario)
    if(nrow(tmp>0)) next
    
    if(sum(sim_results$Isotope_Scenario == isotope_scenario & sim_results$Simulation_Rep == simulation_rep)>0) next
    
    print(paste0("Isotope scenario ",isotope_scenario))
    sim_results[which(sim_results$isotope_scenario == isotope_scenario & sim_results$simulation_rep == simulation_rep)]
    
    years_with_isotopes <- 1
    if (isotope_scenario == 2) years_with_isotopes = seq(1,jags_data_sim$nyear,10)
    
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
    parameters.to.save = c("slope","X")
    
    inits <- NULL
    nsamp <- 1000
    nb <- 1000
    nt <- 5
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
    
    print(out_refit$mcmc.info$elapsed.mins) # about 5 min
    
    # Assess convergence of X estimates
    max_Rhat <- max(unlist(out_refit$Rhat["X"]), na.rm = TRUE)
    
    # Estimated trend in each stratum
    trends_estimated <- 100*((out_refit$sims.list$X[,,jags_data$nyear]/out_refit$sims.list$X[,,1])^(1/(jags_data$nyear-1))-1)
    
    if (file.exists(sim_file)) load(file = sim_file)
    
    # Summarize information for this run
    sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                                Isotope_Scenario = isotope_scenario,
                                                Stratum = "East",
                                                
                                                slope_true = slopes_true[1],
                                                slope_est_q025 = out_refit$q2.5$slope[1],
                                                slope_est_q50 = out_refit$q50$slope[1],
                                                slope_est_q975 = out_refit$q97.5$slope[1],
                                                
                                                trend_true = trends_true[1],
                                                trend_est_q025 = quantile(trends_estimated[,1],0.025),
                                                trend_est_q50 = quantile(trends_estimated[,1],0.50),
                                                trend_est_q975 = quantile(trends_estimated[,1],0.975),
                                                max_Rhat = max_Rhat))
    
    sim_results <- rbind(sim_results,data.frame(Simulation_Rep = simulation_rep,
                                                Isotope_Scenario = isotope_scenario,
                                                Stratum = "West",
                                                
                                                slope_true = slopes_true[2],
                                                slope_est_q025 = out_refit$q2.5$slope[2],
                                                slope_est_q50 = out_refit$q50$slope[2],
                                                slope_est_q975 = out_refit$q97.5$slope[2],
                                                
                                                trend_true = trends_true[2],
                                                trend_est_q025 = quantile(trends_estimated[,2],0.025),
                                                trend_est_q50 = quantile(trends_estimated[,2],0.50),
                                                trend_est_q975 = quantile(trends_estimated[,2],0.975),
                                                max_Rhat = max_Rhat))
    
    save(sim_results, file = sim_file)
    
  }
  
  
  limits <- range(sim_results[,c("slope_est_q025","slope_est_q975","slope_true")])
  result_plot <- ggplot(data = sim_results, aes(x = slope_true, y = slope_est_q50, ymin = slope_est_q025, ymax = slope_est_q975))+
    geom_errorbar(width=0)+
    geom_point()+
    geom_abline(slope=1,intercept=0)+
    facet_grid(Stratum~Isotope_Scenario)+
    coord_cartesian(ylim=limits,xlim=limits)+
    theme_bw()
  
  print(result_plot)
  
}

# ------------------------------------------------------
# Summarize results of repeated simulations
# ------------------------------------------------------

sim_results <- subset(sim_results, max_Rhat <= 1.1)
sim_results$Scenario <- c("Collect Feathers Once","Collect Feathers Every 5 Years", "Collect Feathers Every Year")[sim_results$Isotope_Scenario] %>% factor(levels = c("Collect Feathers Once","Collect Feathers Every 5 Years","Collect Feathers Every Year"))

limits <- range(sim_results[,c("trend_est_q025","trend_est_q975")])
result_plot <- ggplot(data = sim_results, aes(x = trend_true, y = trend_est_q50, ymin = trend_est_q025, ymax = trend_est_q975))+
  geom_abline(trend=1,intercept=0, col = "gray75")+
  geom_errorbar(width=0)+
  geom_point()+
  
  facet_grid(Stratum~Scenario)+
  coord_cartesian(ylim=limits,xlim=limits)+
  theme_bw()+
  xlab("Simulated ('True') Trend\n\n")+
  ylab("Estimated Trend\n\n")+
  ggtitle("Simulation Results")

print(result_plot)

png(file ="analysis/2_output/Figures_Appendix/Appendix_Figure_Simulation.png", units = "in", width = 8, height = 6, res = 600)
print(result_plot)
dev.off()

# Mean bias
bias_coverage <- sim_results %>% 
  group_by(Scenario) %>%
  summarize('Mean Bias in Slope Estimate' = round(mean(slope_est_q50 - slope_true),3),
            '95% Interval Coverage' = round(mean(slope_est_q025<slope_true & slope_est_q975>slope_true),3),
            'Mean Width of 95% CI' = round(mean(slope_est_q975 - slope_est_q025),1))
bias_coverage
write.csv(bias_coverage, file ="analysis/2_output/Figures_Appendix/Simulation_Summary.csv",row.names = FALSE)

# ------------------------------------------------------
# Example figures for appendix, to accompany simulation results
# ------------------------------------------------------

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

set.seed(1)

# Simulate trajectories for each of two strata
X_sim <- matrix(0,nrow=jags_data$nstrata,ncol=jags_data$nyear)
X_sim[,1] <- 1
for (t in 2:jags_data$nyear) X_sim[,t] <- exp(log(X_sim[,t-1]) + rnorm(2,0,0.2))

jags_data_sim <- jags_data

# Remove existing data
jags_data_sim$N_origin <- jags_data$N_origin * NA
jags_data_sim$count <- jags_data$count * NA

# Assume up to 20 isotope samples will be collected per year at each site
jags_data_sim$N_station_sampled <- array(20,dim = dim(jags_data$N_station_sampled), dimnames = dimnames(jags_data$N_station_sampled))
jags_data_sim$X <- X_sim

# Use biologically plausible parameters
jags_data_sim$rho <- array(
  runif(length(out$mean$rho),0,2),
  #rlnorm(length(out$mean$rho),-2,1),
  
  dim=dim(out$mean$rho))

jags_data_sim$sigma_rho <- out$mean$sigma_rho
jags_data_sim$sigma_stationday <- out$mean$sigma_stationday
jags_data_sim$migration_phenology_mean <- out$mean$migration_phenology_mean
jags_data_sim$migration_phenology_sd <- out$mean$migration_phenology_sd
jags_data_sim$site_effect <- out$mean$site_effect
jags_data_sim$rho_fix[,] <- 1

parameters.to.save = c("N_origin","count","rho","T_star","station_composition")

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
trends_true <- 100*((X_sim[,jags_data$nyear]/X_sim[,1])^(1/(jags_data$nyear-1))-1)
trends_true

Year <- 1:jags_data$nyear
Fig_S5_1 <- ggplot()+
  geom_line(aes(x = Year, y = X_sim[1,], col = "Stratum 1 ('West')"), linewidth = 2)+
  geom_line(aes(x = Year, y = X_sim[2,], col = "Stratum 2 ('East')"), linewidth = 2)+
  scale_color_manual(values=strata_colours, name = "Stratum")+
  coord_cartesian(ylim=c(0,max(X_sim)))+
  ylab("Index of abundance")

png(file ="analysis/2_output/Figures_Appendix/Fig_S5_1.png", units = "in", width = 6, height = 4, res = 600)
print(Fig_S5_1)
dev.off()

# ---------------------------------------
# Illustrate daily counts, at each of two stations, across 5 years of data...
# ---------------------------------------
count_df_sim <- count_df %>% 
  mutate(count = out_sim$sims.list$count[1,],
         station_name = paste0("Station ",station_number))

sim_summary <- count_df_sim %>%
  group_by(station_name) %>%
  summarize(mean_count = mean(count)) %>%
  arrange(mean_count)

Fig_S5_2 <- ggplot(data = subset(count_df_sim, station_number %in% c(10,4) & year_abs %in% 2000:2018))+
  geom_point(aes(x = day_number, y = count))+
  facet_grid(station_name~year_abs, scales = "free_y")

Fig_S5_2

png(file ="analysis/2_output/Figures_Appendix/Fig_S5_2.png", units = "in", width = 20, height = 4, res = 600)
print(Fig_S5_2)
dev.off()

# ---------------------------------------
# Illustrate daily counts, at each of two stations, across 5 years of data...
# ---------------------------------------
count_df_sim <- count_df %>% 
  mutate(count = out_sim$sims.list$count[1,],
         station_name = paste0("Station ",station_number))

sim_summary <- count_df_sim %>%
  group_by(station_name) %>%
  summarize(mean_count = mean(count)) %>%
  arrange(mean_count)

Fig_S5_2 <- ggplot(data = subset(count_df_sim, station_number %in% c(10,4) & year_abs %in% 2000:2018))+
  geom_point(aes(x = day_number, y = count))+
  facet_grid(station_name~year_abs, scales = "free_y")

Fig_S5_2

png(file ="analysis/2_output/Figures_Appendix/Fig_S5_2.png", units = "in", width = 20, height = 4, res = 600)
print(Fig_S5_2)
dev.off()

# ---------------------------------------
# Illustrate isotope sampling for scenario 1
# ---------------------------------------
years_with_isotopes = 1
jags_data_sim$N_origin <- jags_data$N_origin * NA
jags_data_sim$N_origin[,,years_with_isotopes] <- out_sim$sims.list$N_origin[1,,,years_with_isotopes]

station_assignments <- jags_data_sim$N_origin %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, Year = Var3, n = value)
station_assignments$Stratum = factor(station_assignments$Stratum, levels = c("West","East"))
station_assignments$Station_Number <- factor(station_assignments$Station) %>% as.numeric() %>% paste0("Station ",.)

Fig_S5_3 <- ggplot(data = station_assignments, 
                   aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station_Number~.)+
  coord_cartesian(xlim=c(1998,2018))+
  ylab("Number of bird samples analyzed")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

png(file ="analysis/2_output/Figures_Appendix/Fig_S5_3.png", units = "in", width = 10, height = 8, res = 600)
print(Fig_S5_3)
dev.off()

# ---------------------------------------
# Illustrate isotope sampling for scenario 2
# ---------------------------------------

years_with_isotopes = seq(1,jags_data_sim$nyear,5)
jags_data_sim$N_origin <- jags_data$N_origin * NA
jags_data_sim$N_origin[,,years_with_isotopes] <- out_sim$sims.list$N_origin[1,,,years_with_isotopes]

station_assignments <- jags_data_sim$N_origin %>%
  reshape2::melt() %>%
  rename(Stratum = Var1, Station = Var2, Year = Var3, n = value)
station_assignments$Stratum = factor(station_assignments$Stratum, levels = c("West","East"))
station_assignments$Station_Number <- factor(station_assignments$Station) %>% as.numeric() %>% paste0("Station ",.)

Fig_S5_4 <- ggplot(data = station_assignments, 
                   aes(x = Year, y = n, fill = Stratum)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = strata_colours, name = "Stratum of origin")+
  facet_wrap(Station_Number~.)+
  coord_cartesian(xlim=c(1998,2018))+
  ylab("Number of bird samples analyzed")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

png(file ="analysis/2_output/Figures_Appendix/Fig_S5_4.png", units = "in", width = 10, height = 8, res = 600)
print(Fig_S5_4)
dev.off()
