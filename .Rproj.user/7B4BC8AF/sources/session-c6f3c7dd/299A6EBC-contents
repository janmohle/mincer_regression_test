# Initial cleaning
rm(list = ls())
if (dev.cur() != 1) {
  dev.off()
}
cat('\14')

# Load libraries
#install.packages('rugarch')
#install.packages('tidyverse')
#install.packages('zoo')
#install.packages('sandwich')
#install.packages('lmtest')
#install.packages('foreach')
#install.packages('doParallel')
#install.packages('car')
#install.packages('progressr')
#library(progressr)
library(rugarch)
library(tidyverse)
library(zoo)
library(sandwich)
library(lmtest)
library(foreach)
library(doParallel)
library(car)

source('functions.R')

#############################################
########### Global parameters ###############
#############################################
tolerance_lvl = 0.05
n_loop_est = 4 #1000
n_loop_fix = 12 #1000
oos_window_est = 5000
oos_window_fix = 10000
est_window = 750
white_adjust = TRUE
cores = 4
seed = 123

###########################################
########### Size parameters ###############
###########################################

# Simulation functional parameter
omega_size = 0.005
alpha1_size = 0.05
beta1_size = 0.94
eta11_size = NA

# Simulation specifications
var_spec_sim_size = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_size = list(armaOrder = c(0,0))
dist_spec_sim_size = 'norm'

# Estimation parameters
var_spec_est_size = var_spec_sim_size
mean_spec_est_size = mean_spec_sim_size
dist_spec_est_size = dist_spec_sim_size

############################################
########### Power parameters ###############
############################################

# Simulation functional parameter
omega_power = 0.005
alpha1_power = 0.02
beta1_power = 0.94
eta11_power = 0.06

# Simulation specifications
var_spec_sim_power = list(model = 'fGARCH', garchOrder = c(1, 1), submodel = 'TGARCH')
mean_spec_sim_power = list(armaOrder = c(0,0))
dist_spec_sim_power = 'norm'

# Estimation parameters
var_spec_est_power = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power = mean_spec_sim_power
dist_spec_est_power = dist_spec_sim_power

##########################################
########### Loop execution ###############
##########################################

# Size with estimated parameters
result_size_est <- estimation_loop_par(n_loop=n_loop_est,
                                       est_window=est_window,
                                       oos_window=oos_window_est,
                                       tolerance_lvl=tolerance_lvl,
                                       estimate=TRUE,
                                       fixed_pars=NA,
                                       var_spec_sim=var_spec_sim_size,
                                       mean_spec_sim=mean_spec_sim_size,
                                       dist_spec_sim=dist_spec_sim_size,
                                       omega=omega_size,
                                       alpha1=alpha1_size,
                                       beta1=beta1_size,
                                       eta11=eta11_size,
                                       var_spec_est=var_spec_est_size,
                                       mean_spec_est=mean_spec_est_size,
                                       dist_spec_est=dist_spec_est_size,
                                       cores=cores,
                                       white_adjust=white_adjust,
                                       seed=seed)
start <- Sys.time()
# Size with fixed parameters
result_size_fix <- estimation_loop_par(n_loop=n_loop_fix,
                                       est_window=est_window,
                                       oos_window=oos_window_fix,
                                       tolerance_lvl=tolerance_lvl,
                                       estimate=FALSE,
                                       fixed_pars=list(mu = 0, omega = omega_size, alpha1 = alpha1_size, beta1 = beta1_size),
                                       var_spec_sim=var_spec_sim_size,
                                       mean_spec_sim=mean_spec_sim_size,
                                       dist_spec_sim=dist_spec_sim_size,
                                       omega=omega_size,
                                       alpha1=alpha1_size,
                                       beta1=beta1_size,
                                       eta11=eta11_size,
                                       var_spec_est=var_spec_est_size,
                                       mean_spec_est=mean_spec_est_size,
                                       dist_spec_est=dist_spec_est_size,
                                       cores=cores,
                                       white_adjust=white_adjust,
                                       seed=seed)
print(Sys.time() - start)
# Power with estimated parameters
result_power_est <- estimation_loop_par(n_loop=n_loop_est,
                                        est_window=est_window,
                                        oos_window=oos_window_est,
                                        tolerance_lvl=tolerance_lvl,
                                        estimate=TRUE,
                                        fixed_pars=NA,
                                        var_spec_sim=var_spec_sim_power,
                                        mean_spec_sim=mean_spec_sim_power,
                                        dist_spec_sim=dist_spec_sim_power,
                                        omega=omega_power,
                                        alpha1=alpha1_power,
                                        beta1=beta1_power,
                                        eta11=eta11_power,
                                        var_spec_est=var_spec_est_power,
                                        mean_spec_est=mean_spec_est_power,
                                        dist_spec_est=dist_spec_est_power,
                                        cores=cores,
                                        white_adjust=white_adjust,
                                        seed=seed)

# Power with fixed parameters
result_power_fix <- estimation_loop_par(n_loop=n_loop_fix,
                                        est_window=est_window,
                                        oos_window=oos_window_fix,
                                        tolerance_lvl=tolerance_lvl,
                                        estimate=FALSE,
                                        fixed_pars=list(mu = 0, omega = 0.005, alpha1 = 0.05, beta1 = 0.94),
                                        var_spec_sim=var_spec_sim_power,
                                        mean_spec_sim=mean_spec_sim_power,
                                        dist_spec_sim=dist_spec_sim_power,
                                        omega=omega_power,
                                        alpha1=alpha1_power,
                                        beta1=beta1_power,
                                        eta11=eta11_power,
                                        var_spec_est=var_spec_est_power,
                                        mean_spec_est=mean_spec_est_power,
                                        dist_spec_est=dist_spec_est_power,
                                        cores=cores,
                                        white_adjust=white_adjust,
                                        seed=seed)