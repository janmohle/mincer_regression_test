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
n_loop_est = 1000
n_loop_fix = 1000
oos_window_est = 5000
oos_window_fix = 10000
est_window = 750
white_adjust = "hc0"
cores = 3#4
seed = 78

###########################################
########### Size parameters ###############
###########################################

# Simulation functional parameter
omega_size = 0.005
alpha1_size = 0.05
beta1_size = 0.94

# Simulation specifications
var_spec_sim_size = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_size = list(armaOrder = c(0,0))
dist_spec_sim_size = 'norm'
fixed_pars_sim_size = list(mu = 0,
                           omega = omega_size,
                           alpha1 = alpha1_size,
                           beta1 = beta1_size)

# Estimation parameters
var_spec_est_size = var_spec_sim_size
mean_spec_est_size = mean_spec_sim_size
dist_spec_est_size = dist_spec_sim_size

############################################
########### Power parameters ###############
############################################

# Simulation functional parameter
omega_power = 0.005
alpha1_power = 0.01
beta1_power = 0.94
gamma1_power = 0.08

# Simulation specifications
var_spec_sim_power = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_power = list(armaOrder = c(0,0))
dist_spec_sim_power = 'norm'
fixed_pars_sim_power = list(mu = 0,
                            omega = omega_power,
                            alpha1 = alpha1_power,
                            beta1 = beta1_power,
                            gamma1 = gamma1_power)

# Estimation parameters
var_spec_est_power = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power = mean_spec_sim_power
dist_spec_est_power = dist_spec_sim_power

###########################################################
########### Mincer Regression Specifications ##############
###########################################################

mincer_spec <- list(formula = list(simple_shortfall = shortfall ~ 1,
                                   simple_return = Return ~ ES,
                                   variance_shortfall = shortfall ~ variance_t_min_1,
                                   variance_return = Return ~ variance_t_min_1 + ES,
                                   residual_sqrt_shortfall = shortfall ~ residual_t_min_1_quadr,
                                   residual_sqrt_return = Return ~ residual_t_min_1_quadr + ES,
                                   residual_sqrt0_shortfall = shortfall ~ residual_t_min_1_quadr_lower_0,
                                   residual_sqrt0_return = Return ~ residual_t_min_1_quadr_lower_0 + ES,
                                   full_shortfall = shortfall ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0,
                                   full_return = Return ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0 + ES,
                                   indicator_lower_0_shortfall = shortfall ~ indicator_residual_lower_0,
                                   indicator_lower_0_return = Return ~ indicator_residual_lower_0 + ES),
                    h0 = list(simple_shortfall = c('(Intercept) = 0'),
                              simple_return = c('(Intercept) = 0', 'ES = 1'),
                              variance_shortfall = c('(Intercept) = 0', 'variance_t_min_1 = 0'),
                              variance_return = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'ES = 1'),
                              residual_sqrt_shortfall = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0'),
                              residual_sqrt_return = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0', 'ES = 1'),
                              residual_sqrt0_shortfall = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0'),
                              residual_sqrt0_return = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'),
                              full_shortfall = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0'),
                              full_return = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'),
                              indicator_lower_0_shortfall = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0'),
                              indicator_lower_0_return = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0', 'ES = 1')))

##########################################
########### Loop execution ###############
##########################################

# Size with estimated parameters
result_size_est <- estimation_loop_par(n_loop=n_loop_est,
                                       est_window=est_window,
                                       oos_window=oos_window_est,
                                       tolerance_lvl=tolerance_lvl,
                                       var_spec_sim=var_spec_sim_size,
                                       mean_spec_sim=mean_spec_sim_size,
                                       dist_spec_sim=dist_spec_sim_size,
                                       fixed_pars_sim=fixed_pars_sim_size,
                                       estimate=TRUE,
                                       var_spec_est=var_spec_est_size,
                                       mean_spec_est=mean_spec_est_size,
                                       dist_spec_est=dist_spec_est_size,
                                       fixed_pars_est=NA,
                                       cores=cores,
                                       white_adjust=white_adjust,
                                       seed=seed,
                                       mincer_spec=mincer_spec)
result_size_est_matrix <- create_result_matrix(result_size_est)

# Size with fixed parameters
result_size_fix <- estimation_loop_par(n_loop=n_loop_fix,
                                       est_window=est_window,
                                       oos_window=oos_window_fix,
                                       tolerance_lvl=tolerance_lvl,
                                       var_spec_sim=var_spec_sim_size,
                                       mean_spec_sim=mean_spec_sim_size,
                                       dist_spec_sim=dist_spec_sim_size,
                                       fixed_pars_sim=fixed_pars_sim_size,
                                       estimate=FALSE,
                                       var_spec_est=var_spec_est_size,
                                       mean_spec_est=mean_spec_est_size,
                                       dist_spec_est=dist_spec_est_size,
                                       fixed_pars_est=fixed_pars_sim_size,
                                       cores=cores,
                                       white_adjust=white_adjust,
                                       seed=seed,
                                       mincer_spec=mincer_spec)
result_size_fix_matrix <- create_result_matrix(result_size_fix)

# Power with estimated parameters
result_power_est <- estimation_loop_par(n_loop=n_loop_est,
                                        est_window=est_window,
                                        oos_window=oos_window_est,
                                        tolerance_lvl=tolerance_lvl,
                                        var_spec_sim=var_spec_sim_power,
                                        mean_spec_sim=mean_spec_sim_power,
                                        dist_spec_sim=dist_spec_sim_power,
                                        fixed_pars_sim=fixed_pars_sim_power,
                                        estimate=TRUE,
                                        var_spec_est=var_spec_est_power,
                                        mean_spec_est=mean_spec_est_power,
                                        dist_spec_est=dist_spec_est_power,
                                        fixed_pars_est=NA,
                                        cores=cores,
                                        white_adjust=white_adjust,
                                        seed=seed,
                                        mincer_spec=mincer_spec)
result_power_est_matrix <- create_result_matrix(result_power_est)

# Power with fixed parameters
result_power_fix <- estimation_loop_par(n_loop=n_loop_fix,
                                        est_window=est_window,
                                        oos_window=oos_window_fix,
                                        tolerance_lvl=tolerance_lvl,
                                        var_spec_sim=var_spec_sim_power,
                                        mean_spec_sim=mean_spec_sim_power,
                                        dist_spec_sim=dist_spec_sim_power,
                                        fixed_pars_sim=fixed_pars_sim_power,
                                        estimate=FALSE,
                                        var_spec_est=var_spec_est_power,
                                        mean_spec_est=mean_spec_est_power,
                                        dist_spec_est=dist_spec_est_power,
                                        fixed_pars_est=fixed_pars_sim_size,
                                        cores=cores,
                                        white_adjust=white_adjust,
                                        seed=seed,
                                        mincer_spec=mincer_spec)
result_power_fix_matrix <- create_result_matrix(result_power_fix)