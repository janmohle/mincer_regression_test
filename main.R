# Initial cleaning
rm(list = ls())
if (dev.cur() != 1) {
  dev.off()
}
cat('\14')

# Load libraries
#install.packages('rugarch') not installed
#install.packages('tidyverse')
#install.packages('zoo')
#install.packages('sandwich')
#install.packages('lmtest')
#install.packages('foreach')
#install.packages('doParallel')
#install.packages('car')
#install.packages('esback')
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
library(esback)

source('functions.R')

#############################################
########### Global parameters ###############
#############################################
result_txt_file = 'results.txt'
tolerance_lvl = 0.05
n_loop_est = 1000
n_loop_fix = 1000
oos_window_est = 5000
oos_window_fix = 10000
est_window = 750
cores = detectCores() * 0.75
seed = 78
power_loops = 1:8
white_adjust = c('hc3')#, FALSE)
empirical=TRUE

################################################
########### Size loop parameters ###############
################################################

# Simulation functional parameters
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

###################################################
########### Power loop 1 parameters ###############
###################################################

# Simulation functional parameters
omega_power_1 = 0.005
alpha1_power_1 = 0.03
beta1_power_1 = 0.94
gamma1_power_1 = 0.04

# Simulation specifications
var_spec_sim_power_1 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_power_1 = list(armaOrder = c(0,0))
dist_spec_sim_power_1 = 'norm'
fixed_pars_sim_power_1 = list(mu = 0,
                              omega = omega_power_1,
                              alpha1 = alpha1_power_1,
                              beta1 = beta1_power_1,
                              gamma1 = gamma1_power_1)

# Estimation parameters
var_spec_est_power_1 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power_1 = mean_spec_sim_power_1
dist_spec_est_power_1 = dist_spec_sim_power_1

###################################################
########### Power loop 2 parameters ###############
###################################################

# Simulation functional parameters
omega_power_2 = 0.005
alpha1_power_2 = 0.01
beta1_power_2 = 0.94
gamma1_power_2 = 0.08

# Simulation specifications
var_spec_sim_power_2 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_power_2 = list(armaOrder = c(0,0))
dist_spec_sim_power_2 = 'norm'
fixed_pars_sim_power_2 = list(mu = 0,
                              omega = omega_power_2,
                              alpha1 = alpha1_power_2,
                              beta1 = beta1_power_2,
                              gamma1 = gamma1_power_2)

# Estimation parameters
var_spec_est_power_2 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power_2 = mean_spec_sim_power_2
dist_spec_est_power_2 = dist_spec_sim_power_2

###################################################
########### Power loop 3 parameters ###############
###################################################

# Simulation functional parameter
omega_power_3 = 0.005
alpha1_power_3 = 0.05
beta1_power_3 = 0.94
shape_power_3 = 4

# Simulation specifications
var_spec_sim_power_3 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_power_3 = list(armaOrder = c(0,0))
dist_spec_sim_power_3 = 'std'
fixed_pars_sim_power_3 = list(mu = 0,
                              omega = omega_power_3,
                              alpha1 = alpha1_power_3,
                              beta1 = beta1_power_3,
                              shape = shape_power_3)

# Estimation parameters
var_spec_est_power_3 = var_spec_sim_power_3
mean_spec_est_power_3 = mean_spec_sim_power_3
dist_spec_est_power_3 = 'norm'

###################################################
########### Power loop 4 parameters ###############
###################################################

# Simulation functional parameter
omega_power_4 = 0.005
alpha1_power_4 = 0.03
beta1_power_4 = 0.94
gamma1_power_4 = 0.04
shape_power_4 = 4

# Simulation specifications
var_spec_sim_power_4 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_power_4 = list(armaOrder = c(0,0))
dist_spec_sim_power_4 = 'std'
fixed_pars_sim_power_4 = list(mu = 0,
                              omega = omega_power_4,
                              alpha1 = alpha1_power_4,
                              beta1 = beta1_power_4,
                              gamma1 = gamma1_power_4,
                              shape = shape_power_4)

# Estimation parameters
var_spec_est_power_4 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power_4 = mean_spec_sim_power_4
dist_spec_est_power_4 = 'norm'

###################################################
########### Power loop 5 parameters ###############
###################################################

# Simulation functional parameter
omega_power_5 = 0.005
alpha1_power_5 = 0.01
beta1_power_5 = 0.94
gamma1_power_5 = 0.08
shape_power_5 = 4

# Simulation specifications
var_spec_sim_power_5 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_power_5 = list(armaOrder = c(0,0))
dist_spec_sim_power_5 = 'std'
fixed_pars_sim_power_5 = list(mu = 0,
                              omega = omega_power_5,
                              alpha1 = alpha1_power_5,
                              beta1 = beta1_power_5,
                              gamma1 = gamma1_power_5,
                              shape = shape_power_5)

# Estimation parameters
var_spec_est_power_5 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_power_5 = mean_spec_sim_power_5
dist_spec_est_power_5 = 'norm'




##############TEST

###################################################
########### Power loop 6 parameters ###############
###################################################

# Simulation functional parameter
omega_power_6 = 0.005
alpha1_power_6 = 0.05
beta1_power_6 = 0.94
shape_power_6 = 6

# Simulation specifications
var_spec_sim_power_6 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_power_6 = list(armaOrder = c(0,0))
dist_spec_sim_power_6 = 'std'
fixed_pars_sim_power_6 = list(mu = 0,
                              omega = omega_power_6,
                              alpha1 = alpha1_power_6,
                              beta1 = beta1_power_6,
                              shape = shape_power_6)

# Estimation parameters
var_spec_est_power_6 = var_spec_sim_power_6
mean_spec_est_power_6 = mean_spec_sim_power_6
dist_spec_est_power_6 = 'norm'

###################################################
########### Power loop 7 parameters ###############
###################################################

# Simulation functional parameter
omega_power_7 = 0.005
alpha1_power_7 = 0.05
beta1_power_7 = 0.94
shape_power_7 = 30

# Simulation specifications
var_spec_sim_power_7 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_power_7 = list(armaOrder = c(0,0))
dist_spec_sim_power_7 = 'std'
fixed_pars_sim_power_7 = list(mu = 0,
                              omega = omega_power_7,
                              alpha1 = alpha1_power_7,
                              beta1 = beta1_power_7,
                              shape = shape_power_7)

# Estimation parameters
var_spec_est_power_7 = var_spec_sim_power_7
mean_spec_est_power_7 = mean_spec_sim_power_7
dist_spec_est_power_7 = 'norm'

###################################################
########### Power loop 8 parameters ###############
###################################################

# Simulation functional parameter
omega_power_8 = 0.005
alpha1_power_8 = 0.05
beta1_power_8 = 0.94
shape_power_8 = 100

# Simulation specifications
var_spec_sim_power_8 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_power_8 = list(armaOrder = c(0,0))
dist_spec_sim_power_8 = 'std'
fixed_pars_sim_power_8 = list(mu = 0,
                              omega = omega_power_8,
                              alpha1 = alpha1_power_8,
                              beta1 = beta1_power_8,
                              shape = shape_power_8)

# Estimation parameters
var_spec_est_power_8 = var_spec_sim_power_8
mean_spec_est_power_8 = mean_spec_sim_power_8
dist_spec_est_power_8 = 'norm'

###################################################
########### Size loop 2 parameters ###############
###################################################

# Simulation functional parameter
omega_size_2 = 0.005
alpha1_size_2 = 0.05
beta1_size_2 = 0.94
shape_size_2 = 4

# Simulation specifications
var_spec_sim_size_2 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_size_2 = list(armaOrder = c(0,0))
dist_spec_sim_size_2 = 'std'
fixed_pars_sim_size_2 = list(mu = 0,
                             omega = omega_size_2,
                             alpha1 = alpha1_size_2,
                             beta1 = beta1_size_2,
                             shape = shape_size_2)

# Estimation parameters
var_spec_est_size_2 = var_spec_sim_size_2
mean_spec_est_size_2 = mean_spec_sim_size_2
dist_spec_est_size_2 = dist_spec_sim_size_2



###########################################################
########### Mincer Regression Specifications ##############
###########################################################

mincer_spec <- list(simple_shortfall = list(formula = shortfall ~ 1,
                                            h0 = c('(Intercept) = 0')),
                    simple_return = list(formula = Return ~ ES,
                                         h0 = c('(Intercept) = 0', 'ES = 1')),
                    variance_shortfall = list(formula = shortfall ~ variance_t_min_1,
                                              h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0')),
                    # variance_return = list(formula = Return ~ variance_t_min_1 + ES,
                    #                        h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'ES = 1')),
                    residual_sqrt_shortfall = list(formula = shortfall ~ residual_t_min_1_quadr,
                                                   h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0')),
                    # residual_sqrt_return = list(formula = Return ~ residual_t_min_1_quadr + ES,
                    #                             h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0', 'ES = 1')),
                    indicator_lower_0_shortfall = list(formula = shortfall ~ indicator_residual_lower_0,
                                                       h0 = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0')),
                    # indicator_lower_0_return = list(formula = Return ~ indicator_residual_lower_0 + ES,
                    #                                 h0 = c('(Intercept) = 0', 'indicator_residual_lower_0 = 0', 'ES = 1')),
                    residual_sqrt0_shortfall = list(formula = shortfall ~ residual_t_min_1_quadr_lower_0,
                                                    h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0'))#,
                    # residual_sqrt0_return = list(formula = Return ~ residual_t_min_1_quadr_lower_0 + ES,
                    #                              h0 = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1')),
                    # full_shortfall = list(formula = shortfall ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0,
                    #                       h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0')),
                    # full_return = list(formula = Return ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0 + ES,
                    #                    h0 = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'))
                    )


##########################################
########### Loop execution ###############
##########################################

# Write date and time of execution start into txt result file
write(paste0('\n\n', Sys.time()), file = result_txt_file, append = TRUE)
write(paste0('\nEmpirical distribution = ', empirical), file = result_txt_file, append = TRUE)

# Size loops without estimated parameters
for(wa in white_adjust){
  result_size_fix_wa <- estimation_loop_par(n_loop=n_loop_fix,
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
                                            white_adjust=wa,
                                            seed=seed,
                                            mincer_spec=mincer_spec,
                                            execute_additional_tsts=TRUE,
                                            empirical=empirical)
  result_size_fix_wa_matrix <- create_result_matrix(result_size_fix_wa)

  assign(paste0('result_size_fix_', wa), result_size_fix_wa)
  assign(paste0('result_size_fix_', wa, '_matrix'), result_size_fix_wa_matrix)

  write_results_to_txt(name = paste0('result_size_fix_', wa), txt_file = result_txt_file)
}
rm(result_size_fix_wa, wa, result_size_fix_wa_matrix)


# Size loops 2 without estimated parameters
# for(wa in white_adjust){
#   result_size_fix_wa <- estimation_loop_par(n_loop=n_loop_fix,
#                                             est_window=est_window,
#                                             oos_window=oos_window_fix,
#                                             tolerance_lvl=tolerance_lvl,
#                                             var_spec_sim=var_spec_sim_size_2,
#                                             mean_spec_sim=mean_spec_sim_size_2,
#                                             dist_spec_sim=dist_spec_sim_size_2,
#                                             fixed_pars_sim=fixed_pars_sim_size_2,
#                                             estimate=FALSE,
#                                             var_spec_est=var_spec_est_size_2,
#                                             mean_spec_est=mean_spec_est_size_2,
#                                             dist_spec_est=dist_spec_est_size_2,
#                                             fixed_pars_est=fixed_pars_sim_size_2,
#                                             cores=cores,
#                                             white_adjust=wa,
#                                             seed=seed,
#                                             mincer_spec=mincer_spec,
#                                             execute_additional_tsts=TRUE,
#                                             empirical=empirical)
#   result_size_fix_wa_matrix <- create_result_matrix(result_size_fix_wa)
# 
#   assign(paste0('result_size_2_fix_', wa), result_size_fix_wa)
#   assign(paste0('result_size_2_fix_', wa, '_matrix'), result_size_fix_wa_matrix)
# 
#   write_results_to_txt(name = paste0('result_size_2_fix_', wa), txt_file = result_txt_file)
# }
# rm(result_size_fix_wa, wa, result_size_fix_wa_matrix)


# Power loop executions
for(i in power_loops){

  var_spec_sim_power_i <- get(paste0('var_spec_sim_power_', i))
  mean_spec_sim_power_i <- get(paste0('mean_spec_sim_power_', i))
  dist_spec_sim_power_i <- get(paste0('dist_spec_sim_power_', i))
  fixed_pars_sim_power_i <- get(paste0('fixed_pars_sim_power_', i))

  var_spec_est_power_i <- get(paste0('var_spec_est_power_', i))
  mean_spec_est_power_i <- get(paste0('mean_spec_est_power_', i))
  dist_spec_est_power_i <- get(paste0('dist_spec_est_power_', i))

  for(wa in white_adjust){
    result_power_i_fix_wa <- estimation_loop_par(n_loop=n_loop_fix,
                                                 est_window=est_window,
                                                 oos_window=oos_window_fix,
                                                 tolerance_lvl=tolerance_lvl,
                                                 var_spec_sim=var_spec_sim_power_i,
                                                 mean_spec_sim=mean_spec_sim_power_i,
                                                 dist_spec_sim=dist_spec_sim_power_i,
                                                 fixed_pars_sim=fixed_pars_sim_power_i,
                                                 estimate=FALSE,
                                                 var_spec_est=var_spec_est_power_i,
                                                 mean_spec_est=mean_spec_est_power_i,
                                                 dist_spec_est=dist_spec_est_power_i,
                                                 fixed_pars_est=fixed_pars_sim_size,
                                                 cores=cores,
                                                 white_adjust=wa,
                                                 seed=seed,
                                                 mincer_spec=mincer_spec,
                                                 execute_additional_tsts=TRUE,
                                                 empirical=empirical)

    result_power_i_fix_wa_matrix <- create_result_matrix(result_power_i_fix_wa)

    assign(paste0('result_power_', i, '_fix_', wa), result_power_i_fix_wa)
    assign(paste0('result_power_', i, '_fix_', wa, '_matrix'), result_power_i_fix_wa_matrix)

    write_results_to_txt(name = paste0('result_power_', i, '_fix_', wa), txt_file = result_txt_file)
  }
}
rm(var_spec_sim_power_i, mean_spec_sim_power_i, dist_spec_sim_power_i, fixed_pars_sim_power_i, var_spec_est_power_i, mean_spec_est_power_i, dist_spec_est_power_i, result_power_i_fix_wa_matrix, result_power_i_fix_wa, wa, i)



# empirical=empirical NEEDS TO BE ADDED

# Power loop 1 without estimated parameters (with HC3)
# result_power_1_fix_hc3 <- estimation_loop_par(n_loop=n_loop_fix,
#                                               est_window=est_window,
#                                               oos_window=oos_window_fix,
#                                               tolerance_lvl=tolerance_lvl,
#                                               var_spec_sim=var_spec_sim_power_1,
#                                               mean_spec_sim=mean_spec_sim_power_1,
#                                               dist_spec_sim=dist_spec_sim_power_1,
#                                               fixed_pars_sim=fixed_pars_sim_power_1,
#                                               estimate=FALSE,
#                                               var_spec_est=var_spec_est_power_1,
#                                               mean_spec_est=mean_spec_est_power_1,
#                                               dist_spec_est=dist_spec_est_power_1,
#                                               fixed_pars_est=fixed_pars_sim_size,
#                                               cores=cores,
#                                               white_adjust='hc3',
#                                               seed=seed,
#                                               mincer_spec=mincer_spec,
#                                               execute_additional_tsts=TRUE)
# result_power_1_fix_hc3_matrix <- create_result_matrix(result_power_1_fix_hc3)
# write_results_to_txt(name = 'result_power_1_fix_hc3', txt_file = result_txt_file)
# 
# 
# # Power loop 1 without estimated parameters (with HC0)
# result_power_1_fix_hc0 <- estimation_loop_par(n_loop=n_loop_fix,
#                                               est_window=est_window,
#                                               oos_window=oos_window_fix,
#                                               tolerance_lvl=tolerance_lvl,
#                                               var_spec_sim=var_spec_sim_power_1,
#                                               mean_spec_sim=mean_spec_sim_power_1,
#                                               dist_spec_sim=dist_spec_sim_power_1,
#                                               fixed_pars_sim=fixed_pars_sim_power_1,
#                                               estimate=FALSE,
#                                               var_spec_est=var_spec_est_power_1,
#                                               mean_spec_est=mean_spec_est_power_1,
#                                               dist_spec_est=dist_spec_est_power_1,
#                                               fixed_pars_est=fixed_pars_sim_size,
#                                               cores=cores,
#                                               white_adjust='hc0',
#                                               seed=seed,
#                                               mincer_spec=mincer_spec,
#                                               execute_additional_tsts=TRUE)
# result_power_1_fix_hc0_matrix <- create_result_matrix(result_power_1_fix_hc0)
# write_results_to_txt(name = 'result_power_1_fix_hc0', txt_file = result_txt_file)
# 
# 
# 
# # Power loop 2 without estimated parameters (with HC3)
# result_power_2_fix_hc3 <- estimation_loop_par(n_loop=n_loop_fix,
#                                               est_window=est_window,
#                                               oos_window=oos_window_fix,
#                                               tolerance_lvl=tolerance_lvl,
#                                               var_spec_sim=var_spec_sim_power_2,
#                                               mean_spec_sim=mean_spec_sim_power_2,
#                                               dist_spec_sim=dist_spec_sim_power_2,
#                                               fixed_pars_sim=fixed_pars_sim_power_2,
#                                               estimate=FALSE,
#                                               var_spec_est=var_spec_est_power_2,
#                                               mean_spec_est=mean_spec_est_power_2,
#                                               dist_spec_est=dist_spec_est_power_2,
#                                               fixed_pars_est=fixed_pars_sim_size,
#                                               cores=cores,
#                                               white_adjust='hc3',
#                                               seed=seed,
#                                               mincer_spec=mincer_spec,
#                                               execute_additional_tsts=TRUE)
# result_power_2_fix_hc3_matrix <- create_result_matrix(result_power_2_fix_hc3)
# write_results_to_txt(name = 'result_power_2_fix_hc3', txt_file = result_txt_file)
# 
# 
# # Power loop 2 without estimated parameters (with HC0)
# result_power_2_fix_hc0 <- estimation_loop_par(n_loop=n_loop_fix,
#                                               est_window=est_window,
#                                               oos_window=oos_window_fix,
#                                               tolerance_lvl=tolerance_lvl,
#                                               var_spec_sim=var_spec_sim_power_2,
#                                               mean_spec_sim=mean_spec_sim_power_2,
#                                               dist_spec_sim=dist_spec_sim_power_2,
#                                               fixed_pars_sim=fixed_pars_sim_power_2,
#                                               estimate=FALSE,
#                                               var_spec_est=var_spec_est_power_2,
#                                               mean_spec_est=mean_spec_est_power_2,
#                                               dist_spec_est=dist_spec_est_power_2,
#                                               fixed_pars_est=fixed_pars_sim_size,
#                                               cores=cores,
#                                               white_adjust='hc0',
#                                               seed=seed,
#                                               mincer_spec=mincer_spec,
#                                               execute_additional_tsts=TRUE)
# result_power_2_fix_hc0_matrix <- create_result_matrix(result_power_2_fix_hc0)
# write_results_to_txt(name = 'result_power_2_fix_hc0', txt_file = result_txt_file)
# 
# 
# 
# # Power loop 3 without estimated parameters
# result_power_3_fix <- estimation_loop_par(n_loop=n_loop_fix,
#                                           est_window=est_window,
#                                           oos_window=oos_window_fix,
#                                           tolerance_lvl=tolerance_lvl,
#                                           var_spec_sim=var_spec_sim_power_3,
#                                           mean_spec_sim=mean_spec_sim_power_3,
#                                           dist_spec_sim=dist_spec_sim_power_3,
#                                           fixed_pars_sim=fixed_pars_sim_power_3,
#                                           estimate=FALSE,
#                                           var_spec_est=var_spec_est_power_3,
#                                           mean_spec_est=mean_spec_est_power_3,
#                                           dist_spec_est=dist_spec_est_power_3,
#                                           fixed_pars_est=fixed_pars_sim_size,
#                                           cores=cores,
#                                           white_adjust=white_adjust,
#                                           seed=seed,
#                                           mincer_spec=mincer_spec,
#                                           execute_additional_tsts=TRUE)
# result_power_3_fix_matrix <- create_result_matrix(result_power_3_fix)
# write_results_to_txt(name = 'result_power_3_fix', txt_file = result_txt_file)
# 
# # Power loop 3 without estimated parameters (with HC0)
# result_power_3_fix_hc0 <- estimation_loop_par(n_loop=n_loop_fix,
#                                               est_window=est_window,
#                                               oos_window=oos_window_fix,
#                                               tolerance_lvl=tolerance_lvl,
#                                               var_spec_sim=var_spec_sim_power_3,
#                                               mean_spec_sim=mean_spec_sim_power_3,
#                                               dist_spec_sim=dist_spec_sim_power_3,
#                                               fixed_pars_sim=fixed_pars_sim_power_3,
#                                               estimate=FALSE,
#                                               var_spec_est=var_spec_est_power_3,
#                                               mean_spec_est=mean_spec_est_power_3,
#                                               dist_spec_est=dist_spec_est_power_3,
#                                               fixed_pars_est=fixed_pars_sim_size,
#                                               cores=cores,
#                                               white_adjust='hc0',
#                                               seed=seed,
#                                               mincer_spec=mincer_spec,
#                                               execute_additional_tsts=TRUE)
# result_power_3_fix_hc0_matrix <- create_result_matrix(result_power_3_fix_hc0)
# write_results_to_txt(name = 'result_power_3_fix_hc0', txt_file = result_txt_file)









# Size loop with estimated parameters
# result_size_est <- estimation_loop_par(n_loop=n_loop_est,
#                                        est_window=est_window,
#                                        oos_window=oos_window_est,
#                                        tolerance_lvl=tolerance_lvl,
#                                        var_spec_sim=var_spec_sim_size,
#                                        mean_spec_sim=mean_spec_sim_size,
#                                        dist_spec_sim=dist_spec_sim_size,
#                                        fixed_pars_sim=fixed_pars_sim_size,
#                                        estimate=TRUE,
#                                        var_spec_est=var_spec_est_size,
#                                        mean_spec_est=mean_spec_est_size,
#                                        dist_spec_est=dist_spec_est_size,
#                                        fixed_pars_est=NA,
#                                        cores=cores,
#                                        white_adjust=white_adjust,
#                                        seed=seed,
#                                        mincer_spec=mincer_spec,
#                                        execute_additional_tsts=FALSE)
# result_size_est_matrix <- create_result_matrix(result_size_est)

# Power loop 1 with estimated parameters
# result_power_1_est <- estimation_loop_par(n_loop=n_loop_est,
#                                           est_window=est_window,
#                                           oos_window=oos_window_est,
#                                           tolerance_lvl=tolerance_lvl,
#                                           var_spec_sim=var_spec_sim_power_1,
#                                           mean_spec_sim=mean_spec_sim_power_1,
#                                           dist_spec_sim=dist_spec_sim_power_1,
#                                           fixed_pars_sim=fixed_pars_sim_power_1,
#                                           estimate=TRUE,
#                                           var_spec_est=var_spec_est_power_1,
#                                           mean_spec_est=mean_spec_est_power_1,
#                                           dist_spec_est=dist_spec_est_power_1,
#                                           fixed_pars_est=NA,
#                                           cores=cores,
#                                           white_adjust=white_adjust,
#                                           seed=seed,
#                                           mincer_spec=mincer_spec,
#                                           execute_additional_tsts=FALSE)
# result_power_1_est_matrix <- create_result_matrix(result_power_1_est)

# Power loop 2 with estimated parameters
# result_power_2_est <- estimation_loop_par(n_loop=n_loop_est,
#                                           est_window=est_window,
#                                           oos_window=oos_window_est,
#                                           tolerance_lvl=tolerance_lvl,
#                                           var_spec_sim=var_spec_sim_power_2,
#                                           mean_spec_sim=mean_spec_sim_power_2,
#                                           dist_spec_sim=dist_spec_sim_power_2,
#                                           fixed_pars_sim=fixed_pars_sim_power_2,
#                                           estimate=TRUE,
#                                           var_spec_est=var_spec_est_power_2,
#                                           mean_spec_est=mean_spec_est_power_2,
#                                           dist_spec_est=dist_spec_est_power_2,
#                                           fixed_pars_est=NA,
#                                           cores=cores,
#                                           white_adjust=white_adjust,
#                                           seed=seed,
#                                           mincer_spec=mincer_spec,
#                                           execute_additional_tsts=FALSE)
# result_power_2_est_matrix <- create_result_matrix(result_power_2_est)