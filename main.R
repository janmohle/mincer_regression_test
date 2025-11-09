# Initial cleaning
rm(list = ls())
if (dev.cur() != 1) {
  dev.off()
}
cat('\14')

# Load libraries
#renv::install('tibble')
#renv::install('dplyr')
#renv::install('purrr')
#renv::install('rugarch')
#renv::install('zoo')
#renv::install('foreach')
#renv::install('doParallel')

library(tibble)
library(dplyr)
library(purrr)
library(rugarch)
library(zoo)
library(foreach)
library(doParallel)

source('functions.R')

# Create process folders if they dont exist
if (!dir.exists('run_logs')){
  dir.create('run_logs',
             recursive = TRUE)
}
if (!dir.exists('results')){
  dir.create('results',
             recursive = TRUE)
}

#############################################
########### Global parameters ###############
#############################################
result_txt_file = 'results/results.txt'
result_latex_file = 'results/results_latex.txt'
tolerance_lvl = 0.05
n_loop_est = 1000
n_loop_fix = 1000
oos_window_est = 5000
oos_window_fix = 1000#0
est_window = 750#2000
cores = 100#detectCores()
seed = 78
sim = 0
variant = 1
white_adjust = c('hc3', FALSE)
execute_additional_tsts = TRUE
lags_es_cc = 5
n_boot_de = 2000
empirical = TRUE
estimation = FALSE


########################################################
########### Simulation specifications ##################
########################################################

# Sim 0
var_spec_sim_0 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_0 = list(armaOrder = c(0,0))
dist_spec_sim_0 = 'norm'
fixed_pars_sim_0 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.05,
                        beta1 = 0.94)

# Sim 1
var_spec_sim_1 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_1 = list(armaOrder = c(0,0))
dist_spec_sim_1 = 'norm'
fixed_pars_sim_1 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.03,
                        beta1 = 0.94,
                        gamma1 = 0.04)

# Sim 2
var_spec_sim_2 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_2 = list(armaOrder = c(0,0))
dist_spec_sim_2 = 'norm'
fixed_pars_sim_2 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.01,
                        beta1 = 0.94,
                        gamma1 = 0.08)

# Sim 3
var_spec_sim_3 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_3 = list(armaOrder = c(0,0))
dist_spec_sim_3 = 'std'
fixed_pars_sim_3 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.05,
                        beta1 = 0.94,
                        shape = 4)

# Sim 4
var_spec_sim_4 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_4 = list(armaOrder = c(0,0))
dist_spec_sim_4 = 'std'
fixed_pars_sim_4 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.03,
                        beta1 = 0.94,
                        gamma1 = 0.04,
                        shape = 4)

# Sim 5
var_spec_sim_5 = list(model = 'gjrGARCH', garchOrder = c(1, 1))
mean_spec_sim_5 = list(armaOrder = c(0,0))
dist_spec_sim_5 = 'std'
fixed_pars_sim_5 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.01,
                        beta1 = 0.94,
                        gamma1 = 0.08,
                        shape = 4)

# Sim 6
var_spec_sim_6 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_6 = list(armaOrder = c(0,0))
dist_spec_sim_6 = 'std'
fixed_pars_sim_6 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.05,
                        beta1 = 0.94,
                        shape = 8)

# Sim 7
var_spec_sim_7 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_sim_7 = list(armaOrder = c(0,0))
dist_spec_sim_7 = 'std'
fixed_pars_sim_7 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.05,
                        beta1 = 0.94,
                        shape = 12)


########################################################
########### Estimation specifications ##################
########################################################

# Variant 1 (also for cases without parameter estimation -> fixed parameters are specified)
var_spec_est_1 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_1 = list(armaOrder = c(0,0))
dist_spec_est_1 = 'norm'
fixed_pars_est_1 = list(mu = 0,
                        omega = 0.005,
                        alpha1 = 0.05,
                        beta1 = 0.94)

# Variant 2
var_spec_est_2 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_est_2 = list(armaOrder = c(0,0))
dist_spec_est_2 = 'norm'
fixed_pars_est_2 = NA

# Variant 3
var_spec_est_3 = list(model = 'sGARCH',garchOrder = c(1,1))
mean_spec_est_3 = list(armaOrder = c(0,0))
dist_spec_est_3 = 'std'
fixed_pars_est_3 = NA

# Variant 4
var_spec_est_4 = list(model = 'gjrGARCH',garchOrder = c(1,1))
mean_spec_est_4 = list(armaOrder = c(0,0))
dist_spec_est_4 = 'std'
fixed_pars_est_4 = NA


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
write(paste0('\n\nExecution started: ', Sys.time(), '\n\n'),
      file = result_txt_file,
      append = TRUE)

# Write date and time of execution start into latex results file
write(paste0('\n\nExecution started: ', Sys.time(), '\n\n'),
      file = result_latex_file,
      append = TRUE)

# Initiate progress_information.txt
write(paste0('Execution started: ', Sys.time(), '\n\n'),
      file = 'run_logs/progress_information.txt',
      append = FALSE)

# Initiate solver_problems.txt
write('Solver problems (Check functions.R for location of problems encoded in numbers):\n\n',
      file = 'run_logs/solver_problems.txt',
      append = FALSE)

# Initiate NA_information.txt
write(paste0('Execution started: ', Sys.time(), '\n\n'),
      file = 'run_logs/NA_information.txt',
      append = FALSE)

# Loop executions
for(i in sim){
  
  # Progress information: Current simulation
  write(paste0('\nSimulation ', i, ':\n'),
        file = 'run_logs/progress_information.txt',
        append = TRUE)
  
  # Current simulation for solver_problems.txt
  write(paste0('\nSimulation ', i, ':\n'),
        file = 'run_logs/solver_problems.txt',
        append = TRUE)
  
  # Current simulation for NA_information.txt
  write(paste0('\nSimulation ', i, ':\n'),
        file = 'run_logs/NA_information.txt',
        append = TRUE)
  
  result_lst <- estimation_loop_par(n_loop = ifelse(estimation, n_loop_est, n_loop_fix),
                                    est_window = est_window,
                                    oos_window = ifelse(estimation, oos_window_est, oos_window_fix),
                                    tolerance_lvl = tolerance_lvl,
                                    var_spec_sim = get(paste0('var_spec_sim_', i)),
                                    mean_spec_sim = get(paste0('mean_spec_sim_', i)),
                                    dist_spec_sim = get(paste0('dist_spec_sim_', i)),
                                    fixed_pars_sim = get(paste0('fixed_pars_sim_', i)),
                                    estimate = estimation,
                                    par_corr = estimation,
                                    var_spec_est = get(paste0('var_spec_est_', variant)),
                                    mean_spec_est = get(paste0('mean_spec_est_', variant)),
                                    dist_spec_est = get(paste0('dist_spec_est_', variant)),
                                    fixed_pars_est = get(paste0('fixed_pars_est_', variant)),
                                    cores = cores,
                                    white_adjust = white_adjust,
                                    seed = seed,
                                    mincer_spec = mincer_spec,
                                    execute_additional_tsts = execute_additional_tsts,
                                    lags_es_cc = lags_es_cc,
                                    n_boot_de = n_boot_de,
                                    empirical = empirical)
  
  # Number of NAs in UC_par_corr and CC_par_corr printed to NA_information.txt
  if(estimation){
    write(paste0('Number of NAs in UC_par_corr: ', sum(is.na(result_lst[['UC_par_corr']][[white_adjust[1]]])), '\n'),
          file = 'run_logs/NA_information.txt',
          append = TRUE)
    write(paste0('Number of NAs in CC_par_corr: ', sum(is.na(result_lst[['CC_par_corr']][[white_adjust[1]]])), '\n'),
          file = 'run_logs/NA_information.txt',
          append = TRUE)
  }

  result_matrix <- create_result_matrix(result_lst = result_lst,
                                        na_rm = TRUE)
  
  assign(paste0('result_sim_', i, '_var_', variant, ifelse(estimation, '_est_', '_fix_'), ifelse(empirical, 'emp', 'par')), result_lst)
  assign(paste0('result_sim_', i, '_var_', variant, ifelse(estimation, '_est_', '_fix_'), ifelse(empirical, 'emp', 'par'), '_matrix'), result_matrix)
  
  write_results_to_txt(name = paste0('result_sim_', i, '_var_', variant, ifelse(estimation, '_est_', '_fix_'), ifelse(empirical, 'emp', 'par')),
                       txt_file = result_txt_file)
  write_results_to_latex(name = paste0('result_sim_', i, '_var_', variant, ifelse(estimation, '_est_', '_fix_'), ifelse(empirical, 'emp', 'par')),
                         txt_file = result_latex_file,
                         double_hline = ifelse(execute_additional_tsts, length(mincer_spec), NA))
  
  saveRDS(result_lst,
          file = paste0('results/result_sim_', i, '_var_', variant, ifelse(estimation, '_est_', '_fix_'), ifelse(empirical, 'emp', 'par'), '.RData'))
}
rm(result_lst, result_matrix, i)