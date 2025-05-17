####################################################################
################### Function calculates VaR and ES #################
####################################################################
VaR_ES_forecast <- function(data_zoo,
                            var_spec,
                            mean_spec,
                            dist_spec,
                            tolerance_lvl,
                            estimate,
                            fixed_pars_est){
  
  # Store last date
  last_date <- tail(index(data_zoo), 1)
  
  if(estimate){
    # GARCH Specification
    garchspec <- ugarchspec(variance.model = var_spec,
                            mean.model = mean_spec,
                            distribution.model = dist_spec)
    
    # Fitting GARCH model
    garchfit <- tryCatch(
      {
        ugarchfit(spec = garchspec,
                  data = data_zoo,
                  solver = 'hybrid')
      }, error = function(e){
        return(NA)
      }
    )
    
  } else {
    # GARCH Specification with all parameters fixed
    garchspec <- ugarchspec(variance.model = var_spec,
                            mean.model = mean_spec,
                            distribution.model = dist_spec,
                            fixed.pars = fixed_pars_est)
    
    # Filter data
    garchfit <- ugarchfilter(spec = garchspec,
                             data = data_zoo)
  }
  
  if(suppressWarnings(!is.na(garchfit))){
    
    # Extracting coefficients
    coef_fit <- coef(garchfit)
    
    
    #TEST
    # mu <<- c(mu, coef_fit['mu'])
    # omega <<- c(omega, coef_fit['omega'])
    # alpha1 <<- c(alpha1, coef_fit['alpha1'])
    # beta1 <<- c(beta1, coef_fit['beta1'])
    
    
    # Extracting last standard deviation from fit
    st_dev_t_min_1 <- tail(sigma(garchfit), 1)

    # Extracting last residual from fit
    residual_t_min_1 <- tail(residuals(garchfit), 1)
    
    # Extracting skewness, shape and lambda parameter (if not prevalent, NA gets assigned)
    skew <- ifelse('skew' %in% names(coef_fit),
                   coef_fit['skew'],
                   NA)
    shape <- ifelse('shape' %in% names(coef_fit),
                    coef_fit['shape'],
                    NA)
    lambda <- ifelse('ghlambda' %in% names(coef_fit),
                     coef_fit['ghlambda'],
                     NA)
    
    # Forecasting
    if(estimate){
      garchforecast <- ugarchforecast(fitORspec = garchfit,
                                      n.ahead = 1)
    } else {
      garchforecast <- ugarchforecast(fitORspec = garchspec,
                                      data = data_zoo,
                                      n.ahead = 1)
    }
    
    # Extracting mu and sigma
    mu <- fitted(garchforecast)
    sigma <- sigma(garchforecast)
    
    # Function returns pth quantile of current distribution
    pth_quantile <- function(p) {
      q <- qdist(distribution = dist_spec,
                 p = p,
                 skew = skew,
                 shape = shape,
                 lambda = lambda)
      return(q)
    }
    
    # Calculate VaR
    VaR <- mu + sigma * pth_quantile(p = tolerance_lvl)
    
    # Calculate ES
    integrated <- integrate(pth_quantile,
                            lower = 0,
                            upper = tolerance_lvl)
    integrated_value <- integrated[['value']]
    ES <- mu + sigma / tolerance_lvl * integrated_value
    
    # Return results
    results <- data.frame(last_date = last_date,
                          mu = as.double(mu),
                          sigma = as.double(sigma),
                          variance = as.double(sigma^2),
                          st_dev_t_min_1 = as.double(st_dev_t_min_1),
                          variance_t_min_1 = as.double(st_dev_t_min_1)^2,
                          residual_t_min_1 = as.double(residual_t_min_1),
                          residual_t_min_1_quadr = as.double(residual_t_min_1)^2,
                          VaR = as.double(VaR),
                          ES = as.double(ES))
  } else {
    results <- data.frame(last_date = last_date,
                          mu=NA,
                          sigma = NA,
                          variance = NA,
                          st_dev_t_min_1 = NA,
                          variance_t_min_1 = NA,
                          residual_t_min_1 = NA,
                          residual_t_min_1_quadr = NA,
                          VaR = NA,
                          ES = NA)
  }
  
  results_zoo <- zoo(results[, -1], order.by = results$last_date)
  return(results_zoo)
}

#############################################################################
###################     Mincer Regression    ################################
#############################################################################
mincer_regression <- function(formula,
                              shortfall_df,
                              h0,
                              white_adjust){
  mincer_reg <- lm(formula = formula, data = shortfall_df)
  mincer_reg_result <- linearHypothesis(mincer_reg, h0, white.adjust = white_adjust)
  p_value <- mincer_reg_result$'Pr(>F)'[2]
  return(p_value)
}

#############################################################################
################### Parallel estimation Loop ################################
#############################################################################
estimation_loop_par <- function(n_loop,
                                est_window,
                                oos_window,
                                tolerance_lvl,
                                var_spec_sim,
                                mean_spec_sim,
                                dist_spec_sim,
                                fixed_pars_sim,
                                estimate,
                                var_spec_est,
                                mean_spec_est,
                                dist_spec_est,
                                fixed_pars_est,
                                cores,
                                white_adjust,
                                seed=NA,
                                mincer_spec){
  
  # Register cores
  registerDoParallel(cores = cores)
  
  # Estimation loop
  result_foreach <- foreach(i=1:n_loop) %dopar% {
    
    # Set seed for reproducability of results
    if(!is.na(seed)){
      set.seed(seed+i^2-i)
    }
    
    # Only displayed if cores=1
    if(cores==1){
      cat(paste0('Iteration Nr.: ', i, '\n'))
    }
    
    # Simulate garch process
    garchspec_sim <- ugarchspec(variance.model = var_spec_sim,
                                mean.model = mean_spec_sim,
                                distribution.model = dist_spec_sim,
                                fixed.pars = fixed_pars_sim)
    
    garchsimulation <- ugarchpath(spec = garchspec_sim,
                                  n.sim = est_window + oos_window,
                                  n.start = 1000,
                                  m.sim = 1)
    
    # Create zoo object of simulation
    garchsimulation_df <- as.data.frame(garchsimulation@path$seriesSim)
    garchsimulation_df <- garchsimulation_df %>%
      mutate(Date = as.Date('2000-01-01') + 1:nrow(garchsimulation_df)) %>%
      rename(Return = V1)
    garchsimulation_zoo <- zoo(garchsimulation_df$Return, garchsimulation_df$Date)
    
    # Calculate GARCH model
    rolling_VaR_ES <- rollapply(garchsimulation_zoo,
                                width = est_window,
                                FUN = function(x) VaR_ES_forecast(data_zoo = x,
                                                                  var_spec = var_spec_est,
                                                                  mean_spec = mean_spec_est,
                                                                  dist_spec = dist_spec_est,
                                                                  tolerance_lvl = tolerance_lvl,
                                                                  estimate = estimate,
                                                                  fixed_pars_est = fixed_pars_est),
                                align = 'right',
                                coredata = FALSE)

    rolling_VaR_ES_lead <- stats::lag(rolling_VaR_ES, k = -1)

    rolling_VaR_ES_df <- as.data.frame(rolling_VaR_ES_lead)
    rolling_VaR_ES_df <- rolling_VaR_ES_df %>%
      mutate(Date = as.Date(rownames(rolling_VaR_ES_df)))
    
    VaR_ES_results_df <- inner_join(garchsimulation_df, rolling_VaR_ES_df, by = 'Date')

    # Filter for shortfalls
    shortfall_df <- VaR_ES_results_df %>%
      filter(Return <= VaR) %>%
      mutate(shortfall = Return - ES) %>%
      mutate(residual_t_min_1_quadr_lower_0 = ifelse(residual_t_min_1 < 0, residual_t_min_1_quadr, 0)) %>%
      mutate(indicator_residual_lower_0 = ifelse(residual_t_min_1 < 0, 1, 0))
    
    # Save number of observations that go into Mincer regression
    n_obs_mincer <- nrow(shortfall_df)

    # Mincer regression execution
    p_values_vec <- vector()
    for(j in 1:length(mincer_spec)){
      p_values_vec[j] <- mincer_regression(formula = mincer_spec[[j]][['formula']],
                                           shortfall_df = shortfall_df,
                                           h0 = mincer_spec[[j]][['h0']],
                                           white_adjust = white_adjust)
    }

    # Save results
    result_lst <- list()
    result_lst_names <- names(mincer_spec)

    for(j in 1:length(mincer_spec)){
      result_lst[[result_lst_names[j]]][['p']] <- p_values_vec[j]
      result_lst[[result_lst_names[j]]][['p0_01']] <- ifelse(p_values_vec[j] < 0.01, 1, 0)
      result_lst[[result_lst_names[j]]][['p0_05']] <- ifelse(p_values_vec[j] < 0.05, 1, 0)
      result_lst[[result_lst_names[j]]][['p0_1']] <- ifelse(p_values_vec[j] < 0.1, 1, 0)
    }
    result_lst[['n_obs_mincer']] <- n_obs_mincer
    result_lst
  }

  # Organize results
  result <- list()
  for(mincer_reg_name in names(mincer_spec)){
    result[[mincer_reg_name]] <- list(p = vector(),
                                      p0_01 = vector(),
                                      p0_05 = vector(),
                                      p0_1 = vector())
  }
  result[['n_obs_mincer']] <- vector()
  
  for(i in 1:n_loop){
    for(method in names(result_foreach[[i]])){
      if(method == 'n_obs_mincer'){
        result[[method]] <- c(result[[method]], result_foreach[[i]][[method]])
      } else {
        for(metrix in names(result_foreach[[i]][[method]])){
          result[[method]][[metrix]] <- c(result[[method]][[metrix]], result_foreach[[i]][[method]][[metrix]])
        }
      }
    }
  }
  return(result)
}

#############################################################################
################### Create result matrix     ################################
#############################################################################
create_result_matrix <- function(result_lst){
  result_lst[['n_obs_mincer']] <- NULL
  row_names <- names(result_lst)
  col_names <- names(result_lst[[1]])[-1]
  
  result_matrix <- matrix(nrow = length(row_names), ncol = length(col_names))
  rownames(result_matrix) <- row_names
  colnames(result_matrix) <- col_names
  
  for(row in row_names){
    for(col in col_names){
      result_matrix[row, col] <- round(mean(result_lst[[row]][[col]]), digits = 3)
    }
  }
  return(result_matrix)
}