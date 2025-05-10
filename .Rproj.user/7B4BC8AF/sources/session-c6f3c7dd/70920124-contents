####################################################################
################### Function calculates VaR and ES #################
####################################################################
VaR_ES_forecast <- function(data_zoo,
                            var_spec,
                            mean_spec,
                            dist_spec,
                            tolerance_lvl,
                            estimate,
                            fixed_pars){
  
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
                            fixed.pars = fixed_pars)
    
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
    
    
    # Extracting variance from fit
    variance_t_min_1 <- tail(sigma(garchfit), 1) ^ 2
    
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
                          mu=as.double(mu),
                          sigma = as.double(sigma),
                          variance = as.double(sigma^2),
                          variance_t_min_1 = as.double(variance_t_min_1),
                          residual_t_min_1 = as.double(residual_t_min_1),
                          residual_t_min_1_quadr = as.double(residual_t_min_1)^2,
                          VaR = as.double(VaR),
                          ES = as.double(ES))
  } else {
    results <- data.frame(last_date = last_date,
                          mu=NA,
                          sigma = NA,
                          variance = NA,
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
                                estimate,
                                fixed_pars,
                                var_spec_sim,
                                mean_spec_sim,
                                dist_spec_sim,
                                omega,
                                alpha1,
                                beta1,
                                eta11,
                                var_spec_est,
                                mean_spec_est,
                                dist_spec_est,
                                cores,
                                white_adjust,
                                seed){
  
  # Register cores
  registerDoParallel(cores = cores)
  
  # Estimation loop
  result_foreach <- foreach(i=1:n_loop) %dopar% {
    
    # Set seed for reproducability of results
    set.seed(seed+i)
    
    # Only displayed if cores=1
    if(cores==1){
      cat(paste0('Iteration Nr.: ', i, '\n'))
    }
    
    # Simulate garch process
    garchspec_sim <- ugarchspec(variance.model = var_spec_sim,
                                mean.model = mean_spec_sim,
                                distribution.model = dist_spec_sim,
                                fixed.pars = list(mu = 0,
                                                  omega = omega,
                                                  alpha1 = alpha1,
                                                  beta1 = beta1,
                                                  eta11 = eta11))
    
    garchsimulation <- ugarchpath(spec = garchspec_sim,
                                  n.sim = est_window + oos_window,
                                  n.start = 0,
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
                                                                  fixed_pars = fixed_pars),
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
      mutate(residual_t_min_1_quadr_lower_0 = ifelse(residual_t_min_1 < 0, residual_t_min_1_quadr, 0))

    # Mincer regression specifications
    mincer_spec <- list(formula = list(simple_shortfall = shortfall ~ 1,
                                       simple_return = Return ~ ES,
                                       variance_shortfall = shortfall ~ variance_t_min_1,
                                       variance_return = Return ~ variance_t_min_1 + ES,
                                       residual_sqrt_shortfall = shortfall ~ residual_t_min_1_quadr,
                                       residual_sqrt_return = Return ~ residual_t_min_1_quadr + ES,
                                       residual_sqrt0_shortfall = shortfall ~ residual_t_min_1_quadr_lower_0,
                                       residual_sqrt0_return = Return ~ residual_t_min_1_quadr_lower_0 + ES,
                                       full_shortfall = shortfall ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0,
                                       full_return = Return ~ variance_t_min_1 + residual_t_min_1_quadr + residual_t_min_1_quadr_lower_0 + ES),
                        h0 = list(simple_shortfall = c('(Intercept) = 0'),
                                  simple_return = c('(Intercept) = 0', 'ES = 1'),
                                  variance_shortfall = c('(Intercept) = 0', 'variance_t_min_1 = 0'),
                                  variance_return = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'ES = 1'),
                                  residual_sqrt_shortfall = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0'),
                                  residual_sqrt_return = c('(Intercept) = 0', 'residual_t_min_1_quadr = 0', 'ES = 1'),
                                  residual_sqrt0_shortfall = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0'),
                                  residual_sqrt0_return = c('(Intercept) = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1'),
                                  full_shortfall = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0'),
                                  full_return = c('(Intercept) = 0', 'variance_t_min_1 = 0', 'residual_t_min_1_quadr = 0', 'residual_t_min_1_quadr_lower_0 = 0', 'ES = 1')))
    
    # Mincer regression execution
    for(j in 1:length(mincer_spec[['formula']])){
      p_value_j <- mincer_regression(formula = mincer_spec[['formula']][[j]],
                                     shortfall_df = shortfall_df,
                                     h0 = mincer_spec[['h0']][[j]],
                                     white_adjust = white_adjust)
      assign(paste0('p_value_', j), p_value_j)
      rm(p_value_j)
    }
    # (Maybe add coefficients in this way: coef_var <- as.double(coef(mincer_reg)['variance_t_min_1']))
    
    # Save results
    result_lst <- list()
    result_lst_names <- c('simple_shortfall', 'simple_return',
                          'variance_shortfall', 'variance_return',
                          'residual_sqrt_shortfall', 'residual_sqrt_return',
                          'residual_sqrt0_shortfall', 'residual_sqrt0_return',
                          'full_shortfall', 'full_return')
    for(j in 1:10){
      p_value <- get(paste0('p_value_', j))
      result_lst[[result_lst_names[j]]][['p']] <- p_value
      result_lst[[result_lst_names[j]]][['p0_01']] <- ifelse(p_value <= 0.01, 1, 0)
      result_lst[[result_lst_names[j]]][['p0_05']] <- ifelse(p_value <= 0.05, 1, 0)
      result_lst[[result_lst_names[j]]][['p0_1']] <- ifelse(p_value <= 0.1, 1, 0)
    }
    result_lst
  }

  # Organize results
  result <- list(simple_shortfall = list(p = vector(),
                                         p0_01 = vector(),
                                         p0_05 = vector(),
                                         p0_1 = vector()),
                 simple_return = list(p = vector(),
                                      p0_01 = vector(),
                                      p0_05 = vector(),
                                      p0_1 = vector()),
                 variance_shortfall = list(p = vector(),
                                           p0_01 = vector(),
                                           p0_05 = vector(),
                                           p0_1 = vector()),
                 variance_return = list(p = vector(),
                                        p0_01 = vector(),
                                        p0_05 = vector(),
                                        p0_1 = vector()),
                 residual_sqrt_shortfall = list(p = vector(),
                                                p0_01 = vector(),
                                                p0_05 = vector(),
                                                p0_1 = vector()),
                 residual_sqrt_return = list(p = vector(),
                                             p0_01 = vector(),
                                             p0_05 = vector(),
                                             p0_1 = vector()),
                 residual_sqrt0_shortfall = list(p = vector(),
                                                 p0_01 = vector(),
                                                 p0_05 = vector(),
                                                 p0_1 = vector()),
                 residual_sqrt0_return = list(p = vector(),
                                              p0_01 = vector(),
                                              p0_05 = vector(),
                                              p0_1 = vector()),
                 full_shortfall = list(p = vector(),
                                       p0_01 = vector(),
                                       p0_05 = vector(),
                                       p0_1 = vector()),
                 full_return = list(p = vector(),
                                    p0_01 = vector(),
                                    p0_05 = vector(),
                                    p0_1 = vector()))
  for(i in 1:n_loop){
    for(method in names(result_foreach[[i]])){
      for(metrix in names(result_foreach[[i]][[method]])){
        result[[method]][[metrix]] <- c(result[[method]][[metrix]], result_foreach[[i]][[method]][[metrix]])
      }
    }
  }
  return(result)
}