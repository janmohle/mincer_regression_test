####################################################################
################### Function calculates VaR and ES #################
####################################################################
VaR_ES_forecast <- function(data_zoo,
                            var_spec,
                            mean_spec,
                            dist_spec,
                            tolerance_lvl,
                            estimate,
                            par_corr,
                            fixed_pars_est,
                            empirical){
  
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
    
    # Extract covariance matrix of parameters and calculate gradients of mu and sigma evaluated at coef_fit
    if(par_corr){
      vcov_matrix_par <- matrix_to_string(vcov(garchfit))
      
      mu_sigma_func <- function(par_vec){
        
        spec_tmp <- ugarchspec(variance.model = var_spec,
                               mean.model = mean_spec,
                               distribution.model = dist_spec,
                               fixed.pars = as.list(par_vec))
        
        forecast_tmp <- ugarchforecast(fitORspec = spec_tmp,
                                       data = data_zoo,
                                       n.ahead = 1)
        
        return(c(mu = as.double(fitted(forecast_tmp)), sigma = as.double(sigma(forecast_tmp))))
      }
      gradients <- jacobian(func = mu_sigma_func,
                            x = coef_fit,
                            method = 'Richardson')
      mu_grad <- vector_to_string(gradients[1,])
      sigma_grad <- vector_to_string(gradients[2,])
    }
    
    # Extracting mu and sigma
    mu <- fitted(garchforecast)
    sigma <- sigma(garchforecast)
    
    if(empirical){
      # Empirical distribution
      emp_dist <- as.vector(residuals(garchfit, standardize = TRUE))
      emp_dist_vec <- vector_to_string(emp_dist)
      #TEST!
      #plot(density(emp_dist))
      # pth quantile of empirical distribution
      q_emp_dist <- quantile(emp_dist,
                             probs = tolerance_lvl,
                             na.rm = TRUE,
                             names = FALSE,
                             type = 8) # Approximately median unbiased regardless of distribution
      #TEST!
      #q_empiricl <<- c(q_empiricl,q_emp_dist)
      # Calculating standardized ES based on empirical distribution
      standardized_ES <- mean(emp_dist[emp_dist <= q_emp_dist])
      #TEST!
      #st_ES <<- c(st_ES, standardized_ES)
      # VaR and ES calculation
      VaR <- mu + sigma * q_emp_dist
      ES <- mu + sigma * standardized_ES
      
    } else{
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
    }
    
    # Return results
    results <- data.frame(last_date = last_date,
                          mu = as.double(mu),
                          sigma = as.double(sigma),
                          variance = as.double(sigma)^2,
                          skew = skew,
                          shape = shape,
                          lambda = lambda,
                          dist_spec = dist_spec,
                          st_dev_t_min_1 = as.double(st_dev_t_min_1),
                          variance_t_min_1 = as.double(st_dev_t_min_1)^2,
                          residual_t_min_1 = as.double(residual_t_min_1),
                          residual_t_min_1_quadr = as.double(residual_t_min_1)^2,
                          VaR = as.double(VaR),
                          ES = as.double(ES),
                          vcov_matrix_par = ifelse(par_corr,vcov_matrix_par,NA),
                          mu_grad = ifelse(par_corr,mu_grad,NA),
                          sigma_grad = ifelse(par_corr,sigma_grad,NA),
                          emp_dist_vec = ifelse(empirical, emp_dist_vec, NA))
  } else {
    results <- data.frame(last_date = last_date,
                          mu = NA,
                          sigma = NA,
                          variance = NA,
                          skew = NA,
                          shape = NA,
                          lambda = NA,
                          dist_spec = NA,
                          st_dev_t_min_1 = NA,
                          variance_t_min_1 = NA,
                          residual_t_min_1 = NA,
                          residual_t_min_1_quadr = NA,
                          VaR = NA,
                          ES = NA,
                          vcov_matrix_par = NA,
                          mu_grad = NA,
                          sigma_grad = NA,
                          emp_dist_vec = NA)
  }
  
  results_zoo <- zoo(dplyr::select(results, -'last_date'), order.by = results[['last_date']])
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

##########################################################################
############# Unconditional coverage test for ES #########################
##########################################################################
ES_uc_backtest <- function(CumVio,
                           tolerance_lvl,
                           par_corr,
                           est_window,
                           Return,
                           mu,
                           sigma,
                           dist_spec,
                           skew,
                           shape,
                           lambda,
                           vcov_matrix_par_code,
                           mu_grad_code,
                           sigma_grad_code){
  
  # Extracting non-NA cumulative violations
  not_na_CumVio <- !is.na(CumVio)
  H_hut <- CumVio[not_na_CumVio]
  
  # Mean and length of H_hut
  H_mean <- mean(H_hut)
  n <- length(H_hut)
  
  # Parameter uncertainty correction
  if(par_corr){
    
    # Create tibble with required data
    tibble_tmp <- tibble(Return = Return,
                         mu = mu,
                         sigma = sigma,
                         dist_spec = dist_spec,
                         skew = skew,
                         shape = shape,
                         lambda = lambda,
                         vcov_matrix_par_code = vcov_matrix_par_code,
                         mu_grad_code = mu_grad_code,
                         sigma_grad_code = sigma_grad_code) %>%
      mutate(eps = (Return - mu) / sigma,
             vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
             mu_grad         = map(mu_grad_code, eval_string_code),
             sigma_grad      = map(sigma_grad_code, eval_string_code))
    
    # Calculate R_ES
    tibble_tmp <- tibble_tmp %>%
      mutate(R_ES_t = pmap(
        list(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma),
        function(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma) {
          
          # Check for any NA values
          if (anyNA(c(eps, dist_spec, mu_grad, sigma_grad, sigma))) {return(NA_real_)}
          
          # Compute the density of eps
          density_eps <- ddist(distribution = dist_spec,
                               y = eps,
                               skew = skew,
                               shape = shape,
                               lambda = lambda)
          
          # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
          indicator <- ifelse(
            eps <= qdist(distribution = dist_spec,
                         p = tolerance_lvl,
                         skew = skew,
                         shape = shape,
                         lambda = lambda),
            1, 0)
          
          # Compute gradient expression
          grad_term <- (mu_grad + eps * sigma_grad) / sigma
          
          # Final result
          return(density_eps * indicator * grad_term)
        }
      ))
    
    # Sum R_ES_t vectors across all rows (ignoring NA entries)
    R_ES_sum <- tibble_tmp %>%
      filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      pull(R_ES_t) %>%
      reduce(`+`)
    
    # Number of R_ES_t without NA
    n_no_NA <- tibble_tmp %>%
      filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      nrow()
    
    # Calculate final R_ES
    R_ES <- (1 / (tolerance_lvl * n_no_NA)) * R_ES_sum
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
    # Final calculation of correction factor
    corr_fac <- (n/est_window) * t(as.vector(R_ES)) %*% as.matrix(W_mean) %*% as.vector(R_ES)
    
    ##TEST
    #tibble_uc <<- tibble_tmp
    
  } else {
    
    # No correction if estimate = FALSE
    corr_fac <- 0
  }
  
  ##TEST
  #corr_fac_save <<- corr_fac
  
  # Test statistic calculation
  U <- (sqrt(n) * (H_mean - (tolerance_lvl / 2))) / sqrt(tolerance_lvl * ((1 / 3) - (tolerance_lvl / 4)) + corr_fac)
  
  # p-value calculation
  p <- 2 * pnorm(q = abs(U),
                 lower.tail = FALSE)

  # Return p-value
  return(as.numeric(p))
}

##########################################################################
############# Conditional coverage test for ES ###########################
##########################################################################
ES_cc_backtest <- function(CumVio,
                           tolerance_lvl,
                           lags,                           
                           par_corr,
                           est_window,
                           Return,
                           mu,
                           sigma,
                           dist_spec,
                           skew,
                           shape,
                           lambda,
                           vcov_matrix_par_code,
                           mu_grad_code,
                           sigma_grad_code){
  
  # Extracting non-NA cumulative violations
  not_na_CumVio <- !is.na(CumVio)
  H_hut <- CumVio[not_na_CumVio]
  
  # Length of H_hut
  n <- length(H_hut)
  
  # Variance of H_huts
  gamma_n0 <- 1 / n * sum((H_hut - (tolerance_lvl / 2)) * (H_hut - (tolerance_lvl / 2)))
  
  # Vector with covariance between H_hut and j-laged H_hut
  gamma_nj <- vector(length = lags)
  
  for(j in 1:lags){
    gamma_nj[j] <- 1 / (n - j) * sum((H_hut[(j+1):n] - (tolerance_lvl / 2)) * (H_hut[1:(n-j)] - (tolerance_lvl / 2)))
  }
  
  # Vector with correlations between H_hut and j-laged H_hut
  rho_nj <- as.vector(gamma_nj / gamma_n0)
  
  # Parameter uncertainty correction
  if(par_corr){
    
    # Create tibble with required data
    tibble_tmp <- tibble(Return = Return,
                         mu = mu,
                         sigma = sigma,
                         dist_spec = dist_spec,
                         skew = skew,
                         shape = shape,
                         lambda = lambda,
                         vcov_matrix_par_code = vcov_matrix_par_code,
                         mu_grad_code = mu_grad_code,
                         sigma_grad_code = sigma_grad_code,
                         CumVio = CumVio) %>%
      mutate(eps = (Return - mu) / sigma,
             vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
             mu_grad         = map(mu_grad_code, eval_string_code),
             sigma_grad      = map(sigma_grad_code, eval_string_code))
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
    # Compute R_j columns
    R_j_list <- list()
    for(j in 1:lags){
      
      # Create lagged CumVio column
      CumVio_lag_j <- paste0('CumVio_lag_', j)
      tibble_tmp <- tibble_tmp %>%
        mutate(!!CumVio_lag_j := dplyr::lag(x = CumVio, n = j))
      
      # Create R_j_t column
      R_j_t <- paste0('R_', j, '_t')
      tibble_tmp <- tibble_tmp %>%
        mutate(!!R_j_t := pmap(
          list(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma, .data[[CumVio_lag_j]]),
          function(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma, CumVio_lag_j) {
            
            # Check for any NA values
            if (anyNA(c(eps, dist_spec, mu_grad, sigma_grad, sigma, CumVio_lag_j))) {return(NA_real_)}
            
            # Compute the density of eps
            density_eps <- ddist(distribution = dist_spec,
                                 y = eps,
                                 skew = skew,
                                 shape = shape,
                                 lambda = lambda)
            
            # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
            indicator <- ifelse(
              eps <= qdist(distribution = dist_spec,
                           p = tolerance_lvl,
                           skew = skew,
                           shape = shape,
                           lambda = lambda),
              1, 0)
            
            # Compute gradient expression
            grad_term <- (mu_grad + eps * sigma_grad) / sigma
            
            # Final result
            (CumVio_lag_j - (tolerance_lvl / 2)) * density_eps * indicator * grad_term
          }
        ))
      
      # Sum R_j_t vectors across all rows (ignoring NA entries)
      R_j_t_sum <- tibble_tmp %>%
        filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
        pull(.data[[R_j_t]]) %>%
        reduce(`+`)
      
      # Number of R_j_t without NA
      n_no_NA <- tibble_tmp %>%
        filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
        nrow()
      
      # Calculate final R_j
      R_j_list[[j]] <- (1/(tolerance_lvl * (1/3 - tolerance_lvl/4))) * (1 / n_no_NA) * R_j_t_sum
    }
    
    # Create sigma matrix
    sigma_mat <- matrix(nrow = lags, ncol = lags)
    for(i in 1:lags){
      for(j in 1:lags){
        sigma_mat[i,j] <- ifelse(i==j,1,0) + (n/est_window) * (t(as.vector(R_j_list[[i]])) %*% as.matrix(W_mean) %*% as.vector(R_j_list[[j]]))
      }
    }
    
    # Compute final correction matrix
    corr_matrix <- solve(sigma_mat)
    
    ##TEST
    #tibble_cc <<- tibble_tmp
    
  } else {
    
    # No correction if estimate = FALSE
    corr_matrix <- diag(x = 1, nrow = lags)
  }
  
  ##TEST
  #corr_matrix_save <<- corr_matrix
  
  # Test statistic calculation
  C <- n * t(rho_nj) %*% corr_matrix %*% rho_nj
  
  # p-value calculation
  p <- pchisq(q = C,
              df = lags,
              lower.tail = FALSE)
  
  # Return p-value
  return(as.numeric(p))
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
                                par_corr,
                                var_spec_est,
                                mean_spec_est,
                                dist_spec_est,
                                fixed_pars_est,
                                cores,
                                white_adjust,
                                seed=NA,
                                mincer_spec,
                                execute_additional_tsts,
                                lags_es_cc,
                                empirical){
  
  # Register cores
  registerDoParallel(cores = cores)
  
  # Estimation loop
  result_foreach <- foreach(i=1:n_loop) %dopar% {
    
    # Set seed for reproducability of results
    if(!is.na(seed)){
      current_seed <- (seed+i^2-i) %% 2147483647
      set.seed(current_seed)
    } else {
      current_seed <- NA
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
    #test <<- plot(garchsimulation_zoo)
    
    # Calculate GARCH model
    rolling_VaR_ES <- rollapply(garchsimulation_zoo,
                                width = est_window,
                                FUN = function(x) VaR_ES_forecast(data_zoo = x,
                                                                  var_spec = var_spec_est,
                                                                  mean_spec = mean_spec_est,
                                                                  dist_spec = dist_spec_est,
                                                                  tolerance_lvl = tolerance_lvl,
                                                                  estimate = estimate,
                                                                  par_corr = par_corr,
                                                                  fixed_pars_est = fixed_pars_est,
                                                                  empirical=empirical),
                                align = 'right',
                                coredata = FALSE)
    
    rolling_VaR_ES_lead <- stats::lag(rolling_VaR_ES, k = -1)
    rolling_VaR_ES_df <- as.data.frame(rolling_VaR_ES_lead)
    
    colvert_col <- setdiff(names(rolling_VaR_ES_df),c('dist_spec', 'vcov_matrix_par', 'mu_grad', 'sigma_grad', 'emp_dist_vec'))
    rolling_VaR_ES_df[colvert_col] <- lapply(rolling_VaR_ES_df[colvert_col], as.numeric)
    
    rolling_VaR_ES_df <- rolling_VaR_ES_df %>%
      mutate(Date = as.Date(rownames(rolling_VaR_ES_df)))
    VaR_ES_results_df <- inner_join(garchsimulation_df, rolling_VaR_ES_df, by = 'Date')
    
    ##TEST
    VaR_ES_results_df_save <<- VaR_ES_results_df
    
    # Filter for shortfalls
    shortfall_df <- VaR_ES_results_df %>%
      filter(Return <= VaR) %>%
      mutate(shortfall = Return - ES) %>%
      mutate(residual_t_min_1_quadr_lower_0 = ifelse(residual_t_min_1 < 0, residual_t_min_1_quadr, 0)) %>%
      mutate(indicator_residual_lower_0 = ifelse(residual_t_min_1 < 0, 1, 0))
    
    # Save number of observations that go into Mincer regression
    n_obs_mincer <- nrow(shortfall_df)
    
    # Loop over specified white adjustments
    for(i_wa in 1:length(white_adjust)){
      wa <- white_adjust[i_wa]
      
      # Mincer regression execution
      p_values_vec <- vector()
      for(j in 1:length(mincer_spec)){
        p_values_vec[j] <- mincer_regression(formula = mincer_spec[[j]][['formula']],
                                             shortfall_df = shortfall_df,
                                             h0 = mincer_spec[[j]][['h0']],
                                             white_adjust = wa)
      }
      
      # Save results
      if(i_wa == 1){
        result_lst <- list()
        result_lst_names <- names(mincer_spec)
      }
      
      for(j in 1:length(mincer_spec)){
        result_lst[[result_lst_names[j]]][[paste0('p_', wa)]] <- p_values_vec[j]
        result_lst[[result_lst_names[j]]][[paste0('p0_01_', wa)]] <- ifelse(p_values_vec[j] < 0.01, 1, 0)
        result_lst[[result_lst_names[j]]][[paste0('p0_05_', wa)]] <- ifelse(p_values_vec[j] < 0.05, 1, 0)
        result_lst[[result_lst_names[j]]][[paste0('p0_1_', wa)]] <- ifelse(p_values_vec[j] < 0.1, 1, 0)
      }
      
      # Calculate results for unconditional coverage test
      if(execute_additional_tsts){
        
        # Calculate cumulative violations
        u <- rugarch::pdist(distribution = dist_spec_est,
                            q = VaR_ES_results_df[['Return']],
                            mu = VaR_ES_results_df[['mu']],
                            sigma = VaR_ES_results_df[['sigma']],
                            skew = VaR_ES_results_df[['skew']],
                            shape = VaR_ES_results_df[['shape']],
                            lambda = VaR_ES_results_df[['lambda']])
        CumVio <- (1 / tolerance_lvl) * (tolerance_lvl - u) * ifelse(u <= tolerance_lvl, 1, 0)
        
        # Additional tests
        add_tests <- c('UC', 'CC', 'ESR')
        p_add <- vector(length = length(add_tests))
        names(p_add) <- add_tests
        
        # Execute unconditional and conditional coverage backtest for ES
        p_add['UC'] <- ES_uc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      par_corr = par_corr,
                                      est_window = est_window,
                                      Return = VaR_ES_results_df[['Return']],
                                      mu = VaR_ES_results_df[['mu']],
                                      sigma = VaR_ES_results_df[['sigma']],
                                      dist_spec = VaR_ES_results_df[['dist_spec']],
                                      skew = VaR_ES_results_df[['skew']],
                                      shape = VaR_ES_results_df[['shape']],
                                      lambda = VaR_ES_results_df[['lambda']],
                                      vcov_matrix_par_code = VaR_ES_results_df[['vcov_matrix_par']],
                                      mu_grad_code = VaR_ES_results_df[['mu_grad']],
                                      sigma_grad_code = VaR_ES_results_df[['sigma_grad']])
        
        p_add['CC'] <- ES_cc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      lags = lags_es_cc,
                                      par_corr = par_corr,
                                      est_window = est_window,
                                      Return = VaR_ES_results_df[['Return']],
                                      mu = VaR_ES_results_df[['mu']],
                                      sigma = VaR_ES_results_df[['sigma']],
                                      dist_spec = VaR_ES_results_df[['dist_spec']],
                                      skew = VaR_ES_results_df[['skew']],
                                      shape = VaR_ES_results_df[['shape']],
                                      lambda = VaR_ES_results_df[['lambda']],
                                      vcov_matrix_par_code = VaR_ES_results_df[['vcov_matrix_par']],
                                      mu_grad_code = VaR_ES_results_df[['mu_grad']],
                                      sigma_grad_code = VaR_ES_results_df[['sigma_grad']])
        
        # Execute ESR backtest
        p_add['ESR'] <- esr_backtest(r = VaR_ES_results_df[['Return']],
                                     q = VaR_ES_results_df[['VaR']],
                                     e = VaR_ES_results_df[['ES']],
                                     alpha = tolerance_lvl,
                                     version = 2)[['pvalue_twosided_asymptotic']]
        
        # Store result
        for(add_tst in add_tests){
          result_lst[[add_tst]][[paste0('p_', wa)]] <- p_add[add_tst]
          result_lst[[add_tst]][[paste0('p0_01_', wa)]] <- ifelse(p_add[add_tst] < 0.01, 1, 0)
          result_lst[[add_tst]][[paste0('p0_05_', wa)]] <- ifelse(p_add[add_tst] < 0.05, 1, 0)
          result_lst[[add_tst]][[paste0('p0_1_', wa)]] <- ifelse(p_add[add_tst] < 0.1, 1, 0)
        }
      }
    }
    
    # Save seed number and number of observations that entered Mincer regressions
    result_lst[['n_obs_mincer']] <- n_obs_mincer
    result_lst[['seed']] <- current_seed
    result_lst
  }
  
  # Organize results
  result <- list()
  
  mincer_spec_names <- if(execute_additional_tsts){
    c(names(mincer_spec), 'UC', 'CC', 'ESR')
  } else {
    names(mincer_spec)
  }
  
  for(mincer_reg_name in mincer_spec_names){
    result[[mincer_reg_name]] <- list()
    
    for(wa in white_adjust){
      keys <- c(paste0('p_', wa),
                paste0('p0_01_', wa),
                paste0('p0_05_', wa),
                paste0('p0_1_', wa))
      
      # Assign empty vectors to dynamically named list elements
      result[[mincer_reg_name]][keys] <- replicate(length(keys), vector(), simplify = FALSE)
    }
  }
  result[['n_obs_mincer']] <- vector()
  result[['seed']] <- vector()
  
  for(i in 1:n_loop){
    for(method in names(result_foreach[[i]])){
      if(method %in% c('n_obs_mincer', 'seed')){
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
  result_lst[['seed']] <- NULL
  
  row_names <- names(result_lst)
  
  col_names <- names(result_lst[[1]])
  col_names <- col_names[!startsWith(col_names, 'p_')]
  
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

#############################################################################
################### Write results into txt file     #########################
#############################################################################
write_results_to_txt <- function(name,
                                 txt_file){
  matrix <- get(paste0(name, '_matrix'))
  lst <- get(name)
  
  write(paste0('\n\n', toupper(gsub('result_','',name)), ':'), file = txt_file, append = TRUE)
  write(paste0('Number of obs. in Mincer-Regressions: ',mean(lst[['n_obs_mincer']])), file = txt_file, append = TRUE)
  capture.output(print(matrix), file = txt_file, append = TRUE)
}

########################################################################################
################### Functions convert matrix or vector into string and read it #########
########################################################################################
matrix_to_string <- function(m) {
  res_str <- paste0('matrix(nrow = ', nrow(m), ', ncol = ', ncol(m), ', byrow = FALSE, data = c(')
  elements <- as.vector(m)
  res_str <- paste0(res_str, paste(elements, collapse = ", "), '))')
  return(res_str)
}

vector_to_string <- function(v){
  res_vec <- paste0('c(', paste(v, collapse = ", "), ')')
  return(res_vec)
}

eval_string_code <- function(str_code){
  str_res <- eval(parse(text = str_code))
  return(str_res)
}

#######################################################################################
################### Function creates LaTeX table body of results              #########
#######################################################################################
matrix_to_latex_body <- function(mat,
                                 double_hline=NA) {
  col_names <- paste(colnames(mat), collapse = ' & ')
  lines <- paste0('No. & ', col_names, ' \\\\\n\\hline')
  
  for (i in 1:nrow(mat)) {
    formatted_row <- formatC(mat[i, ], format = 'f', digits = 3)
    row <- paste(formatted_row, collapse = ' & ')
    row_line <- paste0(i, ' & ', row, ' \\\\')
    
    if (i %in% double_hline) {
      row_line <- paste0(row_line, '\n\\hline\n\\hline')
    } else {
      row_line <- paste0(row_line, '\n\\hline')
    }
    
    lines <- c(lines, row_line)
  }
  
  return(paste(lines, collapse = '\n'))
}

#############################################################################
################### Write results to LaTeX format    ########################
#############################################################################
write_results_to_latex <- function(name,
                                   txt_file,
                                   double_hline = NA){
  
  if(!is.na(txt_file)){
    matrix <- get(paste0(name, '_matrix'))
    lst <- get(name)
    
    write(paste0('\n\n', toupper(gsub('result_','',name)), ':'), file = txt_file, append = TRUE)
    write(paste0('Number of obs. in Mincer-Regressions: ',mean(lst[['n_obs_mincer']])), file = txt_file, append = TRUE)
    capture.output(cat(matrix_to_latex_body(mat = matrix, double_hline = double_hline)), file = txt_file, append = TRUE)
  }
}