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
                            fixed_pars,
                            empirical,
                            seed_opt=NA){
  
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
                  solver = 'hybrid',
                  solver.control = list(n.restarts = 25,                               # gosolnp: Number of restarts
                                        n.sim = 500,                                   # gosolnp: Number of initial parameter simulation -> min objective function is taken for optimization initiation
                                        rseed = (seed_opt+111),                        # gosolnp: Seed for inital parameter simulation
                                        trace = 0))
      }, error = function(e){
        write('1', file = 'run_logs/solver_problems.txt', append = TRUE)
        return(NULL)
      }
    )
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) & garchfit@fit$convergence == 1){
      write('2', file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Try nloptr+AUGLAG+PRAXIS if previous fitting with hybrid solver failed
    if(is.null(garchfit)){
      garchfit <- tryCatch(
        {
          ugarchfit(spec = garchspec,
                    data = data_zoo,
                    solver = 'nloptr',
                    solver.control = list(solver = 8,                                 # AUGLAG+PRAXIS solver
                                          ranseed = (seed_opt+111),                   # Seed for random processes
                                          print_level = 0))
        }, error = function(e){
          write('3', file = 'run_logs/solver_problems.txt', append = TRUE)
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) & garchfit@fit$convergence == 1){
      write('4', file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Try gosolnp alone with high number of simulations and tries if previous fitting with hybrid and nloptr solver failed
    if(is.null(garchfit)){
      garchfit <- tryCatch(
        {
          ugarchfit(spec = garchspec,
                    data = data_zoo,
                    solver = 'gosolnp',
                    solver.control = list(n.restarts = 100,
                                          n.sim = 2000,
                                          rseed = (seed_opt+231),
                                          trace = 0))
        }, error = function(e){
          write('5', file = 'run_logs/solver_problems.txt', append = TRUE)
          return(NULL)
        }
      )
    }
    
    # Assign NULL to garchfit if solver did not converge
    if(!is.null(garchfit) & garchfit@fit$convergence == 1){
      write('6', file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains NAs
    if(!is.null(garchfit) & anyNA(rugarch::coef(garchfit))){
      write('7', file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
    # Assign NULL to garchfit if coef contains infinite numbers
    if(!is.null(garchfit) & any(is.infinite(rugarch::coef(garchfit)))){
      write('8', file = 'run_logs/solver_problems.txt', append = TRUE)
      garchfit <- NULL
    }
    
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
  
  if(!is.null(garchfit)){
    
    # Extracting coefficients
    coef_fit <- rugarch::coef(garchfit)
    
    # Extracting last standard deviation from fit
    st_dev_t_min_1 <- tail(rugarch::sigma(garchfit), 1)

    # Extracting last residual from fit
    residual_t_min_1 <- tail(rugarch::residuals(garchfit), 1)
    
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
    if(par_corr & estimate){
      vcov_matrix_par <- matrix_to_string(rugarch::vcov(garchfit))

      gradients <- numDeriv::jacobian(func = function(x) {
        
        x_named <- setNames(as.list(x), names(coef_fit))
        
        spec_tmp <- ugarchspec(variance.model = var_spec,
                               mean.model = mean_spec,
                               distribution.model = dist_spec,
                               fixed.pars = x_named)
        
        forecast_tmp <- ugarchforecast(fitORspec = spec_tmp,
                                       data = data_zoo,
                                       n.ahead = 1)
        
        return(c(mu = as.double(rugarch::fitted(forecast_tmp)),
                 sigma = as.double(rugarch::sigma(forecast_tmp))))
      },
      x = coef_fit,
      method = 'Richardson')
      
      mu_grad <- vector_to_string(gradients[1,])
      sigma_grad <- vector_to_string(gradients[2,])
    }
    
    # Extracting mu and sigma
    mu <- rugarch::fitted(garchforecast)
    sigma <- rugarch::sigma(garchforecast)
    
    if(empirical){
      
      # Empirical distribution
      emp_dist <- as.vector(rugarch::residuals(garchfit, standardize = TRUE))
      emp_dist_vec <- vector_to_string(emp_dist)

      # pth quantile of empirical distribution
      q_emp_dist <- stats::quantile(emp_dist,
                                    probs = tolerance_lvl,
                                    na.rm = TRUE,
                                    names = FALSE,
                                    type = 1)

      # VaR and ES calculation
      VaR <- mu + sigma * q_emp_dist
      ES <- mu + sigma * mean(emp_dist[emp_dist <= q_emp_dist])
      
    } else{
      
      # Calculate VaR
      VaR <- mu + sigma * rugarch::qdist(distribution = dist_spec,
                                         p = tolerance_lvl,
                                         skew = skew,
                                         shape = shape,
                                         lambda = lambda)
      
      # Calculate ES
      ES <- mu + sigma / tolerance_lvl * stats::integrate(f = function(x) rugarch::qdist(distribution = dist_spec,
                                                                                         p = x,
                                                                                         skew = skew,
                                                                                         shape = shape,
                                                                                         lambda = lambda),
                                                          lower = 0,
                                                          upper = tolerance_lvl)[['value']]
    }
    
    # Return results
    results <- tibble(last_date = last_date,
                      mu = as.double(mu),
                      sigma = as.double(sigma),
                      variance = as.double(sigma)^2,
                      skew = as.double(skew),
                      shape = as.double(shape),
                      lambda = as.double(lambda),
                      dist_spec = dist_spec,
                      st_dev_t_min_1 = as.double(st_dev_t_min_1),
                      variance_t_min_1 = as.double(st_dev_t_min_1)^2,
                      residual_t_min_1 = as.double(residual_t_min_1),
                      residual_t_min_1_quadr = as.double(residual_t_min_1)^2,
                      VaR = as.double(VaR),
                      ES = as.double(ES),
                      vcov_matrix_par = ifelse(par_corr, vcov_matrix_par, NA),
                      mu_grad = ifelse(par_corr, mu_grad, NA),
                      sigma_grad = ifelse(par_corr, sigma_grad, NA),
                      emp_dist_vec = ifelse(empirical, emp_dist_vec, NA))
  } else {
    results <- tibble(last_date = last_date,
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
  
  results_zoo <- zoo(x = dplyr::select(results, -'last_date'),
                     order.by = results[['last_date']])
  return(results_zoo)
}

#############################################################################
###################     Mincer Regression    ################################
#############################################################################
mincer_regression <- function(formula,
                              shortfall_tbl,
                              h0,
                              white_adjust){
  mincer_reg <- lm(formula = formula,
                   data = shortfall_tbl)
  mincer_reg_result <- car::linearHypothesis(model = mincer_reg,
                                             hypothesis.matrix = h0,
                                             test = 'F',
                                             white.adjust = white_adjust)
  return(mincer_reg_result$'Pr(>F)'[2])
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
                           sigma_grad_code,
                           empirical,
                           emp_dist_code){
  
  # Extracting non-NA cumulative violations
  not_na_CumVio <- !is.na(CumVio)
  H_hut <- CumVio[not_na_CumVio]
  
  # Mean and length of H_hut
  H_mean <- mean(H_hut)
  n <- length(H_hut)
  
  # Parameter uncertainty correction
  if(par_corr){
    
    if(empirical){
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code,
                           emp_dist_code = emp_dist_code) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code),
               emp_dist        = map(emp_dist_code, eval_string_code))
      
      # Calculate R_ES
      tibble_tmp <- tibble_tmp %>%
        mutate(R_ES_t = pmap(
          list(eps, mu_grad, sigma_grad, sigma, emp_dist),
          function(eps, mu_grad, sigma_grad, sigma, emp_dist){
            
            # Check for any NA values
            if (anyNA(c(eps, mu_grad, sigma_grad, sigma))) {return(NA_real_)}
            
            # Compute the density of eps
            kernel_dens_eps <- density(emp_dist,
                                       na.rm = TRUE)
            density_eps <- approx(x = kernel_dens_eps[['x']],
                                  y = kernel_dens_eps[['y']],
                                  xout = eps)[['y']]
            
            # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
            indicator <- ifelse(
              eps <= stats::quantile(emp_dist,
                                     probs = tolerance_lvl,
                                     na.rm = TRUE,
                                     names = FALSE,
                                     type = 8), # Approximately median unbiased regardless of distribution
              1, 0)
            
            # Compute gradient expression
            grad_term <- (mu_grad + eps * sigma_grad) / sigma
            
            # Final result
            return(density_eps * indicator * grad_term)
          }
        ))
      
    } else {
      
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
          function(eps, dist_spec, skew, shape, lambda, mu_grad, sigma_grad, sigma){
            
            # Check for any NA values
            if (anyNA(c(eps, dist_spec, mu_grad, sigma_grad, sigma))) {return(NA_real_)}
            
            # Compute the density of eps
            density_eps <- rugarch::ddist(distribution = dist_spec,
                                          y = eps,
                                          skew = skew,
                                          shape = shape,
                                          lambda = lambda)
            
            # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
            indicator <- ifelse(
              eps <= rugarch::qdist(distribution = dist_spec,
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
    }
    
    # Sum R_ES_t vectors across all rows (ignoring NA entries)
    R_ES_sum <- tibble_tmp %>%
      dplyr::filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      pull(R_ES_t) %>%
      purrr::reduce(`+`)
    
    # Number of R_ES_t without NA
    n_no_NA <- tibble_tmp %>%
      dplyr::filter(!map_lgl(R_ES_t, ~ anyNA(.))) %>%
      nrow()
    
    # Calculate final R_ES
    R_ES <- (1 / (tolerance_lvl * n_no_NA)) * R_ES_sum
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      dplyr::filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- purrr::reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
    # Final calculation of correction factor
    corr_fac <- (n/est_window) * t(as.vector(R_ES)) %*% as.matrix(W_mean) %*% as.vector(R_ES)
    
    ##TEST
    #tibble_uc <<- tibble_tmp
    #W_mean_save_uc <<- W_mean
    #R_ES_save <<- R_ES
    
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
                           sigma_grad_code,
                           empirical,
                           emp_dist_code){
  
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
    
    if(empirical){
      
      # Create tibble with required data
      tibble_tmp <- tibble(Return = Return,
                           mu = mu,
                           sigma = sigma,
                           vcov_matrix_par_code = vcov_matrix_par_code,
                           mu_grad_code = mu_grad_code,
                           sigma_grad_code = sigma_grad_code,
                           emp_dist_code = emp_dist_code,
                           CumVio = CumVio) %>%
        mutate(eps = (Return - mu) / sigma,
               vcov_matrix_par = map(vcov_matrix_par_code, eval_string_code),
               mu_grad         = map(mu_grad_code, eval_string_code),
               sigma_grad      = map(sigma_grad_code, eval_string_code),
               emp_dist        = map(emp_dist_code, eval_string_code))
      
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
            list(eps, mu_grad, sigma_grad, sigma, emp_dist, .data[[CumVio_lag_j]]),
            function(eps, mu_grad, sigma_grad, sigma, emp_dist, CumVio_lag_j) {
              
              # Check for any NA values
              if (anyNA(c(eps, mu_grad, sigma_grad, sigma, CumVio_lag_j))) {return(NA_real_)}
              
              # Compute the empirical density of eps
              kernel_dens_eps <- density(emp_dist,
                                         na.rm = TRUE)
              density_eps <- approx(x = kernel_dens_eps[['x']],
                                    y = kernel_dens_eps[['y']],
                                    xout = eps)[['y']]
              
              # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
              indicator <- ifelse(
                eps <= stats::quantile(emp_dist,
                                       probs = tolerance_lvl,
                                       na.rm = TRUE,
                                       names = FALSE,
                                       type = 8), # Approximately median unbiased regardless of distribution
                1, 0)
              
              # Compute gradient expression
              grad_term <- (mu_grad + eps * sigma_grad) / sigma
              
              # Final result
              return((CumVio_lag_j - (tolerance_lvl / 2)) * density_eps * indicator * grad_term)
            }
          ))
        
        # Sum R_j_t vectors across all rows (ignoring NA entries)
        R_j_t_sum <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          pull(.data[[R_j_t]]) %>%
          purrr::reduce(`+`)
        
        # Number of R_j_t without NA
        n_no_NA <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          nrow()
        
        # Calculate final R_j
        R_j_list[[j]] <- (1/(tolerance_lvl * (1/3 - tolerance_lvl/4))) * (1 / n_no_NA) * R_j_t_sum
      }
      
    } else {
      
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
              density_eps <- rugarch::ddist(distribution = dist_spec,
                                            y = eps,
                                            skew = skew,
                                            shape = shape,
                                            lambda = lambda)
              
              # Compute indicator: 1 if eps <= quantile of tolerance_lvl, else 0
              indicator <- ifelse(
                eps <= rugarch::qdist(distribution = dist_spec,
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
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          pull(.data[[R_j_t]]) %>%
          purrr::reduce(`+`)
        
        # Number of R_j_t without NA
        n_no_NA <- tibble_tmp %>%
          dplyr::filter(!map_lgl(.data[[R_j_t]], ~ anyNA(.))) %>%
          nrow()
        
        # Calculate final R_j
        R_j_list[[j]] <- (1/(tolerance_lvl * (1/3 - tolerance_lvl/4))) * (1 / n_no_NA) * R_j_t_sum
      }
    }
    
    # Calculate mean of vcov_matrix_par
    valid_vcov_mat <- tibble_tmp %>%
      dplyr::filter(!map_lgl(vcov_matrix_par, ~ anyNA(.) || is.null(.))) %>%
      pull(vcov_matrix_par)
    W_mean <- purrr::reduce(valid_vcov_mat, `+`) / length(valid_vcov_mat)
    
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
    #W_mean_save_cc <<- W_mean
    #R_j_list_save <<- R_j_list
    
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
                                n_boot_de,
                                empirical){
  
  # Clear started_loops.txt and finished_loops.txt
  write(x = '\n\n',
        file = 'run_logs/started_loops.txt',
        append = FALSE)
  write(x = '\n\n',
        file = 'run_logs/finished_loops.txt',
        append = FALSE)
  
  # Register cores
  registerDoParallel(cores = cores)
  
  # Estimation loop
  result_foreach <- foreach(i=1:n_loop) %dopar% {
    
    # Register started loops
    write(x = paste0(i),
          file = 'run_logs/started_loops.txt',
          append = TRUE)
    
    # Set seed for reproducibility of results
    if(!is.na(seed)){
      current_seed <- (seed+i^2-i) %% 2147483647
      set.seed(current_seed)
    } else {
      current_seed <- NA
    }
    
    # Progress messages
    if(cores==1){
      cat(paste0('Iteration Nr.: ', i, '\n'))
    } else if(i %% cores == 0){
      write(x = paste0('Loop ', i, ' started at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
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
    
    # Create tibble and zoo object of simulation
    sim_matrix <- garchsimulation@path$seriesSim
    colnames(sim_matrix) <- 'Return'
    garchsimulation_tbl <- as_tibble(sim_matrix) %>%
      mutate(Date = base::as.Date('2000-01-01') + 1:nrow(sim_matrix))
    garchsimulation_zoo <- zoo(x = garchsimulation_tbl[['Return']],
                               order.by = garchsimulation_tbl[['Date']])
    #plot_s <<- plot(garchsimulation_zoo)
    
    ##TEST
    #garchsimulation_zoo_s <<- garchsimulation_zoo
    
    # Calculate GARCH model
    rolling_VaR_ES <- rollapply(data = garchsimulation_zoo,
                                width = est_window,
                                FUN = function(x) VaR_ES_forecast(data_zoo = x,
                                                                  var_spec = var_spec_est,
                                                                  mean_spec = mean_spec_est,
                                                                  dist_spec = dist_spec_est,
                                                                  tolerance_lvl = tolerance_lvl,
                                                                  estimate = estimate,
                                                                  par_corr = par_corr,
                                                                  fixed_pars = fixed_pars_est,
                                                                  empirical = empirical,
                                                                  seed_opt = current_seed),
                                align = 'right',
                                coredata = FALSE)
    
    rolling_VaR_ES_lead <- stats::lag(rolling_VaR_ES,
                                      k = -1)
    rolling_VaR_ES_tbl <- as_tibble(rolling_VaR_ES_lead) %>%
      mutate(across(-c(dist_spec, vcov_matrix_par, mu_grad, sigma_grad, emp_dist_vec), as.numeric),
             Date = index(rolling_VaR_ES_lead))
    VaR_ES_results_tbl <- inner_join(garchsimulation_tbl, rolling_VaR_ES_tbl, by = 'Date')
    
    ##TEST
    #VaR_ES_results_tbl_save <<- VaR_ES_results_tbl
    
    # Report number of NAs in NA_infomation.txt
    n_NAs <- sum(is.na(VaR_ES_results_tbl[['ES']]))
    
    if(n_NAs > 0){
      write(paste0('Number of NAs in loop ', i, ': ', n_NAs, '\n'),
            file = 'run_logs/NA_information.txt',
            append = TRUE)
    }
    
    # Filter for shortfalls
    shortfall_tbl <- VaR_ES_results_tbl %>%
      dplyr::filter(Return <= VaR) %>%
      mutate(shortfall = Return - ES,
             residual_t_min_1_quadr_lower_0 = ifelse(residual_t_min_1 < 0, residual_t_min_1_quadr, 0),
             indicator_residual_lower_0 = ifelse(residual_t_min_1 < 0, 1, 0))
    
    # Save number of observations that go into Mincer regression
    n_obs_mincer <- nrow(shortfall_tbl)
    
    # Loop over specified white adjustments
    for(i_wa in 1:length(white_adjust)){
      wa <- white_adjust[i_wa]
      
      # Mincer regression execution
      p_values_vec <- vector()
      for(j in 1:length(mincer_spec)){
        p_values_vec[j] <- mincer_regression(formula = mincer_spec[[j]][['formula']],
                                             shortfall_tbl = shortfall_tbl,
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
      
      # Compute additional tests if execute_additional_tsts = TRUE
      if(execute_additional_tsts){
        
        # Compute cdf value at return
        if(empirical){
          
          # With empirical distribution
          u <- VaR_ES_results_tbl %>%
            mutate(standardized_Return = (Return - mu) / sigma,
                   emp_dist = map(emp_dist_vec, eval_string_code)) %>%
            mutate(u = map2_dbl(.x = emp_dist,
                                .y = standardized_Return,
                                .f = ~ mean(.x <= .y))) %>%
            pull(u)
          
        } else {
          
          # With parametric distribution assumption
          u <- rugarch::pdist(distribution = dist_spec_est,
                              q = VaR_ES_results_tbl[['Return']],
                              mu = VaR_ES_results_tbl[['mu']],
                              sigma = VaR_ES_results_tbl[['sigma']],
                              skew = VaR_ES_results_tbl[['skew']],
                              shape = VaR_ES_results_tbl[['shape']],
                              lambda = VaR_ES_results_tbl[['lambda']])
        }
        
        # Compute cumulative violations
        CumVio <- (1 / tolerance_lvl) * (tolerance_lvl - u) * ifelse(u <= tolerance_lvl, 1, 0)
        
        # Create vector with names of additional tests
        add_tests <- names_of_add_tests(n_boot_de = n_boot_de,
                                        estimate = estimate,
                                        par_corr = par_corr)
        p_add <- vector(length = length(add_tests))
        names(p_add) <- add_tests
        
        # Execute unconditional and conditional coverage backtest for ES
        p_add['UC'] <- ES_uc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      par_corr = FALSE)
        
        p_add['CC'] <- ES_cc_backtest(CumVio = CumVio,
                                      tolerance_lvl = tolerance_lvl,
                                      lags = lags_es_cc,
                                      par_corr = FALSE)
        
        # Execute unconditional and conditional coverage backtest with bootstrap for ES
        if(!is.na(n_boot_de)){
          if(!is.na(seed)){set.seed(current_seed)}
          shortfall_de_test_boot_res <- tstests::shortfall_de_test(x = u[!is.na(u)],
                                                                   alpha = tolerance_lvl,
                                                                   lags = lags_es_cc,
                                                                   boot = TRUE,
                                                                   n_boot = n_boot_de)
          
          p_add['UC_boot'] <- shortfall_de_test_boot_res[['p_value']][1]
          p_add['CC_boot'] <- shortfall_de_test_boot_res[['p_value']][2]
        }
        
        # Execute unconditional and conditional coverage backtest with parameter uncertainty correction for ES
        if(estimate & par_corr){
          p_add['UC_par_corr'] <- ES_uc_backtest(CumVio = CumVio,
                                                 tolerance_lvl = tolerance_lvl,
                                                 par_corr = TRUE,
                                                 est_window = est_window,
                                                 Return = VaR_ES_results_tbl[['Return']],
                                                 mu = VaR_ES_results_tbl[['mu']],
                                                 sigma = VaR_ES_results_tbl[['sigma']],
                                                 dist_spec = VaR_ES_results_tbl[['dist_spec']],
                                                 skew = VaR_ES_results_tbl[['skew']],
                                                 shape = VaR_ES_results_tbl[['shape']],
                                                 lambda = VaR_ES_results_tbl[['lambda']],
                                                 vcov_matrix_par_code = VaR_ES_results_tbl[['vcov_matrix_par']],
                                                 mu_grad_code = VaR_ES_results_tbl[['mu_grad']],
                                                 sigma_grad_code = VaR_ES_results_tbl[['sigma_grad']],
                                                 empirical = empirical,
                                                 emp_dist_code = VaR_ES_results_tbl[['emp_dist_vec']])
          
          p_add['CC_par_corr'] <- ES_cc_backtest(CumVio = CumVio,
                                                 tolerance_lvl = tolerance_lvl,
                                                 lags = lags_es_cc,
                                                 par_corr = TRUE,
                                                 est_window = est_window,
                                                 Return = VaR_ES_results_tbl[['Return']],
                                                 mu = VaR_ES_results_tbl[['mu']],
                                                 sigma = VaR_ES_results_tbl[['sigma']],
                                                 dist_spec = VaR_ES_results_tbl[['dist_spec']],
                                                 skew = VaR_ES_results_tbl[['skew']],
                                                 shape = VaR_ES_results_tbl[['shape']],
                                                 lambda = VaR_ES_results_tbl[['lambda']],
                                                 vcov_matrix_par_code = VaR_ES_results_tbl[['vcov_matrix_par']],
                                                 mu_grad_code = VaR_ES_results_tbl[['mu_grad']],
                                                 sigma_grad_code = VaR_ES_results_tbl[['sigma_grad']],
                                                 empirical = empirical,
                                                 emp_dist_code = VaR_ES_results_tbl[['emp_dist_vec']])
        }
        
        # Create VaR_ES_results_tbl without NAs in Return, VaR and ES
        if(!exists('VaR_ES_results_tbl_noNA')){
          VaR_ES_results_tbl_noNA <- VaR_ES_results_tbl %>%
            tidyr::drop_na(Return, VaR, ES)
        }
        
        # Execute Auxiliary ESR backtest
        if(!is.na(seed)){set.seed(current_seed)}
        p_add['Auxiliary_ESR'] <- esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                       q = VaR_ES_results_tbl_noNA[['VaR']],
                                                       e = VaR_ES_results_tbl_noNA[['ES']],
                                                       alpha = tolerance_lvl,
                                                       version = 2)[['pvalue_twosided_asymptotic']]
        
        # Execute Strict ESR backtest
        if(!is.na(seed)){set.seed(current_seed)}
        p_add['Strict_ESR'] <- esback::esr_backtest(r = VaR_ES_results_tbl_noNA[['Return']],
                                                    e = VaR_ES_results_tbl_noNA[['ES']],
                                                    alpha = tolerance_lvl,
                                                    version = 1)[['pvalue_twosided_asymptotic']]
        
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
    
    # Progress messages
    if(cores==1){
      cat(paste0('Finished iteration no.: ', i, '\n'))
    } else if(i %% cores == 0){
      write(x = paste0('Loop ', i, ' finished at ', Sys.time(), '\n'),
            file = 'run_logs/progress_information.txt',
            append = TRUE)
    }
    
    # Register finished loops
    write(x = paste0(i),
          file = 'run_logs/finished_loops.txt',
          append = TRUE)
    
    # Result list
    result_lst
  }
  
  ##TEST
  # saveRDS(object = result_foreach,
  #         file = 'result_foreach.RData')
  
  # Organize results
  result <- list()
  
  # Create vector with names of executed tests
  if(execute_additional_tsts){
    add_tests <- names_of_add_tests(n_boot_de = n_boot_de,
                                    estimate = estimate,
                                    par_corr = par_corr)
    test_names <- c(names(mincer_spec), add_tests)
  } else {
    test_names <- names(mincer_spec)
  }
  
  for(test_name in test_names){
    result[[test_name]] <- list()
    
    for(wa in white_adjust){
      keys <- c(paste0('p_', wa),
                paste0('p0_01_', wa),
                paste0('p0_05_', wa),
                paste0('p0_1_', wa))
      
      # Assign empty vectors to dynamically named list elements
      result[[test_name]][keys] <- replicate(length(keys), vector(), simplify = FALSE)
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

####################################################################################
################### For parallel loop: Create vector with names of test ############
####################################################################################
names_of_add_tests <- function(n_boot_de,
                               estimate,
                               par_corr){
  
  add_tst <- c('UC', 'CC', 'Auxiliary_ESR', 'Strict_ESR')
  if(!is.na(n_boot_de)){
    add_tst <- c(add_tst, 'UC_boot', 'CC_boot')
  }
  if(estimate & par_corr){
    add_tst <- c(add_tst, 'UC_par_corr', 'CC_par_corr')
  }
  
  add_tst <- c(sort(add_tst[grepl('^UC', add_tst)]), sort(add_tst[grepl('^CC', add_tst)]), 'Auxiliary_ESR', 'Strict_ESR')
  return(add_tst)
}

#############################################################################
################### Create result matrix     ################################
#############################################################################
create_result_matrix <- function(result_lst,
                                 na_rm = FALSE){
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
      result_matrix[row, col] <- round(mean(result_lst[[row]][[col]], na.rm = na_rm), digits = 3)
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
  
  write(paste0('\n\n', toupper(gsub('result_', '', name)), ':'), file = txt_file, append = TRUE)
  write(paste0('Number of obs. in Mincer-Regressions: ', mean(lst[['n_obs_mincer']], na.rm = TRUE)), file = txt_file, append = TRUE)
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
                                 double_hline = NA) {
  col_names <- paste(colnames(mat), collapse = ' & ')
  lines <- paste0('No. & ', col_names, ' \\\\\n\\hline')
  
  for (i in 1:nrow(mat)){
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
    
    write(paste0('\n\n', toupper(gsub('result_', '', name)), ':'), file = txt_file, append = TRUE)
    write(paste0('Number of obs. in Mincer-Regressions: ', mean(lst[['n_obs_mincer']])), file = txt_file, append = TRUE)
    capture.output(cat(matrix_to_latex_body(mat = matrix, double_hline = double_hline)), file = txt_file, append = TRUE)
  }
}

###################################################################################
######## Create barplot from result matrix for specified white adjustment #########
###################################################################################
create_result_barplot <- function(result_mat,
                                  white_adjust_name){
  
  # Create vector with rejection proportions
  rejection_prop <- as_tibble(result_mat) %>%
    select(contains(white_adjust_name)) %>%
    mutate(p0_01 = .data[[paste0('p0_01_', white_adjust_name)]]) %>%
    mutate(p0_05 = .data[[paste0('p0_05_', white_adjust_name)]] - p0_01) %>%
    mutate(p0_1 = .data[[paste0('p0_1_', white_adjust_name)]] - (p0_01 + p0_05)) %>%
    select(p0_01, p0_05, p0_1) %>%
    as.matrix() %>%
    t() %>%
    as.vector()
  
  # Create tibble with data for plot
  data_tbl <- tibble(category = rep(paste('No.', 1:nrow(result_mat)), each = 3),
                     part = rep(c('Size 0.01', 'Size 0.05', 'Size 0.1'), times = nrow(result_mat)),
                     value = rejection_prop) %>%
    mutate(category = factor(category, levels = rev(paste('No.', 1:nrow(result_mat)))),
           part = factor(part, levels = c('Size 0.1', 'Size 0.05', 'Size 0.01')))
  
  
  # Plot
  ggplot(data_tbl, aes(x = category, y = value, fill = part)) +
    geom_col(width = 0.4, color = 'black') +
    geom_hline(yintercept = 0.01, color = '#c6dbef', linewidth = 1) +
    geom_hline(yintercept = 0.05, color = '#6baed6', linewidth = 1) +
    geom_hline(yintercept = 0.1,  color = '#2171b5', linewidth = 1) +
    scale_fill_manual(values = c('Size 0.01' = '#c6dbef',
                                 'Size 0.05' = '#6baed6',
                                 'Size 0.1' = '#2171b5'),
                      guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(breaks = c(0.01, 0.05, 0.1, seq(0.15, 1, by = 0.05)),
                       labels = scales::number_format(accuracy = 0.01)) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(x = '',
         y = 'Rejection Proportion',
         title = 'Empirical Rejection Proportions at Nominal Levels',
         fill = 'Nominal Levels') +
    theme(axis.text.y = element_text(size = 9),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = 'grey70'),
          panel.grid.minor.x = element_line(color = 'grey85'))
}