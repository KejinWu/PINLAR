require(doSNOW)
require(parallel)
require(foreach)
require(doParallel)
require(robustbase)

##########################################################################################  
# This script shows a simple example of predicting a type of NLAR model in the format Nonlinear_fun = "a + log(b + abs(X_0))".
# L_1 and L_2 optimal point predictions will be generated.
# The QPI with fitted or predictive residuals will be built.
# The PPI centered at L_1 or L_2 point predictions with fitted or predictive residuals will be built.
# The code also provides a way to run replications of simulation studies parallelly, which saves the running time when a large number of replication are required.

# WARNING: this script only covers one specific type of NLAR format (i.e., "a + log(b + abs(X_0))"). 
#          To perform the simulation studies of other NLAR models in the paper, please change the parameter Nonlinear_fun. 
#          Also, please replcace the command "nls(y~a + log(b^2 + abs(x))"  by other model accordingly. 
#          Lastly, we should also use other Simulate_data functions for other types of NLAR models; all functions to simulate different NLAR model data can be found in a separate R file namely "All_simulated_data"
    

##########################################################################################  



##########################################################################################
# Function to simulate NLAR time series data according to the format "a + log(b + abs(X_0))".
##########################################################################################


#' @param burn_in the length of data would be discarded to generate stationary simulated data.
#' @param data_len the length of data at hand.
#' @param Nonlinear_fun the non-linear model form.
#' @param Dis_initial_x the distribution to generate the initial value to ignite the generating process
#' @param Dis_error the error distribution.
#' @note The output will be simulated data for a specifc NLAR format indicated by Nonlinear_fun. 

Simulate_data = function(burn_in = 5000,data_len = 2000 , Nonlinear_fun = "sin(X_0)" ,Dis_initial_x = "runif(1,min = -1, max = 1)",Dis_error = "runif(1,min = 0, max = 2*pi)"){
  Simulated_X = c()
  Total_len = burn_in + data_len
  X_0 = eval(parse(text = Dis_initial_x))
  for (i in c(1:Total_len)) {
    X_0 = eval(parse(text = Nonlinear_fun)) + eval(parse(text = Dis_error))
    Simulated_X = c(Simulated_X,X_0)
  }
  return(Simulated_X)
}


##########################################################################################
# Bootstrap/simulation prediction
##########################################################################################


#' @param Pre_step the prediction horizon.
#' @param quantiles indicates the nominal confidence-level, i.e., quantiles[2] - quantiles[1].
#' @param M the number of pseudo values we generate for future variables (e.g., X_{T+1}). Then, the L_2 and L_1 optimal point predictions can be computed.
#' @param B the number of the forward bootstrap series generated to determine the Pertinent Prediction Interval (i.e., the K parameter in the paper)
#' @param data_len the length of data at hand.
#' @param burn_in the length of data would be discarded to generate stationary simulated data.
#' @param Nonlinear_fun the non-linear model form.
#' @param Dis_error the error distribution.
#' @param residuals_type indicates which residual type will be used to build the prediction interval; choices are "fitted" and "predictive".
#' @note The output will be the true future simulated value, the bootstrap-based L_1 and L_2 point predictions; the QPI with fitted or predictive residuals, and the PPI centered at L_2 and L_1 optimal point with fitted or predictive residuals predictions. 

# Bootstrap/simulation prediction
Boot_pre_interval = function(Pre_step = 5, quantiles = c(0.025,0.975), M  = 50000, B = 500, data_len = 2000, burn_in = 5000, Nonlinear_fun = "log(X_0^2)", Dis_error = "rnorm(1)", residuals_type = "fitted"){
  Pre_Pre_interval = matrix(ncol = Pre_step, nrow = 9) 
  rownames(Pre_Pre_interval) = c("True simulated future","L_2 Bootstrap point pre", "L_1 Bootstrap point pre", "Low quantile PI", "High quantile PI", "Low L_2 Pertinent PI",  "High L_2 Pertinent PI", "Low L_1 Pertinent PI", "High L_1 Pertinent PI" )
  
  # Simulate data
  Boot_mat = matrix(ncol = Pre_step, nrow = M)
  if(residuals_type == "fitted"){ # Get fitted residuals
    
    sim_data = Simulate_data(data_len = data_len,burn_in = burn_in, Nonlinear_fun = Nonlinear_fun,Dis_error = Dis_error)
    x_T = sim_data[length(sim_data)]
    
    True_real = c()
    X_0 = x_T
    # Find true future value by simulation
    for (h in c(1:Pre_step)) {
      X_0 = eval(parse(text = Nonlinear_fun)) + eval(parse(text = Dis_error))
      True_real = c(True_real, X_0)
      #X_0 = log(X_0^2) + rnorm(1)
      #print(X_0)
    }
    Pre_Pre_interval[1,] = True_real
    
    # Find quantile PI and corresponding center value by bootstrap
    # Do estimation first
    
    data = sim_data[(burn_in+1):(burn_in + data_len)]
    n = length(data)
    y = data[2:n]
    x = data[1:(n-1)]
    try_nls = try(nls(y~a + log(b^2 + abs(x)), start = list(a=1,b=1)),silent = TRUE)
    summary_res = summary(try_nls)
    
    while(try(summary(try_nls)[2], silent = TRUE)
          == "try-error"){
      cat("One simulation fails in estimation procedure, do it again \n")
      sim_data = Simulate_data(data_len = data_len,burn_in = burn_in, Nonlinear_fun = Nonlinear_fun,Dis_error = Dis_error)
      x_T = sim_data[length(sim_data)]
      
      True_real = c()
      X_0 = x_T
      # Find true future value by simulation
      for (h in c(1:Pre_step)) {
        X_0 = eval(parse(text = Nonlinear_fun)) + eval(parse(text = Dis_error))
        True_real = c(True_real, X_0)
        #X_0 = log(X_0^2) + rnorm(1)
        #print(X_0)
      }
      Pre_Pre_interval[1,] = True_real
      # Find quantile PI and corresponding center value by bootstrap
      # Do estimation first
      #Boot_mat = matrix(ncol = Pre_step, nrow = M)
      data = sim_data[(burn_in+1):(burn_in + data_len)]
      n = length(data)
      y = data[2:n]
      x = data[1:(n-1)]
      try_nls = try(nls(y~a + log(b^2 + abs(x)), start = list(a=1,b=1)),silent = TRUE)
      summary_res = summary(try_nls)
    }
    
    #nls_res = nls(y~log(a*(x^2)), start = list(a = 1))
    if(try_nls$convInfo$stopMessage != "converged"){
      cat("Caution: the non-linear estimation process is not converged")
    }
    summary_res = summary(try_nls)
    estimate_a = summary_res$parameters[1]
    estimate_b = summary_res$parameters[2]
    #estimate_c = summary_res$parameters[3]
    residuals = summary_res$residuals
    residuals = residuals - mean(residuals)
    Est_Nonlinear_fun = paste( estimate_a,"+log((",estimate_b,")^2+abs(X_0))",sep = "")
    
  }else{ # Get the predictive residuals, do all nls on delete data set
    rerun = TRUE 
    predictive_residual = c()
    while(rerun){
      sim_data = Simulate_data(data_len = data_len,burn_in = burn_in, Nonlinear_fun = Nonlinear_fun,Dis_error = Dis_error)
      x_T = sim_data[length(sim_data)]
      
      True_real = c()
      X_0 = x_T
      # Find true future value by simulation
      for (h in c(1:Pre_step)) {
        X_0 = eval(parse(text = Nonlinear_fun)) + eval(parse(text = Dis_error))
        True_real = c(True_real, X_0)
        #X_0 = log(X_0^2) + rnorm(1)
        #print(X_0)
      }
      Pre_Pre_interval[1,] = True_real
      
      # Find quantile PI and corresponding center value by bootstrap
      # Do estimation first
      
      data = sim_data[(burn_in+1):(burn_in + data_len)]
      n = length(data)
      y = data[2:n]
      x = data[1:(n-1)]
      for(i in c(1:n)){ # For loop to do all nls estimations
        if(i == n){ # When i == n, we do nls estimation on the whole dataset
          try_nls = try(nls(y~a + log(b^2 + abs(x)), start = list(a=1,b=1)),silent = TRUE)
          if(try(summary(try_nls)[2], silent = TRUE)
             == "try-error"){
            rerun = TRUE
            break
          }  
        }else{
          y_d = y[-i]
          x_d = x[-i]
          try_nls_d = try(nls(y_d~a + log(b^2 + abs(x_d)), start = list(a=1,b=1)),silent = TRUE)
          if(try(summary(try_nls_d)[2], silent = TRUE)
             == "try-error"){
            rerun = TRUE
            break
          }
          summary_res_d = summary(try_nls_d)
          estimate_a_d = summary_res_d$parameters[1]
          estimate_b_d = summary_res_d$parameters[2]
          predictive_residual = c(predictive_residual, (y[i] - estimate_a_d - log(estimate_b_d^2 + abs(x[i]))))
          #estimate_c = summary_res$parameters[3]   
        }
        
        if(i == n){
          rerun = FALSE
        }
        
      }     
    }
    
    summary_res = summary(try_nls)
    estimate_a = summary_res$parameters[1]
    estimate_b = summary_res$parameters[2]
    #estimate_c = summary_res$parameters[3]
    
    residuals = predictive_residual - mean(predictive_residual) # Use the predictive residuals
    # Get the estimated non-linear model
    Est_Nonlinear_fun = paste( estimate_a,"+log((",estimate_b,")^2+abs(X_0))",sep = "")
    
  }   
  
  # Do prediction
  for (m in c(1:M)) {
    X_0 = x_T
    for (h in c(1:Pre_step)) {
      X_0 = eval(parse(text = Est_Nonlinear_fun)) + sample(residuals,1,replace = T)
      Boot_mat[m,h] = X_0
      #X_0 = log(X_0^2) + rnorm(1)
      #print(X_0)
    }
  }
  
  Pre_Pre_interval[2,] = colMeans(Boot_mat)
  Pre_Pre_interval[3,] = colMedians(Boot_mat)
  
  
  PI_quantile = apply(Boot_mat, 2, function(x){quantile(x,probs = c(quantiles[1],quantiles[2]))})
  Pre_Pre_interval[4,] = PI_quantile[1,]
  Pre_Pre_interval[5,] = PI_quantile[2,]
  
  
  # Do pertinent prediction interval; We do bootstrap B times to compute predictive root
  Boot_predictive_root_L2 = matrix(nrow = B, ncol = Pre_step)
  Boot_predictive_root_L1 = matrix(nrow = B, ncol = Pre_step)
  for( b in c(1:B)){
    #print(b)
    X_0 = sample(data,1,replace = T)
    Boot_simulated_x = c()
    for (i in c(1:data_len)) {
      X_0 = eval(parse(text = Est_Nonlinear_fun)) + sample(residuals,1,replace = T)
      Boot_simulated_x = c(Boot_simulated_x,X_0)
    }
    #Boot_simulated_x[data_len] = x_T
    X_0 = x_T
    # Find the '' TRUE'' value in the bootstrap world
    True_bootstrap = c()
    for (h in c(1:Pre_step)) {
      X_0 = eval(parse(text = Est_Nonlinear_fun)) + sample(residuals,1,replace = T)
      True_bootstrap = c(True_bootstrap,X_0)
      #X_0 = log(X_0^2) + rnorm(1)
      #print(X_0)
    }
    # Find the bootstrap estimation in the bootstrap world
    Boot_Boot_mat = matrix(ncol = Pre_step, nrow = M)
    n = length(Boot_simulated_x)
    y = Boot_simulated_x[2:n]
    x = Boot_simulated_x[1:(n-1)]
    try_boot_nls = try(nls(y~a + log(b^2 + abs(x)), start = list(a=1,b=1)),silent = TRUE)
    #boot_summary_res = summary(try_boot_nls)
    #print(boot_summary_res)
    while(try(summary(try_boot_nls)[2], silent = TRUE)
          == "try-error"){
      cat("One simulation fails in estimation procedure, do it again \n")
      X_0 = sample(data,1,replace = T)
      Boot_simulated_x = c()
      for (i in c(1:data_len)) {
        X_0 = eval(parse(text = Est_Nonlinear_fun)) + sample(residuals,1,replace = T)
        Boot_simulated_x = c(Boot_simulated_x,X_0)
      }
      n = length(Boot_simulated_x)
      y = Boot_simulated_x[2:n]
      x = Boot_simulated_x[1:(n-1)]
      try_boot_nls = try(nls(y~a + log(b^2 + abs(x)), start = list(a=1,b=1)),silent = TRUE)
      #boot_summary_res = summary(try_boot_nls)
      
    }
    
    if(try_boot_nls$convInfo$stopMessage != "converged"){
      cat("Caution: the non-linear estimation process is not converged")
    }
    boot_summary_res = summary(try_boot_nls)
    Boot_estimate_a = boot_summary_res$parameters[1]
    Boot_estimate_b = boot_summary_res$parameters[2]
    #Boot_estimate_c = boot_summary_res$parameters[3]
    Boot_Est_Nonlinear_fun = paste( Boot_estimate_a,"+log((",Boot_estimate_b,")^2+abs(X_0))",sep = "")
    # Do Bootstrap prediction in the bootstrap world
    for (m in c(1:M)) {
      X_0 = x_T
      for (h in c(1:Pre_step)) {
        X_0 = eval(parse(text = Boot_Est_Nonlinear_fun)) + sample(residuals,1,replace = T)
        #print(X_0)
        Boot_Boot_mat[m,h] = X_0
        #X_0 = log(X_0^2) + rnorm(1)
        #print(X_0)
      }
    }
    Boot_boot_pre_L2 = colMeans(Boot_Boot_mat)
    Boot_boot_pre_L1 = colMedians(Boot_Boot_mat)
    Boot_predictive_root_L2[b,] = True_bootstrap - Boot_boot_pre_L2
    Boot_predictive_root_L1[b,] = True_bootstrap - Boot_boot_pre_L1
  }
  # Get the pertinent prediction interval
  PR_quantile_L2 = apply(Boot_predictive_root_L2, 2, function(x){quantile(x,probs = c(quantiles[1],quantiles[2]))})
  PR_quantile_L1 = apply(Boot_predictive_root_L1, 2, function(x){quantile(x,probs = c(quantiles[1],quantiles[2]))})
  
  Pre_Pre_interval[6,] = PR_quantile_L2[1,] + Pre_Pre_interval[2,]
  Pre_Pre_interval[7,] = PR_quantile_L2[2,] + Pre_Pre_interval[2,]
  Pre_Pre_interval[8,] = PR_quantile_L1[1,] + Pre_Pre_interval[3,]
  Pre_Pre_interval[9,] = PR_quantile_L1[2,] + Pre_Pre_interval[3,]
  
  
  
  return(Pre_Pre_interval)
}


##########################################################################################
# A parallel way to run replications of simulations
##########################################################################################

#' @param Rep the total number of replications we wanted to perform for a simulation study.
#' @param Pre_step the prediction horizon.
#' @param residuals_type indicates which residual type will be used to build the prediction interval; choices are "fitted" and "predictive".
#' @param quantiles indicates the nominal confidence-level, i.e., quantiles[2] - quantiles[1].
#' @param M the number of pseudo values we generate for future variables (e.g., X_{T+1}). Then, the L_2 and L_1 optimal point predictions can be computed.
#' @param B the number of the forward bootstrap series generated to determine the Pertinent Prediction Interval.
#' @param data_len the length of data at hand.
#' @param burn_in the length of data would be discarded to generate stationary simulated data.
#' @param Nonlinear_fun the non-linear model form.
#' @param Dis_error the error distribution.
#' @param n_core the number of CPUs used to run all replications of simulations.
#' @note The output will include all simulation results for all replications.
#' @note To determine the appropriate value for the parameter n_core, we can use the command detectCores() to find how many CPUs our PCs have.


Do_replication_parallel = function(Rep = 200,Pre_step = 5,  residuals_type = "fitted", quantiles = c(0.025,0.975), M  = 5000, B = 500, data_len = 2000, burn_in = 5000, Nonlinear_fun = "log(X_0^2)", Dis_error = "rnorm(1)", n_core = 4){
  cl = makeCluster(n_core)    
  registerDoSNOW(cl)
  #data = data
  #Training_window = Training_window
  #T = fsteps
  #alpha = allalpha[aaa]
  iterations = Rep
  print(paste("repetition times = ", iterations))
  pb = txtProgressBar(max = iterations, style = 3)
  progress = function(j) setTxtProgressBar(pb, j)
  opts = list(progress = progress) 
  All_rep_results = foreach(iii=1:iterations, .combine="rbind",.options.snow = opts,.export=c("Boot_pre_interval","Simulate_data","colMedians")) %dopar%{
    Boot_pre_interval(Pre_step = Pre_step, quantiles = quantiles, M = M , B = B, residuals_type = residuals_type, data_len = data_len, burn_in = burn_in, Nonlinear_fun = Nonlinear_fun, Dis_error = Dis_error)
  }
  stopCluster(cl)
  return(All_rep_results)
}




##########################################################################################
# Do a single prediction simulation of the NLAR model in the format "0.2 + log(0.5 + abs(X_0))"
##########################################################################################

# With fitted residuals
test = Boot_pre_interval(Pre_step = 5, quantiles = c(0.025,0.975), M  = 5000, B = 500, data_len = 50, burn_in = 1000, Nonlinear_fun = "0.2 + log(0.5 + abs(X_0))", Dis_error = "rnorm(1)",residuals_type = "fitted")

# With predictive residuals
test_p = Boot_pre_interval(Pre_step = 5, quantiles = c(0.025,0.975), M  = 5000, B = 500, data_len = 50, burn_in = 1000, Nonlinear_fun = "0.2 + log(0.5 + abs(X_0))", Dis_error = "rnorm(1)",residuals_type = "predictive")


##########################################################################################
# Do a simulation of NLAR model in the format "0.2 + log(0.5 + abs(X_0))" with four replications
##########################################################################################


test_rep = Do_replication_parallel(Rep = 4,Pre_step = 5, ,residuals_type = "fitted", quantiles = c(0.025,0.975), M  = 5000, B = 1000, data_len = 50, burn_in = 1000, Nonlinear_fun = "0.2 + log(0.5 + abs(X_0))", Dis_error = "rnorm(1)",n_core = 4)

# Note: once we have a large number of replications, we can compare the performance of different types of PIs by considering the average performance among these replications.


