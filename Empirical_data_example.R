require(doSNOW)
require(parallel)
require(foreach)
require(doParallel)
require(astsa) # it contains the used UnempRate data

# Bootstrap/simulation prediction /Empirical analysis
Boot_pre_Unemploy = function(Pre_step = 5, quantiles = c(0.025,0.975),  M  = 500, B = 1000, Y = NULL){
  
  
  # Pre_Pre_interval = matrix(ncol = Pre_step, nrow = 19) # 1st row We return the true L_1/L_2 optimal furture value by simulation
  # rownames(Pre_Pre_interval) = c("True simulated future","L_2 Bootstrap point fitres", "L_1 Bootstrap point fitres","L_2 Bootstrap point preres", "L_1 Bootstrap point preres", "Low quantile PI fitres", "High quantile PI fitres", "Low quantile PI preres", "High quantile PI preres", "Low L_2 Pertinent PI fitres",  "High L_2 Pertinent PI fitres", "Low L_1 Pertinent PI fitres",  "High L_1 Pertinent PI fitres", "Low L_2 Pertinent PI preres", "High L_2 Pertinent PI preres","Low L_1 Pertinent PI preres", "High L_1 Pertinent PI preres","Low True sim PI","High True sim PI" )
  # # 2-4 rows we return the qunatile Pre_interval, and its corresponding L2/L1 predictor
  # # 5-6 rows we return the centered Pre_interval, and its corresponding center
  # # Simulate data
  
  Boot_mat = matrix(ncol = Pre_step, nrow = M)
  Boot_mat_pre = matrix(ncol = Pre_step, nrow = M) # Bootstrap prediction with predictive residuals
  
  # Find quantile PI and corresponding center value by bootstrap
  # Do estimation first
  #diff_y = diff(data,1)
  Y_n = length(Y)
  y = Y[3:Y_n]
  y_lag1 = Y[2:(Y_n-1)]
  y_lag2 = Y[1:(Y_n - 2)]
  
  data_x = matrix(ncol = 2, nrow = length(y))
  data_x[,1] = exp(-y_lag1^2)*y_lag1
  data_x[,2] = y_lag2
  data_x = as.data.frame(data_x)
  colnames(data_x) = c("V1","V2")
  
  PPI_bootstrap = vector("list", length = 4) 
  m_o = 2
  lm_reg <- lm(y ~ V1 + V2, data = data_x )
  
  
  fit_residuals = lm_reg$residuals
  fit_residuals = fit_residuals - mean(fit_residuals)
  
  predictive_residual = numeric(length(fit_residuals))
  
  n = length(y)
  
  
  for(i in c((m_o+1):n)){
    lm_reg_d = lm(y[-i] ~ V1 + V2  , data = data_x[-i,] )
    predictive_residual[i] = y[i] - predict(lm_reg_d, newdata = data.frame(data_x[i,]))
    
  }
  
  predictive_residual = predictive_residual - mean(predictive_residual)
  
  # Do prediction naive method within package
  
  
  
  Naive_pre_Interval = matrix(NA,ncol = 2,nrow = Pre_step)
  Naive_pre_point = vector(length = 5)
  
  #new_data = data_x[n,]
  
  latest_data = data.frame(c( exp(-y[n]^2)*y[n] , data_x[n,][1] )  )
  colnames(latest_data) = c("V1","V2")
  
  new_data = latest_data
  for (i in c(1:Pre_step)){
    Naive_pre = predict(lm_reg, newdata = new_data, se.fit = TRUE, interval = c("prediction"))
    Naive_pre_point[i] = Naive_pre$fit[1]
    Naive_pre_Interval[i,1] = Naive_pre$fit[2]
    Naive_pre_Interval[i,2] = Naive_pre$fit[3]
    new_data[2] = new_data[1]
    new_data[1] = exp(-Naive_pre$fit[1]^2)*Naive_pre$fit[1]
  }
  
  
  # Do pertinent prediction interval; We do bootstrap B times to compute predictive root
  
  # First compute the optimal point prediction in the real-world by bootstrap
  
  Pre_mat_fit = matrix(ncol = Pre_step, nrow = M)
  Pre_mat_pre = matrix(ncol = Pre_step, nrow = M)
  
  #Pre_real = numeric(m_o+Pre_step)
  coef = lm_reg$coefficients[1:(m_o+1)]
  
  
  
  # Do Bootstrap prediction in the bootstrap world
  for (m in c(1:M)){
    new_data = latest_data
    for (h in c(1:Pre_step)) {
      
      Pre_real_single = sum(coef[2:(m_o+1)] * new_data) + coef[1] + sample(fit_residuals,1,replace = T)
      Pre_mat_fit[m,h] = Pre_real_single
      new_data[2] = new_data[1]
      new_data[1] = exp(-Pre_real_single^2)*Pre_real_single
      
      
    }
  }
  
  
  
  # Do Bootstrap prediction in the bootstrap world
  for (m in c(1:M)) {
    
    new_data = latest_data
    
    for (h in c(1:Pre_step)) {
      
      Pre_real_single = sum(coef[2:(m_o+1)] * new_data) + coef[1] + sample(predictive_residual,1,replace = T)
      Pre_mat_pre[m,h] = Pre_real_single
      new_data[2] = new_data[1]
      new_data[1] = exp(-Pre_real_single^2)*Pre_real_single
      
      
      
      
    }
  }
  
  
  
  PPI_bootstrap[[1]] = colMeans(Pre_mat_fit)
  PPI_bootstrap[[2]] = colMeans(Pre_mat_pre)
  
  
  
  Boot_predictive_root_L2 = matrix(nrow = B, ncol = Pre_step)
  #Boot_predictive_root_L1 = matrix(nrow = B, ncol = Pre_step)
  
  for( b in c(1:B)){
    
    
    # Find the '' TRUE'' value in the bootstrap world
    True_bootstrap = numeric(Pre_step)
    new_data = latest_data
    
    for (h in c(1:Pre_step)) {
      
      X_0 = sum(coef[2:(m_o+1)] * new_data)  + coef[1] + sample(fit_residuals,1,replace = T)
      True_bootstrap[h] = X_0
      new_data[2] = new_data[1]
      new_data[1] = exp(-X_0^2)*X_0
      
    }
    
    
    
    
    
    
    Boot_simulated_Y = numeric(Y_n+m_o)
    loc = sample(c(1:(Y_n-m_o+1)),1,replace = T)
    Boot_simulated_Y[1:m_o] = Y[loc:(loc+m_o-1)]
    
    
    #Boot_simulated_x_pre = c()
    for (i in c(1:Y_n)){
      
      X_0 = sum(coef[2:(m_o+1)] * c(exp(-Boot_simulated_Y[i+m_o-1]^2)*Boot_simulated_Y[i+m_o-1], Boot_simulated_Y[i]) ) + coef[1] + sample(fit_residuals,1,replace = T)
      Boot_simulated_Y[i+m_o] = X_0
      
      
    }
    
    Boot_simulated_Y = Boot_simulated_Y[(m_o+1):length(Boot_simulated_Y)]
    
    
    Boot_y = Boot_simulated_Y[3:Y_n]
    Boot_y_lag1 = Boot_simulated_Y[2:(Y_n-1)]
    Boot_y_lag2 = Boot_simulated_Y[1:(Y_n-2)]
    
    Boot_data_x = matrix(ncol = 2, nrow = length(Boot_y))
    Boot_data_x[,1] = exp(-Boot_y_lag1^2)*Boot_y_lag1
    Boot_data_x[,2] = Boot_y_lag2
    Boot_data_x = as.data.frame(Boot_data_x)
    colnames(Boot_data_x) = c("V1","V2")
    
    lm_res_boot <- try(lm(Boot_y ~ V1 + V2  , data = Boot_data_x ), silent = TRUE)
    
    
    while(try(summary(lm_res_boot)[2], silent = TRUE) == "try-error" | if(try(summary(lm_res_boot)[2], silent = TRUE) != "try-error"){sum(is.na(lm_res_boot$coefficients))>=1}else{TRUE} ){
      
      Boot_simulated_Y = numeric(Y_n+m_o)
      loc = sample(c(1:(Y_n-m_o+1)),1,replace = T)
      Boot_simulated_Y[1:m_o] = Y[loc:(loc+m_o-1)]
      
      
      #Boot_simulated_x_pre = c()
      for (i in c(1:Y_n)){
        
        X_0 = sum(coef[2:(m_o+1)] * c(exp(-Boot_simulated_Y[i+m_o-1]^2)*Boot_simulated_Y[i+m_o-1], Boot_simulated_Y[i]) ) + coef[1] + sample(fit_residuals,1,replace = T)
        Boot_simulated_Y[i+m_o] = X_0
        
        
      }
      
      Boot_simulated_Y = Boot_simulated_Y[(m_o+1):length(Boot_simulated_Y)]
      
      
      Boot_y = Boot_simulated_Y[3:Y_n]
      Boot_y_lag1 = Boot_simulated_Y[2:(Y_n-1)]
      Boot_y_lag2 = Boot_simulated_Y[1:(Y_n-2)]
      
      Boot_data_x = matrix(ncol = 2, nrow = length(Boot_y))
      Boot_data_x[,1] = exp(-Boot_y_lag1^2)*Boot_y_lag1
      Boot_data_x[,2] = Boot_y_lag2
      Boot_data_x = as.data.frame(Boot_data_x)
      colnames(Boot_data_x) = c("V1","V2")
      
      lm_res_boot <- try(lm(Boot_y ~ V1 + V2  , data = Boot_data_x ), silent = TRUE)
      
    }
    
    # Find the bootstrap estimation in the bootstrap world
    Boot_Boot_mat = matrix(ncol = Pre_step, nrow = M)
    #n = length(Boot_simulated_x)
    
    
    coef_boot = lm_res_boot$coefficients[1:(m_o+1)]
    
    # Do Bootstrap prediction in the bootstrap world
    for (m in c(1:M)){
      new_data = latest_data
      for (h in c(1:Pre_step)) {
        
        Pre_real_single = sum(coef_boot[2:(m_o+1)] * new_data) + coef_boot[1] + sample(fit_residuals,1,replace = T)
        Boot_Boot_mat[m,h] = Pre_real_single
        new_data[2] = new_data[1]
        new_data[1] = exp(-Pre_real_single^2)*Pre_real_single
        
        
      }
    }
    
    
    Boot_boot_pre_L2 = colMeans(Boot_Boot_mat)
    #print(Boot_boot_pre_L2)
    #Boot_boot_pre_L1 = colMedians(Boot_Boot_mat)
    
    Boot_predictive_root_L2[b,] = True_bootstrap - Boot_boot_pre_L2
    #Boot_predictive_root_L1[b,] = True_bootstrap - Boot_boot_pre_L1
    
  }
  
  
  
  # Get the pertinent prediction interval
  PR_quantile_L2 = apply(Boot_predictive_root_L2, 2, function(x){quantile(x,probs = c(quantiles[1],quantiles[2]))})
  
  PPI_low_fit = PR_quantile_L2[1,] + PPI_bootstrap[[1]]
  PPI_up_fit = PR_quantile_L2[2,] + PPI_bootstrap[[1]]
  
  
  #
  # Apply the predictive residuals #############################################
  # #############################################
  # #############################################
  #############################################
  
  
  Boot_predictive_root_L2 = matrix(nrow = B, ncol = Pre_step)
  #Boot_predictive_root_L1 = matrix(nrow = B, ncol = Pre_step)
  
  for( b in c(1:B)){
    
    
    # Find the '' TRUE'' value in the bootstrap world
    True_bootstrap = numeric(Pre_step)
    new_data = latest_data
    
    for (h in c(1:Pre_step)) {
      
      X_0 = sum(coef[2:(m_o+1)] * new_data)  + coef[1] + sample(predictive_residual,1,replace = T)
      True_bootstrap[h] = X_0
      new_data[2] = new_data[1]
      new_data[1] = exp(-X_0^2)*X_0
      
    }
    
    
    
    
    
    
    Boot_simulated_Y = numeric(Y_n+m_o)
    loc = sample(c(1:(Y_n-m_o+1)),1,replace = T)
    Boot_simulated_Y[1:m_o] = Y[loc:(loc+m_o-1)]
    
    
    #Boot_simulated_x_pre = c()
    for (i in c(1:Y_n)){
      
      X_0 = sum(coef[2:(m_o+1)] * c(exp(-Boot_simulated_Y[i+m_o-1]^2)*Boot_simulated_Y[i+m_o-1], Boot_simulated_Y[i]) ) + coef[1] + sample(predictive_residual,1,replace = T)
      Boot_simulated_Y[i+m_o] = X_0
      
      
    }
    
    Boot_simulated_Y = Boot_simulated_Y[(m_o+1):length(Boot_simulated_Y)]
    
    
    Boot_y = Boot_simulated_Y[3:Y_n]
    Boot_y_lag1 = Boot_simulated_Y[2:(Y_n-1)]
    Boot_y_lag2 = Boot_simulated_Y[1:(Y_n-2)]
    
    Boot_data_x = matrix(ncol = 2, nrow = length(Boot_y))
    Boot_data_x[,1] = exp(-Boot_y_lag1^2)*Boot_y_lag1
    Boot_data_x[,2] = Boot_y_lag2
    Boot_data_x = as.data.frame(Boot_data_x)
    colnames(Boot_data_x) = c("V1","V2")
    
    lm_res_boot <- try(lm(Boot_y ~ V1 + V2  , data = Boot_data_x ), silent = TRUE)
    
    
    while(try(summary(lm_res_boot)[2], silent = TRUE) == "try-error" | if(try(summary(lm_res_boot)[2], silent = TRUE) != "try-error"){sum(is.na(lm_res_boot$coefficients))>=1}else{TRUE} ){
      
      Boot_simulated_Y = numeric(Y_n+m_o)
      loc = sample(c(1:(Y_n-m_o+1)),1,replace = T)
      Boot_simulated_Y[1:m_o] = Y[loc:(loc+m_o-1)]
      
      
      #Boot_simulated_x_pre = c()
      for (i in c(1:Y_n)){
        
        X_0 = sum(coef[2:(m_o+1)] * c(exp(-Boot_simulated_Y[i+m_o-1]^2)*Boot_simulated_Y[i+m_o-1], Boot_simulated_Y[i]) ) + coef[1] + sample(predictive_residual,1,replace = T)
        Boot_simulated_Y[i+m_o] = X_0
        
        
      }
      
      Boot_simulated_Y = Boot_simulated_Y[(m_o+1):length(Boot_simulated_Y)]
      
      
      Boot_y = Boot_simulated_Y[3:Y_n]
      Boot_y_lag1 = Boot_simulated_Y[2:(Y_n-1)]
      Boot_y_lag2 = Boot_simulated_Y[1:(Y_n-2)]
      
      Boot_data_x = matrix(ncol = 2, nrow = length(Boot_y))
      Boot_data_x[,1] = exp(-Boot_y_lag1^2)*Boot_y_lag1
      Boot_data_x[,2] = Boot_y_lag2
      Boot_data_x = as.data.frame(Boot_data_x)
      colnames(Boot_data_x) = c("V1","V2")
      
      lm_res_boot <- try(lm(Boot_y ~ V1 + V2  , data = Boot_data_x ), silent = TRUE)
      
    }
    
    # Find the bootstrap estimation in the bootstrap world
    Boot_Boot_mat = matrix(ncol = Pre_step, nrow = M)
    #n = length(Boot_simulated_x)
    
    
    coef_boot = lm_res_boot$coefficients[1:(m_o+1)]
    
    # Do Bootstrap prediction in the bootstrap world
    for (m in c(1:M)){
      new_data = latest_data
      for (h in c(1:Pre_step)) {
        
        Pre_real_single = sum(coef_boot[2:(m_o+1)] * new_data) + coef_boot[1] + sample(predictive_residual,1,replace = T)
        Boot_Boot_mat[m,h] = Pre_real_single
        new_data[2] = new_data[1]
        new_data[1] = exp(-Pre_real_single^2)*Pre_real_single
        
        
      }
    }
    
    
    Boot_boot_pre_L2 = colMeans(Boot_Boot_mat)
    #print(Boot_boot_pre_L2)
    #Boot_boot_pre_L1 = colMedians(Boot_Boot_mat)
    
    Boot_predictive_root_L2[b,] = True_bootstrap - Boot_boot_pre_L2
    #Boot_predictive_root_L1[b,] = True_bootstrap - Boot_boot_pre_L1
    
  }
  
  
  # Get the pertinent prediction interval
  PR_quantile_L2 = apply(Boot_predictive_root_L2, 2, function(x){quantile(x,probs = c(quantiles[1],quantiles[2]))})
  
  PPI_low_pre = PR_quantile_L2[1,] + PPI_bootstrap[[2]]
  PPI_up_pre = PR_quantile_L2[2,] + PPI_bootstrap[[2]]
  
  PPI_fit = list("PPI_bootstrap" = PPI_bootstrap[[1]], "PPI_low_fit" = PPI_low_fit, "PPI_up_fit" = PPI_up_fit )
  PPI_pre = list("PPI_bootstrap" = PPI_bootstrap[[2]], "PPI_low_pre" = PPI_low_pre, "PPI_up_pre" = PPI_up_pre )
  
  return(list("Naive_pre_point" = Naive_pre_point,"Naive_pre_Interval" = Naive_pre_Interval, "PPI_fit" = PPI_fit, "PPI_pre" = PPI_pre ))
}





### Do rolling window prediction


Run_rolling_prediction_Unemploy = function(Training_window = 100, data_all = NULL, H = 5,quantiles = c(0.025,0.975), M  = 500,B = 1000  ){
  L = length(data_all)
  Rep = L -  Training_window - H + 1

  
  All_coverage_rate_PI_naive = matrix(NA,ncol = 5,nrow = Rep)
  All_coverage_rate_PPI_fit = matrix(NA,ncol = 5,nrow = Rep)
  All_coverage_rate_PPI_pre = matrix(NA,ncol = 5,nrow = Rep)
  
  All_average_PI_naive_length = matrix(NA,ncol = 5,nrow = Rep)
  All_average_PPI_fit_length = matrix(NA,ncol = 5,nrow = Rep)
  All_average_PPI_pre_length = matrix(NA,ncol = 5,nrow = Rep)
  
  for (i in c(1:Rep)){
    print(paste("Do the",i,"th rolling window prediction"))
    train_data = data_all[i:(i+Training_window-1)]
    Pre_res = Boot_pre_Unemploy(Pre_step = H, quantiles = quantiles,  M  = M, B = B, Y = train_data)
    Obes_future =  data_all[(i+Training_window):(i+Training_window+H-1)]
    
    All_coverage_rate_PI_naive[i,] = Obes_future > Pre_res$Naive_pre_Interval[,1] & Obes_future < Pre_res$Naive_pre_Interval[,2]
    All_coverage_rate_PPI_fit[i,] = Obes_future > Pre_res$PPI_fit$PPI_low_fit & Obes_future < Pre_res$PPI_fit$PPI_up_fit
    All_coverage_rate_PPI_pre[i,] = Obes_future > Pre_res$PPI_pre$PPI_low_pre & Obes_future < Pre_res$PPI_pre$PPI_up_pre
    
    All_average_PI_naive_length[i,] = Pre_res$Naive_pre_Interval[,2] - Pre_res$Naive_pre_Interval[,1]
    All_average_PPI_fit_length[i,] = Pre_res$PPI_fit$PPI_up_fit - Pre_res$PPI_fit$PPI_low_fit 
    All_average_PPI_pre_length[i,] = Pre_res$PPI_pre$PPI_up_pre - Pre_res$PPI_pre$PPI_low_pre 
    
  }
  
  return(list("All_coverage_rate_PI_naive" = All_coverage_rate_PI_naive, "All_coverage_rate_PPI_fit" =  All_coverage_rate_PPI_fit,"All_coverage_rate_PPI_pre" = All_coverage_rate_PPI_pre,
              "All_average_PI_naive_length" = All_average_PI_naive_length, "All_average_PPI_fit_length" = All_average_PPI_fit_length, "All_average_PPI_pre_length" = All_average_PPI_pre_length  ))
  
}


# Do data manipulation
y_u = UnempRate[1:372]
y_u_q = vector(length = 124)
for (i in c(1:124)) {
  y_u_q[i] = mean(y_u[((i-1)*3+1):(i*3)])
}
log_linear_detrent <- detrend(log(y_u_q), 5)
plot(log_linear_detrent,type="b",col = 1,cex = 0.7)


# Do rolling window prediction with window size 50
Unemploy_EAR = Run_rolling_prediction_Unemploy(Training_window = 50, data_all = log_linear_detrent, H = 5,quantiles = c(0.025,0.975), M  = 200,B = 1000  )


# Note: the result Unemploy_EAR contains "All_coverage_rate_PI_naive"  "All_coverage_rate_PPI_fit"   "All_coverage_rate_PPI_pre"   "All_average_PI_naive_length" "All_average_PPI_fit_length" 
# "All_average_PPI_pre_length" for all rolling window predictions. Then, it is easy to apply the function colMeans to get the average performance of each type of PIs. 




