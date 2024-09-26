

# Function to simulate all order 1 NLAR model; we just need to replace Nonlinear_fun = "sin(X_0)"  by any appropriate NLAR functions.

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

# Function to simulate all order 3 NLAR model; we just need to replace Nonlinear_fun = "log( 4*exp(0.9*X_0[3]) + 5*exp(0.9*X_0[2]) +6*exp(0.9*X_0[1]))"  by any appropriate NLAR functions.

Simulate_data_order3 = function(burn_in = 5000,data_len = 2000 , Nonlinear_fun = "log( 4*exp(0.9*X_0[3]) + 5*exp(0.9*X_0[2]) +6*exp(0.9*X_0[1]))" ,Dis_initial_x = "rnorm(1)",Dis_error = "rnorm(0,1)"){
  Simulated_X = c()
  Total_len = burn_in + data_len
  X_0 = c(eval(parse(text = Dis_initial_x)), eval(parse(text = Dis_initial_x)),eval(parse(text = Dis_initial_x)))
  for (i in c(1:Total_len)) {
    X_latest = eval(parse(text = Nonlinear_fun)) + eval(parse(text = Dis_error))
    Simulated_X = c(Simulated_X,X_latest)
    X_0 = c(X_0[2:3],X_latest)
  }
  return(Simulated_X)
}



# Function to simulate order(1,1) threshold model.

Sim_slef_threshold_1_1_r_0<- function(data_len, alpha, beta, r, burn_in=100)
{
  # Generate noise
  e = rnorm(data_len+burn_in, 0, 1)
  # Create space for y
  y = numeric(data_len+burn_in)
  # Generate time series
  y[1] = runif(1,min = -1,max = 1)
  for(i in 2:(data_len+burn_in))
  {
    if(y[i-1] <= r)
      y[i] = alpha*y[i-1] + e[i]
    else
      y[i] = beta*y[i-1] +  e[i]
  }
  # Throw away first burnin values
  y <- y[-(1:burn_in)]
  # Return result
  return(y)
}

# Function to estimate order(1,1) threshold model.

fitnlts <- function(x, y,r = 0)
{
  ss <- function(par,x = x, y = y)
  {
    alpha <- par[1]
    beta <- par[2]
    n <- length(x)
    e1 <- y - alpha*x
    e2 <- y - beta*x
    regime1 <- (x <= r)
    e <- e1*(regime1) + e2*(!regime1)
    return(sum(e^2))
  }
  fit <- optim(c(0,0),ss,x=x, y=y, control=list(maxit=1000))
  return(fit)
}


# Function to simulate order(3,1) threshold model.

Sim_slef_threshold_3_1_r_0<- function(data_len, alpha, beta, r, burn_in=100)
{
  # Generate noise
  e = rnorm(data_len+burn_in, 0, 1)
  # Create space for y
  y = numeric(data_len+burn_in)
  # Generate time series
  y[1:3] = runif(3,min = -1,max = 1)
  for(i in 4:(data_len+burn_in))
  {
    if(y[i-1] <= r)
      y[i] = alpha[1]*y[i-1] + alpha[2]*y[i-2]  + alpha[3]*y[i-3]  + e[i]
    else
      y[i] = beta*y[i-1] +  e[i]
  }
  # Throw away first burnin values
  y <- y[-(1:burn_in)]
  # Return result
  return(y)
}

# Function to estimate order(3,1) threshold model.

fitnlts_3_1 <- function(x, o,z, y,r = 0)
{
  ss <- function(par,x = x, o = o, z = z, y = y)
  {
    alpha_1 <- par[1]
    alpha_2 = par[2]
    alpha_3 = par[3]
    beta <- par[4]
    n <- length(x)
    e1 <- y - alpha_1*x - alpha_2*o - alpha_3*z
    e2 <- y - beta*x
    regime1 <- (x <= r)
    e <- e1*(regime1) + e2*(!regime1)
    return(sum(e^2))
  }
  fit <- optim(c(1,1,1,1),ss,x=x, o=o,z=z, y=y, control=list(maxit=1000))
  return(fit)
}

# Function to simulate order(3,1) threshold model with heteroscedastic error

Sim_slef_threshold_1_1_r_0_heter<- function(data_len, alpha, beta, heter = 0.5, r, burn_in=100)
{
  # Generate noise
  e = rnorm(data_len+burn_in, 0, 1)
  # Create space for y
  y = numeric(data_len+burn_in)
  # Generate time series
  y[1] = runif(1,min = -1,max = 1)
  for(i in 2:(data_len+burn_in))
  {
    if(y[i-1] <= r)
      y[i] = alpha*y[i-1] + e[i]*exp(-(y[i-1])^2)*heter
    else
      y[i] = beta*y[i-1] +  e[i]*exp(-(y[i-1])^2)*heter
  }
  # Throw away first burnin values
  y <- y[-(1:burn_in)]
  # Return result
  return(y)
}

# The estimation stage 1 of order(3,1) threshold model with heteroscedastic error

fitnlts_heter <- function(x, y,r = 0)
{
  ss <- function(par,x = x, y = y)
  {
    alpha <- par[1]
    beta <- par[2]
    n <- length(x)
    e1 <- y - alpha*x
    e2 <- y - beta*x
    regime1 <- (x <= r)
    e <- e1*(regime1) + e2*(!regime1)
    return(sum(e^2)/n)
  }
  fit <- optim(c(1,1),ss,x=x, y=y, control=list(maxit=1000))
  return(fit)
}

# The estimation stage 2 of order(3,1) threshold model with heteroscedastic error

fitnlts_heter_stage2 <- function(x, y, est_a, est_b, r = 0)
{
  ss <- function(par,x, y, est_a, est_b)
  {
    heter <- par[1]
    n <- length(x)
    e1 <- (y - est_a*x)/(heter*exp(-(x)^2))
    e2 <- (y - est_b*x)/(heter*exp(-(x)^2))
    regime1 <- (x <= r)
    e <- e1*(regime1) + e2*(!regime1)
    return(abs(sum(e^2)/n-1))
  }
  fit <- optim(c(1),ss,x=x, y=y, est_a = est_a, est_b = est_b, method = "Brent", lower = -100, upper = 100,  control=list(maxit=1000))
  return(fit)
}

