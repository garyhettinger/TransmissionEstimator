# Code to get specified parameters

#############################################################################################
## get_pre_simulation_params: get manually-specified parameters for simulation by trend pattern
## and data scenario
### Param @tp: specified trend pattern
### Param @ds_params: list of paramters specific to data scenario
### Param @jitter: indicator of if to jitter covariates or not
### Returns: list of specified simulation parameters
#############################################################################################

get_pre_simulation_params = function(tp, short_run, ds_params, jitter) {
  N_and_I0 = get_N_and_I0(short_run)
  N = N_and_I0$N
  I.0 = N_and_I0$I.0
  x_param_func = get_x_param_func(tp)
  x_params = x_param_func(N, I.0, jitter)
  serial_interval_params = get_serial_interval_params(ds_params$measurement_error_specification)
  measurement_error_params = get_measurement_error_params(ds_params$measurement_error_specification, N)
  transmission_params = get_transmission_params(ds_params$measurement_error_specification, 
    ds_params$autoregression_present)
  pre_sim_params = c(N_and_I0, x_params, serial_interval_params,
    measurement_error_params, transmission_params)
  return(pre_sim_params)
}


#############################################################################################
## Helper functions
#############################################################################################
get_N_and_I0 = function(short_run) {
  # Helper function returns the starting values for a given trend pattern
  if (short_run) {
    N = 50
    I.0 = 200
  }
  else {
    N = 100
    I.0 = 200
  }
  return(list(N=N, I.0=I.0))
}

get_x_param_func = function(tp) {
  # Helper function selects the function to generate the covariates for a given trend pattern
  if (tp == 1) {
    x_func = get_x_params_tp1
  }
  if (tp == 2) {
    x_func = get_x_params_tp2
  }
  if (tp == 3) {
    x_func = get_x_params_tp3
  }
  return(x_func)
}

get_x_params_tp1 = function(N, I.0, jitter) {
  # First simulation involves 3 periods: 
  # A period of high unrestricted growth
  # A period of lockdown and significant decay
  # A period of equilibrium near R=1
  t1 = N/5 
  t2 = 2*N/5
  t3 = N
  
  q = 2
  beta.0.x = 1.2106
  beta.x = 2.2397
  
  ### Phase 1: Initial no changes
  X1 = rep(-0.2, t1) # R corresponds to 2.15
  
  # Phase 2: Lockdown drops rate
  X2 = rep(-0.85, t2-t1) # R here corresponds to 0.5
  
  # Phase 3: Long term equilibrium
  X3 = rep(-beta.0.x / beta.x, t3-t2) # so R = 1.0
  
  X = c(X1, X2, X3)
  if (jitter) {
    X = jitter(X, amount=0.05)
  }
  X = cbind(rep(1,N), X)
  return(list(N=N, betas.x=c(beta.0.x,beta.x), q=q, X=X, I.0=I.0))
}

get_x_params_tp2 = function(N, I.0, jitter) {
  # Second simulation involves I.0 / 10 periods that have a 5 period repeat of: 
  # A period of high unrestricted growth
  # A period of lockdown and significant decay
  # A second period of high unrestricted growth
  # A second period of lockdown and significant decay
  # A period of equilibrium near R=1
  
  ### social distancing score of x = -0.05 -> R = 3, x = -0.85 -> R = 0.5 
  q = 2
  beta.0.x = 1.2106
  beta.x = 2.2397
  
  ### Goal: form I.0 / 10 phases
  P = (N/10)
  X = c()
  for (p in 1:P) {
    if ((p %% 5) == 0) {
      x.val = -beta.0.x / beta.x
    }
    else if ((p %% 2) == 1) {
      x.val = -0.2
    }
    else {
      x.val = -0.85
    }
    Xp = rep(x.val, 10)
    X = c(X, Xp)
  }
  if (jitter) {
    X = jitter(X, amount=0.05)
  }
  X = cbind(rep(1,N), X)
  return(list(N=N, betas.x=c(beta.0.x,beta.x), q=q, X=X, I.0=I.0))
}


get_x_params_tp3 = function(N, I.0, jitter) {
  # Third simulation involves one continuous period: 
  # A flipped sigmoid type curve starting around 3 and dropping to around 0.8
  t = 1:N
  
  ### social distancing score of x = -0.1 -> R = 2.5, x = -0.8 -> R = 0.7 
  q = 2
  beta.x = log(2.5/0.7)/0.7
  beta.0.x = log(2.5) + 0.1*beta.x
  
  X = -0.7/(1+exp(-.25*(t - 15))) - 0.1
  
  if (jitter) {
    X = jitter(X, amount=0.05)
  }
  X = cbind(rep(1,N), X)
  return(list(N=N, betas.x=c(beta.0.x,beta.x), q=q, X=X, I.0=I.0))
}

get_serial_interval_params = function(measurement_error_specification) {
  s.max = 4 # max number of days
  if (measurement_error_specification != 4) {
    # Create w from a mixture of studies
    K = 3
    w.1 = c(.8, .1, .075, .025)
    w.2 = c(.1, .4, .3, .2)
    w.3 = c(.05, .15, .15, .65)
    w.k = cbind(w.1, w.2, w.3)
    lambdas = c(1/10, 7/10, 2/10)
  }
  else {
    # Create w so only one day infection
    K = 1
    w.k = c(1, 0, 0, 0)
    lambdas = c(1)
  }
  return(list(w.k=w.k, K=K, lambdas=lambdas, s.max=s.max))
}

get_measurement_error_params = function(measurement_error_specification, N) {
  # Create thetas depending on specification
  t = 1:N
  tau = 3
  weekly = c(1,1,2,2,3,3,3)
  tau.index = ifelse(t %% 7 != 0, t %% 7, 7)
  tau.t = weekly[tau.index]
  if (measurement_error_specification == 0) {
    thetas = rep(1, tau)
  }
  else if (measurement_error_specification == 1) {
    thetas = c(0.5, 1.5, 1.0)
  }
  else if (measurement_error_specification == 2) { # changing thetas to null halfway
      theta.1 = 0.5
      theta.4 = 1.0
      
      theta.2 = 1.5
      theta.5 = 1.0
      
      theta.3 = 1.0
      theta.6 = 1.0
      
      tau.t = c(rep(0, N/2), rep(3, N/2)) + tau.t
      thetas = c(theta.1, theta.2, theta.3, theta.4, theta.5, theta.6)
    }
  else if (ds_params$measurement_error_specification == 3) { # Dirichlet violated
    thetas = c(0.5, 1.25, 0.8)
  }
  else { # data scenario 3C
    # Create thetas based on R
    theta.1 = 0.5
    R.1 = 2.15 # roughly R in first period
    theta.2 = (R.1 + theta.1) / R.1
    R.2 = 0.5 # roughly R in second period
    theta.3 = (R.2 + theta.1) / R.2
    R.3 = 1.0 # roughly R in third period
    theta.4 = (R.3 + theta.1) / R.3
    tau = 2
    tau.t = c(rep(c(1,2), N/10), rep(c(1,3), N/10), rep(c(1,4), 3*N/10))
    thetas = c(theta.1, theta.2, theta.3, theta.4)
  }
  return(list(thetas=thetas, tau=tau, tau.t=tau.t))
}

get_transmission_params = function(measurement_error_specification, autoregression_present) {
  beta.R = .2 # Rt dispersion parameter about its mean
  if (measurement_error_specification == 4) {
    beta.R = .05
  }
  ### Declare AR Portion of model
  m.phi = 2
  if (autoregression_present) {
    phis = c(.4, -1/6)
  }
  else {
    phis = rep(0, m.phi)
  }
  return(list(beta.R=beta.R, m.phi=m.phi, phis=phis))
}




