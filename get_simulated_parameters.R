# Code to generate simulations

#############################################################################################
## simulate_outbreak(): creates observed data from simulated outbreak
### Param @pre_sim_params: list with pre-specified parameters for simulation
### Param @ds_params: list with parameters to specify data scenario
### Returns: list of observed data, specified parameters, and 'unknown' true parameter values
#############################################################################################

simulate_outbreak = function(pre_sim_params, ds_params) {
  measurement_variation = (ds_params$measurement_error_specification > 0)
  sim_params = run_simulation(pre_sim_params, measurement_variation, ds_params$parametric_misspecification)

  # Observe white noise covariate if predictors not specified
  sim_params = get_observed_predictor(ds_params$predictor_specified, sim_params)

  # Adjust observed measurement error periods to constant and weekly if not
  sim_params = get_observed_measurement_error(ds_params$measurement_error_specification, sim_params)
  
  return(sim_params)
}

#############################################################################################
## run_simulation(): adds simulated time-series data (R, I.star, I) to specified parameters
### Param @pre_sim_params: list with pre-specified parameters for simulation
### Param @measurement_variation: indicator of if to simulate reporting variation
### Param @parametric misspecification: indicator of which parametric specification to use
### Returns: list of observed data, specified parameters, and 'unknown' true parameter values
#############################################################################################

run_simulation = function(pre_sim_params, measurement_variation, parametric_misspecification) {
  # These variables are set manually according to simulation specs
  N = pre_sim_params$N
  X = pre_sim_params$X
  betas.x = pre_sim_params$betas.x
  q = pre_sim_params$q
  m.phi = pre_sim_params$m.phi
  phis = pre_sim_params$phis
  beta.R = pre_sim_params$beta.R
  s.max = pre_sim_params$s.max
  w.k = pre_sim_params$w.k
  lambdas = pre_sim_params$lambdas
  tau = pre_sim_params$tau
  tau.t = pre_sim_params$tau.t
  thetas = pre_sim_params$thetas
  I.0 = pre_sim_params$I.0
  K = pre_sim_params$K

  # Set distribution functions for correct/incorrect specifications
  if (!parametric_misspecification) {
    distributions = list(R=function(mean, scale) rgamma(n=1, shape=mean/scale, scale=scale),
                         I.star=function(mean) rpois(n=1, lambda=mean), 
                         I=function(mean) rpois(n=1, lambda=mean))
  }
  else {
    distributions = list(R=function(mean, scale) rlnorm(n=1, meanlog=(log(mean) - 1/2*scale^2), sdlog=scale),
                         I.star=function(mean, scale=10) rnbinom(n=1, mu=mean, size=scale), 
                         I=function(mean, scale=0.1) round(rlnorm(n=1, meanlog=(log(mean) - 1/2*scale^2), sdlog=scale)))
  }
  
  R = simulate_R(N, X, betas.x, beta.R, phis, m.phi, distributions)
  w = if (K > 1) c(w.k %*% lambdas) else w.k # Specify true serial interval
  I.star = simulate_true_I(N, I.0, thetas, tau.t, R, w, s.max, measurement_variation, distributions)
  I = simulate_I(N, I.0, thetas, tau.t, I.star, measurement_variation, distributions)
  
  return(list(
    study.length=N,
    X=X,
    betas.x=betas.x,
    beta.0.x=betas.x[1],
    beta.x=betas.x[2],
    ar.order=m.phi,
    phis=phis, 
    beta.R=beta.R, 
    max.days.serial=s.max,
    serial.intervals=w.k, 
    serial.interval.weights=lambdas, 
    num.me.periods=tau, 
    me.period.mappings=tau.t,
    weekday=sapply(1:N, function(x) if (x %% 7 > 0) x %% 7 else 7),
    thetas=thetas, 
    newcase=I,
    newcase.adj = I.star,
    R = R,
    w = w,
    w.k = w.k,
    date = seq.Date(as.Date('2020-01-01'),(as.Date('2020-01-01')+N-1),by=1)
  ))
}

#############################################################################################
## Helper Functions
#############################################################################################

simulate_R = function(N, X, betas.x, beta.R, phis, m.phi, distributions) {
  # Simulate R
  R = rep(0.0, N)
  R.means = rep(0.0, N)
  R.ar = rep(0.0, N)
  R.regs = X %*% betas.x # R from predictors only
  R.ar[1] = 0
  R.means[1] = exp(R.regs[1])  # baseline R
  R[1] = distributions$R(mean=R.means[1], scale=beta.R)
  for (i in 2:N) {
    R.ar[i] = t(phis[1:min(i-1, m.phi)]) %*% rev(log(R[max(i-m.phi,1):(i-1)]) - R.regs[max(i-m.phi, 1):(i-1)])
    R.means[i] = exp(R.regs[i] + R.ar[i])
    R[i] = distributions$R(mean=R.means[i], scale = beta.R)
  }
  return(R)
}

simulate_true_I = function(N, I.0, thetas, tau.t, R, w, s.max, measurement_variation, distributions) {
  # Simulate I.star
  I.star = rep(0.0, N)
  I.star.mean = rep(0.0, N)
  I.star.mean[1] = I.0 / thetas[tau.t[1]]
  if (measurement_variation) {
    I.star[1] = pmax(1, distributions$I.star(mean=I.star.mean[1])) # baseline number of cases based off observed (for simulation purposes)
  }
  else {
    I.star[1] = I.star.mean[1]
  }
  for (i in 2:N) {
    I.star.mean[i] = R[i] * (t(w[1:min(i-1, s.max)]) %*% rev(I.star[max(i-s.max, 1):(i-1)]))
    I.star[i] = pmax(1, distributions$I.star(mean=I.star.mean[i]))
  }
  return(I.star)
}

simulate_I = function(N, I.0, thetas, tau.t, I.star, measurement_variation, distributions) {
  ### Simulate I
  I = rep(0, N)
  I[1] = I.0
  for (i in 2:N) {
    if (measurement_variation) {
      I[i] = pmax(1, distributions$I(mean=(I.star[i] * thetas[tau.t[i]])))
    }
    else {
      I[i] = I.star[i] * thetas[tau.t[i]]
    }
  }
  return(I)
}


get_observed_predictor = function(predictor_specified, sim_params) {
  if (!predictor_specified) {
    known_X = cbind(rep(1,sim_params$study.length), runif(sim_params$study.length, -.1, .1))
    known_beta.0.x = 0
    known_beta.x = 0
    sim_params$X = known_X
    sim_params$beta.0.x = known_beta.0.x
    sim_params$beta.x = known_beta.x
  }
  return(sim_params)
}

get_observed_measurement_error = function(measurement_error_specification, sim_params) {
  if (measurement_error_specification > 1) {
    t = 1:sim_params$study.length
    curr_thetas = sim_params$thetas[sim_params$me.period.mappings]
    if (measurement_error_specification == 4) { # in this case, tau's are not weekly
      sim_params$me.period.mappings = ifelse(t %% sim_params$num.me.periods != 0, 
                                t %% sim_params$num.me.periods, sim_params$num.me.periods)
    }
    else {
      weekly = c(1,1,2,2,3,3,3)
      tau.index = ifelse(t %% 7 != 0, t %% 7, 7)
      sim_params$me.period.mappings = weekly[tau.index]
    }
    known_thetas = rep(1, sim_params$num.me.periods)
    for (i in 1:sim_params$num.me.periods) {
      # pass back mean thetas if want to compare to
      known_thetas[i] = mean(curr_thetas[which(sim_params$me.period.mappings == i)]) 
    }
    sim_params$thetas = known_thetas
  }
  return(sim_params)
}
  


