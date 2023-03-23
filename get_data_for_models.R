# Code to get data for models
source("estimate_reporting_periods.R")

#############################################################################################
## get_data_to_pass_to_jags(): creates observed data from simulated outbreak
### Param @study_params: list with relevant parameters for study
### Param @si_known: if TRUE, assumes serial interval is known and will not estimate weights
### Param @me_model_on: if TRUE, assumes no measurement error and will use I=I.star
### Param @me_model_known: if TRUE, assumes reporting periods are known and will not estimate
### Returns: list of data for jags model
#############################################################################################

get_data_to_pass_to_jags = function(study_params, si_known, me_model_on, me_model_known=FALSE) {
    # This function serves to aggregate data and priors to pass to model

    # Known data
    N = study_params$study.length
    m.phi = study_params$ar.order
    s.max = study_params$max.days.serial
    X = study_params$X
    q = ncol(X)
    
    hyperparams = get_noninformative_hyperparameters()
    
    data_for_jags = list(
      N=N,
      s.max=s.max,
      beta.R.scale=hyperparams$beta.R.scale, 
      beta.R.shape=hyperparams$beta.R.shape,
      beta.0.scale=hyperparams$beta.0.scale, 
      beta.0.shape=hyperparams$beta.0.shape, 
      betas.x.scale=hyperparams$betas.x.scale, 
      betas.x.shape=hyperparams$betas.x.shape, 
      m.phi=m.phi, 
      phis.shape=hyperparams$phis.shape,
      phis.scale=hyperparams$phis.scale,
      X=X, 
      q=q
    )

    me_model_data_for_jags = get_me_model_data_for_jags(me_model_on, me_model_known, study_params)
    si_data_for_jags = get_serial_interval_data_for_jags(si_known, study_params)
    data_for_jags = c(data_for_jags, me_model_data_for_jags, si_data_for_jags)
    return(data_for_jags)
}

get_theta_priors_me = function(county_data, study_length, dirichlet.conc=3, prior_var=1) {
  # Retrieves the priors for the measurement error model weights

  study_data = tail(county_data, study_length)
  me_model_info = get_me_model_info(study_data)
  tau.t = sapply(study_data$weekday, function(x) me_model_info$map.info$mapping[x])
  tau = max(me_model_info$map.info$mapping)
  priors = list(thetas.dirch.alpha=dirichlet.conc*me_model_info$map.info$prior,
    tau.t=tau.t, tau=tau)
  return(list(theta.priors=priors, me_model_info=me_model_info))
}

get_noninformative_hyperparameters = function() {
  return(list(beta.R.scale = .01,
    beta.R.shape = .01,
    beta.0.scale = 1/(10^2),
    beta.0.shape = 0,
    betas.x.scale = 1/(10^2),
    betas.x.shape = 0,
    phis.scale = 1/(10^2),
    phis.shape = 0))
}

get_me_model_data_for_jags = function(me_model_on, me_model_known, study_params) {
  I = study_params$newcase
  N = study_params$study.length
  # Measurement error model
  me_model_data_for_jags = list()
  if (me_model_on) {
    if (me_model_known) {
        tau = study_params$num.me.periods 
        tau.t = study_params$me.period.mappings
        theta.priors = list(thetas.dirch.alpha=rep(1, study_params$num.me.periods))
        theta.priors[["tau"]] = tau
        theta.priors[["tau.t"]] = tau.t
        me_model_data_for_jags = theta.priors
    }
    else {
        theta_prior_info = get_theta_priors_me(list(newcase=I, weekday=study_params$weekday), 
          N)
        I = theta_prior_info$me_model_info$outlier.data$outlier.free.cases
        me_model_data_for_jags = theta_prior_info$theta.priors
    }
  }
  me_model_data_for_jags$I = I
  return(me_model_data_for_jags)
} 

get_serial_interval_data_for_jags = function(si_known, study_params) {
  si_data_for_jags = list()
  if (si_known) {
    si_data_for_jags$w = study_params$w
  }
  else {
    si_data_for_jags$w.k = study_params$serial.intervals
    si_data_for_jags$K = ncol(si_data_for_jags$w.k)
    si_data_for_jags$lambda.dirch.alpha = if ("serial_interval_weights" %in% names(
      study_params)) 3*study_params$serial.interval.weights else rep(1, si_data_for_jags$K)
  }
  return(si_data_for_jags)
}


