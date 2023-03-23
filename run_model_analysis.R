# Code to analyze data

#############################################################################################
## run_proposed_model: run model and collect results for proposed approach
### Param @ds_params: scenario specific parameter list
### Param @data_for_jags: list of data to pass to jags models
### Param @params: vector of parameters to keep results for
### Param @burnin: number of burnin iterations
### Param @iter: number of total iterations (including burnin)
### Param @num.chains: number of MCMC chains
### Param @thin: thinning parameter
### Param @measurement_on: boolean to include measurement error model in JAGS model
### Param @keep_convergence_data: boolean to also store data on convergence
### Returns: list of resulting data
#############################################################################################

run_proposed_model = function(ds_params, data_for_jags, params, burnin, iter, num.chains, thin, 
                         measurement_on=TRUE, keep_convergence_data=TRUE) {
  model_data = get_proposed_model_file(ds_params$serial_known, measurement_on, params)
  params = model_data$params
  model_file = model_data$model
  jags_model = R2jags::jags(data=data_for_jags, 
                    inits = function() generate_initial_vals(data_for_jags, measurement_on), 
                    parameters.to.save=params, model.file=model_file, 
                    n.burnin=burnin, n.iter=iter, n.chains=num.chains, n.thin=thin)
  jags.cis = list()
  R.list = jags_model$BUGSoutput$sims.list$R
  for (param in params) {
    jags.cis[[param]] = get_jags_CIs(jags_model$BUGSoutput$sims.list[[param]])
  }
  jags.data.list = list(mean=jags_model$BUGSoutput$mean, CI=jags.cis, tau.t=data_for_jags$tau.t)
  if (keep_convergence_data) {
    jags.data.list[["convergence"]] = jags_model$BUGSoutput$summary[,"Rhat"]
  }
  return(jags.data.list)
}

#############################################################################################
## run_baseline_model: run model and collect results for baseline approach
### Param @ds_params: scenario specific parameter list
### Param @study_params: list of data for study to pass to baseline model
### Param @baseline_uncertain: boolean to use uncertain_si or non_parametric_si approach
### Param @agg.windows: vector of smoothing window (bandwidth) parameters to run model with
### Returns: list of resulting data
#############################################################################################
run_baseline_model = function(ds_params, study_params, baseline_uncertain, agg.windows=c(1,3,7)) {
  # This function calls the baseline model across a few possible time windows

  baseline.data.list = list()
  baseline.w = if(ds_params$serial_known) study_params$w else rowMeans(study_params$w.k)
  baseline.w = c(0, baseline.w)
  for (agg in agg.windows) {
    if (baseline_uncertain) {
      mean_si = t(baseline.w) %*% (0:study_params$s.max)
      diff_mean_si = min(mean_si - 1, study_params$s.max - mean_si - 1)
      baseline.data = EpiEstim::estimate_R(study_params$I,
                             method = "uncertain_si",
                             config=EpiEstim::make_config(
                               mean_si = mean_si, std_mean_si = 0.7,
                               min_mean_si = mean_si - diff_mean_si, max_mean_si = mean_si + diff_mean_si,
                               std_si = 0.55, std_std_si = 0.25,
                               min_std_si = 0.1, max_std_si = 1,
                               t_start = c(rep(2, agg), 3:(study_params$N - agg + 1)),
                               t_end = 2:study_params$N
                             ))
    }
    else {
      baseline.data = EpiEstim::estimate_R(study_params$I,
                             method = "non_parametric_si",
                             config=EpiEstim::make_config(
                               si_distr = baseline.w,
                               t_start = c(rep(2, agg), 3:(study_params$N - agg + 1)),
                               t_end = 2:study_params$N
                             ))
    }
    baseline.data.list[[paste0("run_", agg)]] = list(mean=baseline.data$R[,3],
                                                 lower=baseline.data$R[,5],
                                                 upper=baseline.data$R[,11])
  }
  return(baseline.data.list)
}

get_jags_CIs = function(param_sims) {
  # Helper function to retrieve Confidence interval from JAGS
  return(apply(param_sims, 2, quantile, probs=c(0.025, 0.975)))
}

generate_initial_vals = function(data_for_jags, measurement_on=TRUE) {
  # This function is called by jags with each chain to initialize a new
  # set of starting values

  lambdas = rep(1/data_for_jags$K, data_for_jags$K)
  beta.R = runif(1, 0, 2)
  betas.x = runif(data_for_jags$q, -1, 1)
  phis = runif(data_for_jags$m.phi, -1, 1)
  R = runif(data_for_jags$N, 0, 5)
  I.star = round(pmax(1, runif(data_for_jags$N, -100, 100) + data_for_jags$I))
  vals_list = list(beta.R = beta.R,
                   betas.x = betas.x,
                   phis = phis,
                   R = R,
                   I.star = I.star)
  if (measurement_on) {
    thetas.dirch = rep(1/data_for_jags$tau, data_for_jags$tau)
    vals_list$thetas.dirch = thetas.dirch
  }
  if ("lambda.dirch.alpha" %in% names(data_for_jags)) vals_list$lambdas = lambdas
  return(vals_list)
}

# Select proposed model
get_proposed_model_file = function(serial_known, measurement_on, params) {
  # The actual call to run the proposed model in JAGS
  serial_tag = if (serial_known) "_noserial" else ""
  measurement_tag = if (measurement_on) "" else "_nomeasurement"
  if (!measurement_on) {
    params = params[(params != "I.star") & (params != "thetas")]
  }
  model = paste0("proposed_model", serial_tag, measurement_tag, ".bug")
  return(list(model=model, params=params))
}




