source("call_specific_simulation.R")
source("get_data_for_models.R")
source("run_model_analysis.R")
source("select_me_model_approach.R")

# Simulate Data
sim_data = call_specific_simulation(seed=1, ds="2a", tp=2, short_run=F, remove_low_case_study=TRUE)
saveRDS(sim_data, file="Example_DS2A_TP2_Long.RData")

# See simulated data
plot(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, type="l") # R
plot(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$newcase.adj, type="l") # I.star
plot(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$newcase, type="l") # I

# Run proposed model with measurement error model on, reporting periods known
jags_data_me_model = get_data_to_pass_to_jags(study_params=sim_data$sim_params, si_known=sim_data$scenario_params$serial_known, 
                                              me_model_on=T, me_model_known=TRUE) # me_model_known=FALSE would estimate reporting periods
proposed_results_me_model = run_proposed_model(ds_params=sim_data$scenario_params, 
                                               data_for_jags=jags_data_me_model, 
                                               params=c("R", "I.star", "betas.x", "thetas", "phis", "beta.R", "lambdas"), 
                                               burnin=2000, iter=4000, num.chains=8, thin=2, 
                                               measurement_on=TRUE, keep_convergence_data=FALSE)
plot(x=1:sim_data$sim_params$study.length, y=proposed_results_me_model$mean$R, type="l") # proposed with measurement error model
lines(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, col="red") # true

# Run proposed model with case smoothing, no measurement error model
jags_data_me_smooth = get_data_to_pass_to_jags(study_params=sim_data$sim_params, si_known=sim_data$scenario_params$serial_known, 
                                               me_model_on=F)
proposed_results_me_smooth = run_proposed_model(ds_params=sim_data$scenario_params, 
                                                data_for_jags=jags_data_me_model, 
                                                params=c("R", "I.star", "betas.x", "thetas", "phis", "beta.R", "lambdas"), 
                                                burnin=2000, iter=4000, num.chains=8, thin=2, 
                                                measurement_on=FALSE, keep_convergence_data=FALSE)
plot(x=1:sim_data$sim_params$study.length, y=proposed_results_me_smooth$mean$R, type="l") # proposed with case smoothing
lines(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, col="red") # true

# Check if model reporting pattern (TRUE) or smooth cases (FALSE)
select_me_model(proposed_results_me_model$mean$R, proposed_results_me_smooth$mean$R)

# Run baseline approaches
baseline_results = run_baseline_model(ds_params=sim_data$scenario_params,
                                      study_params=jags_data_me_model, 
                                      baseline_uncertain=FALSE, agg.windows=c(1,3,7)) 
plot(x=2:sim_data$sim_params$study.length, y=baseline_results$run_1$mean, type="l") # baseline W=1
lines(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, col="red") # true

plot(x=2:sim_data$sim_params$study.length, y=baseline_results$run_3$mean, type="l") # baseline W=3 (for evaluation we center these estimates)
lines(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, col="red") # true

plot(x=2:sim_data$sim_params$study.length, y=baseline_results$run_7$mean, type="l") # baseline W=7 (for evaluation we center these estimates)
lines(x=1:sim_data$sim_params$study.length, y=sim_data$sim_params$R, col="red") # true
