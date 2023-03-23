# Code to call specific simulation
source("get_specified_parameters.R")
source("get_simulated_parameters.R")

#############################################################################################
## call_specific_simulation: function to get data for a specific simulation
### Param @seed: random seed 
### Param @ds: data scenario from associated paper ("0", "1a", "1b", "1c", "2a", "2b", "2c", "3a", "3b", "3c")
### Param @tp: specified trend pattern from associated paper (1,2,3)
### Param @short_run: If TRUE: 50 day time-series, else 100 day
### Param @jitter: indicator of if to jitter covariates or not
### Param @remove_low_case_study: if TRUE: return NULL if more than 20 percent of simulation has less than 10 cases
### Returns: list of specified simulation parameters
#############################################################################################

call_specific_simulation = function(seed, ds, tp, short_run, jitter=TRUE, remove_low_case_study=TRUE) {
  set.seed(seed)
  scenario_params = get_scenario_based_params(ds)
  pre_sim_params = get_pre_simulation_params(tp, short_run, scenario_params, jitter)
  sim_params = simulate_outbreak(pre_sim_params, scenario_params)
  if (sum(sim_params$newcase < 10) >= (sim_params$study.length / 5)) { # Throw away uninteresting cases
      warning("Low Case Study Warning: More than 20 percent of study has less than 10 cases")
      if (remove_low_case_study) return(NULL)
    }
  return(list(sim_params=sim_params, scenario_params=scenario_params))
}

get_scenario_based_params = function(scenario) {
  # This function sets the specifics for a given scenario for the simulations
  # autoregression_present: {TRUE/FALSE} determines if simulation uses autoregression term
  # serial_known: {TRUE/FALSE} determines if true serial distribution is unknown mixture or known constant
  # predictor_specified: {TRUE/FALSE} determines if predictors are specified or considered unknown (pass back white noise)
  # parametric_mispecification: {TRUE/FALSE} determines if simulation generates from same distributions as model
  # measurement_error_specification: {0,1,2,3,4} 0=simulation doesn't use measurement error, 1=simulation uses model specified,
  # 2=simulation uses time-dependent variation, 3=simulation uses non-normalized weights, 4=simulation has weights depend on Rt

  if (scenario == "0") { # Baseline Assumptions
    autoregression_present = FALSE
    serial_known = TRUE
    predictor_specified = FALSE
    parametric_misspecification = FALSE
    measurement_error_specification = 0
  }

  if (scenario == "1a") { # Baseline Assumptions + Serial Interval
    autoregression_present = FALSE
    serial_known = FALSE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 0
  }

  if (scenario == "1b") { # Baseline Assumptions + Predictor Model
    autoregression_present = TRUE
    serial_known = TRUE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 0
  }

  if (scenario == "1c") { # Baseline Assumptions + Measurement Error Model
    autoregression_present = FALSE
    serial_known = TRUE
    predictor_specified = FALSE
    parametric_misspecification = FALSE
    measurement_error_specification = 1
  }

  if (scenario == "2a") { # Baseline Assumptions + Serial Interval + Predictor + Measurement Error
    autoregression_present = TRUE
    serial_known = FALSE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 1
  }

  if (scenario == "2b") { # Proposed model but parametric misspecification
    autoregression_present = TRUE
    serial_known = FALSE
    predictor_specified = TRUE
    parametric_misspecification = TRUE
    measurement_error_specification = 1
  }

  if (scenario == "3a") { # Proposed model but measurement error time-dependent
    autoregression_present = TRUE
    serial_known = FALSE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 2
  }

  if (scenario == "3b") { # Proposed model but measurement error non-normalized
    autoregression_present = TRUE
    serial_known = FALSE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 3
  }

  if (scenario == "3c") { # Proposed model but measurement error depends on Rt
    autoregression_present = TRUE
    serial_known = TRUE
    predictor_specified = TRUE
    parametric_misspecification = FALSE
    measurement_error_specification = 4
  }
  
  return(list(autoregression_present = autoregression_present,
              serial_known = serial_known,
              predictor_specified = predictor_specified,
              parametric_misspecification = parametric_misspecification,
              measurement_error_specification = measurement_error_specification))
}

