# TransmissionEstimator

This file implements code relevant to "Estimating the Instantaneous Reproduction Number With Imperfect Data: A Method to Account for Case-Reporting Variation and Serial Interval Uncertainty" by Hettinger, Rubin, and Huang 2023 (https://arxiv.org/abs/2302.12078). These Bayesian models estimate the instantaneous reproduction number, Rt, while incorporating predictor information and adjusting for weekly case-reporting variation and uncertainty in selecting a serial interval estimate.

## Sample script

* example_analysis_run.R: Provides example code to run proposed simulations and estimators

## Model Code

* run_model_analysis.R: Contains functions relevant to fitting MCMC models.
* get_data_for_models.R: Contains functions relevant to organizing data to pass to MCMC models
* proposed_model.bug: Proposed model used when uncertainty in serial interval and reporting patterns
* proposed_model_noserial.bug: Proposed model used when uncertainty in reporting patterns but not serial interval
* proposed_model_nomeasurement.bug: Proposed model used when uncertainty in serial interval but not reporting patterns
* proposed_model_noserial_nomeasurement.bug: Proposed model used when no uncertainty in serial interval or reporting patterns

## Reporting Period Functions 

* estimate_reporting_periods.R: contains functions relevant to estimating day-of-week reporting periods when model parsimony is of interest. Otherwise, one can set a different reporting period per day of week
* select_me_model_approach.R: code to implement test to decide between selecting measurement error modeling approach or case smoothing approach

## Simulating Data

* call_specific_simulation.R: functions to simulate and collect data for a given seed, data scenario, and trend pattern as presented in associated paper
* get_specified_parameters.R: functions to define and collect parameters necessary to specify for a given simulation setting
* get_simulated_parameters.R: functions to simulate and collect random data for a given simulation setting
* Example_DS2A_TP2_Long.RData: example data from Data Scenario 2A, Trend Pattern 2, Long Run (100 days)
