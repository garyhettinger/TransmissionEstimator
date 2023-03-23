# Code to decide on smoothing approach
#############################################################################################
## run_proposed_model: run model and collect results for proposed approach
### Param @me_model_on_R: vector of estimated R time-series using measurement error model approach
### Param @case_smoothing_R: vector of estimated R time-series using case smoothing approach
### Returns: TRUE if should use measurement error model approach and FALSE if should use case smoothing
#############################################################################################

select_me_model = function(me_model_on_R, case_smoothing_R) {
  a1 = acf(me_model_on_R, plot=F)
  se1 = 1.96*sqrt((1/a1$n.used)*(1+2*a1$acf[2,1,1]^2))
  a2 = acf(case_smoothing_R, plot=F)
  CI = a2$acf[2,1,1]-(a1$acf[2,1,1] + 1.96*c(1)*sqrt(se1^2))
  if (CI[1] <= 0) return(T)
  return(F)
}