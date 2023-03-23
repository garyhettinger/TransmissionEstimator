# Code to estimate measurement error periods

#############################################################################################
## get_me_model_info(): gets data for estimated reporting period patterns
### Param @study_data: list with relevant parameters for study
### Returns: list of data for measurement error model
#############################################################################################

get_me_model_info = function(study_data) {
  orig.cases = study_data$newcase
  outlier.data = remove_outliers(orig.cases)
  outlier.cases = outlier.data$outlier.free.cases
  ma.cases = get_rolling_mean(outlier.cases, 7)
  detrended.cases = (outlier.cases - ma.cases) / ma.cases
  map.data = data.frame(weekday=study_data$weekday, newcase=detrended.cases)
  map.info = get_theta_mapping_me(map.data)
  return(list(map.info=map.info, outlier.data=outlier.data, detrended.data=list(ma.cases=ma.cases, detrended.cases=detrended.cases)))
}

remove_outliers = function(cases) {
  num.outliers = length(forecast::tsoutliers(cases, lambda="auto")$index)
  outlier.free.cases = as.integer(round(forecast::tsclean(cases, lambda="auto")))
  return(list(num.outliers=num.outliers, outlier.free.cases=outlier.free.cases))
}

remove_loess = function(cases) {
  t = 1:length(cases)
  county.loess = loess(cases ~ t, span=0.2)
  detrended.cases = cases - county.loess$fitted
  detrended.cases = detrended.cases - min(detrended.cases) + 1
  return(list(loess.cases=county.loess$fitted, detrended.cases=detrended.cases))
}

adjust_cases_by_weights = function(cases, map.info, weekday.data) {
  return(sapply(1:length(weekday.data), function(x) cases[x] / 
         map.info$prior[map.info$map[weekday.data[x]]]))
}

get_rolling_averages = function(study_data, pred_data, agg_pred=TRUE, window=7) {
  outlier.free.cases = remove_outliers(study_data$newcase)$outlier.free.cases
  rolling_data = data.frame(date=study_data$date)
  rolling_data[["newcase"]] = as.integer(round(get_rolling_mean(outlier.free.cases, window)))
  for (pred in 1:ncol(pred_data)) {
    rolling_data[[paste0("X", pred)]] = if (agg_pred) get_rolling_mean(pred_data[,pred], window) else pred_data[,pred]
  }
  return(rolling_data)
}

get_rolling_mean = function(x, n=7, align="center") {
  roll_means = data.table::frollmean(x, n:1, algo="exact", align=align)
  final_mean = roll_means[[1]]
  for (i in 2:n) {
    final_mean[is.na(final_mean)] = roll_means[[i]][is.na(final_mean)]
    if (sum(is.na(final_mean)) == 0) break
  }
  return(final_mean)
}

get_prospective_me_mapping = function(day_of_week_data, nweeks, week.means, nclusters) {
    # Gets the information on a theta mapping based on number of clusters

    if (nclusters == 1) {
      mapping = rep(1,7)
      form = as.formula("X3 ~ 0") #as.formula("X3 ~ offset(X1) + 0")
    }
    else {
      day_of_week_means = rowMeans(day_of_week_data)
      mapping = if (nclusters == 7) 1:7 else kmeans(day_of_week_means, nclusters, iter.max=50, nstart=50)$cluster
      form = as.formula("X3 ~ as.factor(X2) + 0") # as.formula("X3 ~ offset(X1) + as.factor(X2):X1 + 0")
    }
    lin_reg_data = matrix(nrow=nweeks*7, ncol=3)
    for (i in 1:7) {
      for (j in 1:nweeks) {
        lin_reg_data[(i-1)*nweeks + j,] = c(week.means[j], mapping[i], day_of_week_data[i,j]) 
      }
    }
    lin_reg_data = data.frame(lin_reg_data)
    m1 = glm(form,data=lin_reg_data,family=gaussian(link="identity"))
    return(list(m1=m1, mapping=mapping))
}

get_theta_mapping_me = function(map_data) {
  # Gets the best theta mapping and prior based on AIC

  nweeks = as.integer(length(map_data$newcase)/7)
  day_of_week_data = matrix(nrow=7, ncol=nweeks)
  for (i in 1:7) {
    day_of_week_data[i,] = map_data[map_data$weekday == i, c("newcase")][1:nweeks]
  }
  week.means = colMeans(day_of_week_data)
  best.vals = list(aic=Inf, mapping=rep(1,7), prior=c(1))
  for (nclusters in 1:7) {
    pros_mapping_info = get_prospective_me_mapping(day_of_week_data, nweeks, 
                                                   week.means, nclusters)
    aic = AIC(pros_mapping_info$m1)
    # print(paste0("K: ", nclusters, " and aic: ", aic))
    # print(pros_mapping_info$mapping)
    # print(1+tail(pros_mapping_info$m1$coefficients, nclusters))
    if (aic < best.vals$aic) {
      best.vals$aic = aic
      best.vals$mapping = pros_mapping_info$mapping
      best.vals$prior = 1 + (if (length(pros_mapping_info$m1$coefficients) == 0) c(0) 
        else tail(pros_mapping_info$m1$coefficients, nclusters))
    }
  }
  return(best.vals)
}

