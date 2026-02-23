# Simulation summaries

summarize_simulation <- function(results_df, truth_df) {
  merged <- merge(results_df, truth_df, by = "term", all.x = TRUE)
  merged$bias <- merged$estimate - merged$truth

  agg <- aggregate(cbind(estimate, bias) ~ estimator + term, merged, function(v) c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE)))
  tidy <- data.frame(
    estimator = agg$estimator,
    term = agg$term,
    mean_estimate = agg$estimate[, "mean"],
    sd_estimate = agg$estimate[, "sd"],
    mean_bias = agg$bias[, "mean"],
    rmse = NA_real_
  )

  rm <- aggregate(bias ~ estimator + term, merged, function(v) sqrt(mean(v^2, na.rm = TRUE)))
  tidy$rmse <- rm$bias[match(paste(tidy$estimator, tidy$term), paste(rm$estimator, rm$term))]
  tidy
}

collect_fail_rates <- function(diag_df) {
  aggregate(cbind(fail_pols, fail_fgls, fail_s2sls, fail_gmm) ~ 1, diag_df, mean)
}
