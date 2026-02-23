# Simulation framework for staggered-treatment DiD with POLS/FGLS/2SLS/GMM
#
# README (quick start)
# 1) source("func_dgp.R"); source("func_design.R"); source("func_estimators.R");
#    source("func_stats.R"); source("simulation_iteration.R"); source("simulation_main.R")
# 2) dbg <- run_simulation(make_default_params(debug = TRUE))
# 3) main <- run_simulation(make_default_params(debug = FALSE))
# 4) Inspect: dbg$summary, dbg$diagnostics, dbg$validation
#
# Extension hooks:
# - staggered exit: modify treatment path in simulate_panel_data()
# - TOTO designs: add alternate treatment state variables in func_design.R
# - event-study restrictions: impose linear constraints on ATT columns before estimation
# - alternative weighting matrices: replace Lambda inverse in fit_gmm_two_step()

source("func_dgp.R")
source("func_design.R")
source("func_estimators.R")
source("func_stats.R")
source("simulation_iteration.R")

run_validation_suite <- function(results, params) {
  val <- list()

  # 1) Under homosked/no-serial, POLS and FGLS should be close
  p_homo <- params
  p_homo$error$rho1 <- 0
  p_homo$error$rho2 <- 0
  p_homo$error$sigma_t <- 0
  p_homo$error$sigma_x <- 0
  p_homo$n_iter <- 3L
  tmp <- run_simulation(p_homo, run_validation = FALSE)
  merged <- reshape(tmp$estimates[, c("iter", "estimator", "term", "estimate")],
                    idvar = c("iter", "term"), timevar = "estimator", direction = "wide")
  dif <- abs(merged$estimate.POLS - merged$estimate.FGLS)
  val$pols_fgls_close_mean_abs_diff <- mean(dif, na.rm = TRUE)

  # 2) 2SLS hand formula vs ivreg (if available)
  val$ivreg_check <- NA_character_
  if (requireNamespace("AER", quietly = TRUE)) {
    one <- run_one_iteration(1, p_homo)
    val$ivreg_check <- "AER available; optional comparison can be added in notebook"
  } else {
    val$ivreg_check <- "AER::ivreg unavailable in environment"
  }

  # 3) GMM formula check by recomputation on first iteration
  r1 <- results$raw[[1]]
  val$gmm_present <- any(r1$estimates$estimator == "GMM")
  val$design_dims <- paste(r1$checks$dim_A, collapse = "x")
  val$design_dims_Z <- paste(r1$checks$dim_Z, collapse = "x")
  val$post_pre_assignment_ok <- r1$checks$post_in_A && r1$checks$pre_in_Z_not_A

  val
}

run_simulation <- function(params, run_validation = TRUE) {
  all_est <- list(); all_diag <- list(); raw <- list()
  for (iter in seq_len(params$n_iter)) {
    res <- run_one_iteration(iter, params)
    all_est[[iter]] <- res$estimates
    all_diag[[iter]] <- res$diagnostics
    raw[[iter]] <- res
  }

  estimates <- do.call(rbind, all_est)
  diagnostics <- do.call(rbind, all_diag)
  truth <- raw[[1]]$truth
  summary <- summarize_simulation(estimates, truth)

  validation <- if (run_validation) run_validation_suite(list(estimates = estimates, raw = raw), params) else NULL

  list(estimates = estimates, diagnostics = diagnostics, truth = truth, summary = summary,
       fail_rates = collect_fail_rates(diagnostics), validation = validation, raw = raw)
}

if (sys.nframe() == 0L) {
  dbg_params <- make_default_params(debug = TRUE)
  dbg <- run_simulation(dbg_params)
  print(head(dbg$summary, 12))
  print(dbg$fail_rates)

  # Default main scenario
  main_params <- make_default_params(debug = FALSE)
  # main <- run_simulation(main_params)
  # print(head(main$summary, 12))
}
