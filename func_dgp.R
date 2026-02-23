# DGP helpers for staggered-treatment DiD simulations

make_default_params <- function(debug = FALSE) {
  if (debug) {
    list(
      N = 100L,
      T = 6L,
      T0 = 3L,
      n_iter = 5L,
      seed = 123,
      x1_shape = 2,
      x1_scale = 1,
      x2_a = 0.5,
      x2_noise_sd = 1,
      x2_threshold = 1,
      cohort_intercept = -0.2,
      cohort_coef_c = 0.8,
      cohort_coef_x1 = 0.25,
      cohort_coef_x2 = 0.5,
      never_prob_floor = 0.10,
      include_x_heter_att = FALSE,
      include_pretrend = FALSE,
      pretrend_type = "none", # none|single_shock|trend
      pretrend_size = 0.25,
      error = list(
        ar_order = 1L,
        rho1 = 0.4,
        rho2 = 0.0,
        sigma_base = 1.0,
        sigma_t = 0.08,
        sigma_x = 0.15
      ),
      dgp = list(
        beta_x1 = 0.6,
        beta_x2 = 0.35,
        beta_tx1 = 0.06,
        beta_gx1 = 0.08,
        beta_gx2 = 0.06,
        c_sd = 0.9,
        att_base = 0.5,
        att_dyn = 0.12,
        att_cohort = 0.08,
        att_x1 = 0.00,
        att_x2 = 0.00
      ),
      ridge = 1e-8,
      kappa_warn = 1e10,
      debug_checks = TRUE
    )
  } else {
    list(
      N = 500L,
      T = 10L,
      T0 = 5L,
      n_iter = 200L,
      seed = 2026,
      x1_shape = 2,
      x1_scale = 1,
      x2_a = 0.5,
      x2_noise_sd = 1,
      x2_threshold = 1,
      cohort_intercept = -0.2,
      cohort_coef_c = 0.8,
      cohort_coef_x1 = 0.25,
      cohort_coef_x2 = 0.5,
      never_prob_floor = 0.15,
      include_x_heter_att = FALSE,
      include_pretrend = FALSE,
      pretrend_type = "none",
      pretrend_size = 0.25,
      error = list(
        ar_order = 1L,
        rho1 = 0.6,
        rho2 = 0.1,
        sigma_base = 1.0,
        sigma_t = 0.10,
        sigma_x = 0.20
      ),
      dgp = list(
        beta_x1 = 0.6,
        beta_x2 = 0.35,
        beta_tx1 = 0.06,
        beta_gx1 = 0.08,
        beta_gx2 = 0.06,
        c_sd = 1.0,
        att_base = 0.5,
        att_dyn = 0.12,
        att_cohort = 0.08,
        att_x1 = 0.00,
        att_x2 = 0.00
      ),
      ridge = 1e-8,
      kappa_warn = 1e10,
      debug_checks = FALSE
    )
  }
}

name_time <- function(t) sprintf("f%02d", t)
name_cohort <- function(g) sprintf("d%04d", g)

simulate_errors <- function(N, T, x1, error_params) {
  U <- matrix(0, nrow = N, ncol = T)
  innov <- matrix(rnorm(N * T), nrow = N, ncol = T)
  x1s <- as.numeric(scale(x1))
  sigma_it <- matrix(0, nrow = N, ncol = T)

  for (tt in seq_len(T)) {
    sigma_it[, tt] <- error_params$sigma_base * exp(error_params$sigma_t * (tt - 1) + error_params$sigma_x * x1s)
    e_tt <- innov[, tt] * sigma_it[, tt]
    if (error_params$ar_order == 1L) {
      lag1 <- if (tt > 1L) U[, tt - 1L] else 0
      U[, tt] <- error_params$rho1 * lag1 + e_tt
    } else {
      lag1 <- if (tt > 1L) U[, tt - 1L] else 0
      lag2 <- if (tt > 2L) U[, tt - 2L] else 0
      U[, tt] <- error_params$rho1 * lag1 + error_params$rho2 * lag2 + e_tt
    }
  }
  list(U = U, sigma_it = sigma_it)
}

assign_cohorts <- function(c_i, x1, x2, T, T0, pars) {
  treated_cohorts <- seq.int(T0 + 1L, T)
  latent <- pars$cohort_intercept + pars$cohort_coef_c * as.numeric(scale(c_i)) +
    pars$cohort_coef_x1 * as.numeric(scale(x1)) + pars$cohort_coef_x2 * x2 + rnorm(length(c_i), sd = 0.8)

  never_cut <- quantile(latent, probs = pars$never_prob_floor)
  cohort <- rep(Inf, length(c_i))
  tr_idx <- which(latent > never_cut)
  if (length(tr_idx) > 0L) {
    probs <- seq(0, 1, length.out = length(treated_cohorts) + 1L)
    brks <- quantile(latent[tr_idx], probs = probs)
    bins <- cut(latent[tr_idx], breaks = unique(brks), include.lowest = TRUE, labels = FALSE)
    bins[is.na(bins)] <- length(treated_cohorts)
    bins <- pmin(pmax(bins, 1L), length(treated_cohorts))
    cohort[tr_idx] <- treated_cohorts[bins]
  }
  cohort
}

simulate_panel_data <- function(params, iter_seed = NULL) {
  if (!is.null(iter_seed)) set.seed(iter_seed)
  N <- params$N; T <- params$T; T0 <- params$T0
  gset <- seq.int(T0 + 1L, T)

  x1 <- rgamma(N, shape = params$x1_shape, scale = params$x1_scale)
  x2 <- as.integer(params$x2_a * x1 + rnorm(N, sd = params$x2_noise_sd) > params$x2_threshold)
  c_i <- rnorm(N, sd = params$dgp$c_sd)
  cohort <- assign_cohorts(c_i, x1, x2, T, T0, params)

  err <- simulate_errors(N, T, x1, params$error)
  U <- err$U
  alpha_t <- seq_len(T) * 0.2
  x1s <- as.numeric(scale(x1))

  cohort_center <- ifelse(is.finite(cohort), cohort - mean(gset), 0)
  lambda_g <- 0.15 * cohort_center

  Y0 <- matrix(0, N, T)
  ATT <- matrix(0, N, T)
  D <- matrix(0, N, T)

  for (tt in seq_len(T)) {
    Y0[, tt] <- alpha_t[tt] + lambda_g + params$dgp$beta_x1 * x1 + params$dgp$beta_x2 * x2 +
      params$dgp$beta_tx1 * tt * x1s + params$dgp$beta_gx1 * cohort_center * x1s +
      params$dgp$beta_gx2 * cohort_center * x2 + c_i + U[, tt]

    treated <- is.finite(cohort) & (tt >= cohort)
    rel <- pmax(tt - cohort, 0)
    rel[!is.finite(rel)] <- 0
    D[, tt] <- as.numeric(treated)

    ATT[, tt] <- treated * (params$dgp$att_base + params$dgp$att_dyn * rel +
      params$dgp$att_cohort * (T - cohort) + params$dgp$att_x1 * x1s + params$dgp$att_x2 * x2)

    if (params$include_pretrend && params$pretrend_type != "none") {
      pre <- is.finite(cohort) & (tt < cohort)
      if (params$pretrend_type == "single_shock") {
        ATT[, tt] <- ATT[, tt] + pre * as.numeric(tt == (cohort - 1L)) * params$pretrend_size
      } else if (params$pretrend_type == "trend") {
        ATT[, tt] <- ATT[, tt] + pre * params$pretrend_size * (tt - params$T0)
      }
    }
  }

  Y <- Y0 + ATT
  panel <- data.frame(
    id = rep(seq_len(N), each = T),
    t = rep(seq_len(T), times = N),
    cohort = rep(cohort, each = T),
    x1 = rep(x1, each = T),
    x2 = rep(x2, each = T),
    y = as.vector(t(Y)),
    u = as.vector(t(U)),
    d = as.vector(t(D))
  )

  true_att <- list()
  for (g in gset) {
    idxg <- which(cohort == g)
    for (tt in g:T) {
      nm <- sprintf("int_d%04d_f%02d", g, tt)
      true_att[[nm]] <- if (length(idxg) == 0L) NA_real_ else mean(ATT[idxg, tt])
    }
  }

  list(panel = panel, Y = Y, U = U, cohort = cohort, true_att = unlist(true_att), gset = gset)
}
