# Monte Carlo for rank-one covariance transformation in FGLS staggered DiD / ETWFE
#
# How to run:
# 1) source("rank1_fgls_did_mc.R")
# 2) quick <- run_mc(n_rep = 100, seed = 123)
# 3) quick$summary_tables$scenario_overview
# 4) full  <- run_mc(n_rep = 1000, seed = 123)
#
# Dependencies: base R + Matrix

suppressPackageStartupMessages(library(Matrix))

safe_symmetrize <- function(M) {
  0.5 * (M + t(M))
}

safe_eigen <- function(M) {
  eigen(safe_symmetrize(M), symmetric = TRUE, only.values = TRUE)$values
}

safe_condition_number <- function(M, tol = 1e-12) {
  ev <- safe_eigen(M)
  ev_abs <- abs(ev)
  if (length(ev_abs) == 0L) return(Inf)
  ev_max <- max(ev_abs)
  ev_min <- min(ev_abs)
  if (!is.finite(ev_max) || !is.finite(ev_min) || ev_min < tol) return(Inf)
  ev_max / ev_min
}

safe_inverse <- function(M, require_spd = FALSE) {
  M <- safe_symmetrize(M)
  if (require_spd) {
    out <- tryCatch({
      R <- chol(M)
      chol2inv(R)
    }, error = function(e) NULL)
    return(out)
  }
  out <- tryCatch({
    R <- chol(M)
    chol2inv(R)
  }, error = function(e) {
    tryCatch(solve(M), error = function(e2) NULL)
  })
  out
}

build_treatment_index <- function(TT = 10L, treat_start_min = 6L, treat_start_max = 10L) {
  cohorts <- seq.int(treat_start_min, treat_start_max)
  idx <- vector("list", length = 0L)
  nm <- character(0L)
  ccount <- 1L
  for (g in cohorts) {
    for (t in g:TT) {
      idx[[ccount]] <- c(g, t)
      nm[ccount] <- paste0("ATT_g", g, "_t", t)
      ccount <- ccount + 1L
    }
  }
  mat <- do.call(rbind, idx)
  colnames(mat) <- c("g", "t")
  list(map = mat, names = nm)
}

generate_data <- function(
    N = 1000L,
    TT = 10L,
    scenario = c("A", "B", "C"),
    seed = NULL,
    treat_start_min = 6L,
    treat_start_max = 10L,
    p_never = 0.20) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  pars <- switch(
    scenario,
    A = list(rho1 = 0.00, rho2 = 0.00, rho_het = 0.00, sigma_c = 1.0, sigma_e = 1.0),
    B = list(rho1 = 0.55, rho2 = 0.25, rho_het = 0.10, sigma_c = 1.5, sigma_e = 1.0),
    C = list(rho1 = 0.92, rho2 = 0.05, rho_het = 0.18, sigma_c = 6.0, sigma_e = 0.4)
  )

  tvec <- seq_len(TT)
  treat_grid <- seq.int(treat_start_min, treat_start_max)

  c_i <- rnorm(N, mean = 0, sd = pars$sigma_c)

  latent <- 0.8 * scale(c_i)[, 1] + rnorm(N, sd = 0.6)
  q <- quantile(latent, probs = seq(0, 1, length.out = length(treat_grid) + 2L))
  cohort <- rep(Inf, N)
  never_cut <- quantile(latent, probs = p_never)
  never_idx <- latent <= never_cut
  cohort[never_idx] <- Inf
  tr_idx <- which(!never_idx)
  if (length(tr_idx) > 0L) {
    bins <- cut(latent[tr_idx], breaks = q, include.lowest = TRUE, labels = FALSE)
    bins[is.na(bins)] <- 1L
    bins <- pmin(pmax(bins - 1L, 1L), length(treat_grid))
    cohort[tr_idx] <- treat_grid[bins]
  }

  cohort_shift <- ifelse(is.finite(cohort), (cohort - mean(treat_grid)) / sd(treat_grid), 0)
  X1 <- as.numeric(0.7 * scale(c_i)[, 1] + 0.5 * cohort_shift + rnorm(N))
  p_x2 <- plogis(-0.1 + 0.7 * scale(c_i)[, 1] + 0.6 * cohort_shift)
  X2 <- rbinom(N, size = 1, prob = p_x2)

  U <- matrix(0, nrow = N, ncol = TT)
  e <- matrix(rnorm(N * TT, sd = pars$sigma_e), nrow = N, ncol = TT)
  sigma_star <- sqrt(1 + pars$rho1^2 + pars$rho2^2)
  for (tt in seq_len(TT)) {
    u_l1 <- if (tt > 1L) U[, tt - 1L] else 0
    u_l2 <- if (tt > 2L) U[, tt - 2L] else 0
    scale_t <- (1 + pars$rho_het * tt) / sigma_star
    U[, tt] <- (pars$rho1 * u_l1 + pars$rho2 * u_l2 + e[, tt]) * scale_t
  }
  V <- U + c_i

  alpha_t <- 0.3 * tvec + 0.08 * (tvec - 5)^2 / 5
  lambda_g <- ifelse(is.finite(cohort), 0.15 * (cohort - 7), 0.0)

  Y0 <- matrix(0, nrow = N, ncol = TT)
  tau <- matrix(0, nrow = N, ncol = TT)
  D <- matrix(0, nrow = N, ncol = TT)

  for (tt in seq_len(TT)) {
    Y0[, tt] <- alpha_t[tt] +
      lambda_g +
      0.9 * X1 + 0.5 * X2 +
      0.08 * tt * X1 +
      0.12 * cohort_shift * X2 +
      V[, tt]

    is_treated_t <- is.finite(cohort) & (tt >= cohort)
    D[, tt] <- as.numeric(is_treated_t)
    rel_t <- pmax(tt - cohort, 0)
    rel_t[!is.finite(rel_t)] <- 0
    tau[, tt] <- is_treated_t * (
      0.6 + 0.12 * rel_t + 0.08 * (10 - cohort) + 0.10 * X1 + 0.12 * X2
    )
  }

  Y <- Y0 + tau

  treat_index <- build_treatment_index(TT, treat_start_min, treat_start_max)
  true_att <- setNames(rep(NA_real_, nrow(treat_index$map)), treat_index$names)
  for (k in seq_len(nrow(treat_index$map))) {
    g <- treat_index$map[k, "g"]
    tt <- treat_index$map[k, "t"]
    idx <- which(cohort == g)
    if (length(idx) > 0L) {
      true_att[k] <- mean(tau[idx, tt])
    }
  }

  list(
    N = N, TT = TT,
    id = rep(seq_len(N), each = TT),
    time = rep(seq_len(TT), times = N),
    cohort = cohort,
    X1 = X1,
    X2 = X2,
    Y = as.vector(t(Y)),
    Y_mat = Y,
    treat_index = treat_index,
    true_att = true_att,
    D_mat = D,
    tau_mat = tau,
    scenario = scenario,
    pars = pars
  )
}

build_design_matrix <- function(dat, att_type = c("cohort_time")) {
  att_type <- match.arg(att_type)
  if (att_type != "cohort_time") {
    stop("Only att_type = 'cohort_time' is supported in this specification.")
  }

  N <- dat$N
  TT <- dat$TT
  n_obs <- N * TT
  y <- dat$Y
  time <- dat$time
  cohort_rep <- rep(dat$cohort, each = TT)

  X1_rep <- rep(dat$X1, each = TT)
  X2_rep <- rep(dat$X2, each = TT)

  # Time FE: t2,...,tT (drop t1)
  time_fe <- model.matrix(~ factor(time) - 1)
  time_fe <- time_fe[, -1, drop = FALSE]
  colnames(time_fe) <- paste0("t", 2:TT)

  # Cohort FE over treated cohorts only; drop baseline (g=6)
  treated_cohorts <- sort(unique(dat$treat_index$map[, "g"]))
  g_baseline <- min(treated_cohorts)
  cohort_levels <- setdiff(treated_cohorts, g_baseline)

  cohort_fe_unit <- sapply(cohort_levels, function(g) as.numeric(dat$cohort == g))
  if (is.null(dim(cohort_fe_unit))) {
    cohort_fe_unit <- matrix(cohort_fe_unit, ncol = 1)
  }
  colnames(cohort_fe_unit) <- paste0("g", cohort_levels)
  cohort_fe <- cohort_fe_unit[rep(seq_len(N), each = TT), , drop = FALSE]

  # Nuisance matrix W
  intercept <- matrix(1, nrow = n_obs, ncol = 1, dimnames = list(NULL, "intercept"))
  x_main <- cbind(x1 = X1_rep, x2 = X2_rep)

  cohort_x1 <- cohort_fe * X1_rep
  cohort_x2 <- cohort_fe * X2_rep
  colnames(cohort_x1) <- paste0(colnames(cohort_fe), "_x1")
  colnames(cohort_x2) <- paste0(colnames(cohort_fe), "_x2")

  time_x1 <- time_fe * X1_rep
  time_x2 <- time_fe * X2_rep
  colnames(time_x1) <- paste0(colnames(time_fe), "_x1")
  colnames(time_x2) <- paste0(colnames(time_fe), "_x2")

  W <- cbind(intercept, time_fe, cohort_fe, x_main, cohort_x1, cohort_x2, time_x1, time_x2)

  # Treatment matrix Z: post-treatment cohort x time indicators + covariate interactions
  idx_map <- dat$treat_index$map
  K <- nrow(idx_map)
  Z_base <- matrix(0, nrow = n_obs, ncol = K)
  z_names <- character(K)
  for (k in seq_len(K)) {
    g <- idx_map[k, "g"]
    tt <- idx_map[k, "t"]
    Z_base[, k] <- as.numeric((cohort_rep == g) & (time == tt))
    z_names[k] <- paste0("ATT_g", g, "_t", tt)
  }
  colnames(Z_base) <- z_names

  Z_x1 <- Z_base * X1_rep
  Z_x2 <- Z_base * X2_rep
  colnames(Z_x1) <- paste0(z_names, "_x1")
  colnames(Z_x2) <- paste0(z_names, "_x2")

  Z <- cbind(Z_base, Z_x1, Z_x2)

  # Truth vectors aligned with treatment cells and treatment interactions
  true_base <- dat$true_att[z_names]
  names(true_base) <- z_names
  true_x1 <- setNames(rep(NA_real_, K), paste0(z_names, "_x1"))
  true_x2 <- setNames(rep(NA_real_, K), paste0(z_names, "_x2"))
  true_treat <- c(true_base, true_x1, true_x2)

  X <- cbind(W, Z)

  stopifnot(
    nrow(X) == length(y),
    all(is.finite(X)),
    all(is.finite(y)),
    qr(W)$rank == ncol(W)
  )

  list(
    y = y,
    X = X,
    W = W,
    Z = Z,
    treat_names = colnames(Z),
    true_treat = true_treat,
    true_treat_base = true_base,
    true_treat_x1 = true_x1,
    true_treat_x2 = true_x2,
    has_intercept = any(colnames(W) == "intercept")
  )
}


estimate_ols_residuals <- function(y, X) {
  fit <- lm.fit(x = X, y = y)
  list(resid = fit$residuals, coef = fit$coefficients)
}

estimate_Omega_hat <- function(resid, N, TT) {
  Vhat <- matrix(resid, nrow = N, ncol = TT, byrow = TRUE)
  crossprod(Vhat) / N
}

choose_a_transform <- function(
    Omega_hat,
    method = c("grid", "adaptive"),
    require_spd = FALSE,
    kappa_threshold = 1e8,
    grid_length = 401L,
    grid_scale = 2.0) {
  method <- match.arg(method)
  TT <- nrow(Omega_hat)
  one <- matrix(1, nrow = TT, ncol = 1)
  J <- one %*% t(one)

  kappa0 <- safe_condition_number(Omega_hat)

  eval_candidate <- function(a) {
    Oa <- Omega_hat - a * J
    eigv <- safe_eigen(Oa)
    min_ev <- min(eigv)
    inv <- safe_inverse(Oa, require_spd = require_spd)
    valid <- !is.null(inv)
    list(valid = valid, kappa = if (valid) safe_condition_number(Oa) else Inf,
         min_eig = min_ev, Oa = Oa, inv = inv)
  }

  if (method == "adaptive" && is.finite(kappa0) && kappa0 <= kappa_threshold) {
    base <- eval_candidate(0)
    return(list(a = 0, Omega_a = Omega_hat, Omega_a_inv = base$inv,
                kappa_before = kappa0, kappa_after = base$kappa,
                min_eig_before = min(safe_eigen(Omega_hat)),
                min_eig_after = base$min_eig,
                search_method = method))
  }

  scale_a <- grid_scale * mean(diag(Omega_hat)) / TT
  if (!is.finite(scale_a) || scale_a <= 0) scale_a <- 1
  avec <- seq(-scale_a, scale_a, length.out = grid_length)

  best <- list(kappa = Inf, a = NA_real_, inv = NULL, Oa = NULL, min_eig = NA_real_)
  for (a in avec) {
    cand <- eval_candidate(a)
    if (cand$valid && is.finite(cand$kappa) && cand$kappa < best$kappa) {
      best <- list(kappa = cand$kappa, a = a, inv = cand$inv, Oa = cand$Oa, min_eig = cand$min_eig)
    }
  }

  if (is.null(best$inv)) {
    base <- eval_candidate(0)
    return(list(a = 0, Omega_a = Omega_hat, Omega_a_inv = base$inv,
                kappa_before = kappa0, kappa_after = Inf,
                min_eig_before = min(safe_eigen(Omega_hat)),
                min_eig_after = NA_real_, search_method = paste0(method, "_fallback")))
  }

  list(
    a = best$a,
    Omega_a = best$Oa,
    Omega_a_inv = best$inv,
    kappa_before = kappa0,
    kappa_after = best$kappa,
    min_eig_before = min(safe_eigen(Omega_hat)),
    min_eig_after = best$min_eig,
    search_method = method
  )
}

fgls_estimate <- function(y, X, N, TT, Omega_inv, p_treat, require_spd = FALSE) {
  if (is.null(Omega_inv)) {
    return(list(success = FALSE, coef_treat = rep(NA_real_, p_treat), coef_all = NULL))
  }

  p <- ncol(X)
  XtWX <- matrix(0, p, p)
  XtWy <- numeric(p)

  ok <- TRUE
  for (i in seq_len(N)) {
    idx <- ((i - 1L) * TT + 1L):(i * TT)
    Xi <- X[idx, , drop = FALSE]
    yi <- y[idx]
    XiW <- t(Xi) %*% Omega_inv
    XtWX <- XtWX + XiW %*% Xi
    XtWy <- XtWy + as.vector(XiW %*% yi)
  }

  b <- tryCatch({
    solve(XtWX, XtWy)
  }, error = function(e) {
    ok <<- FALSE
    rep(NA_real_, p)
  })

  list(success = ok, coef_treat = tail(b, p_treat), coef_all = b)
}

run_one_rep <- function(
    scenario = c("A", "B", "C"),
    N = 1000L,
    TT = 10L,
    seed = NULL,
    att_type = c("cohort_time"),
    a_method = c("grid", "adaptive"),
    tol_invariance = 1e-8,
    kappa_threshold = 1e8) {

  scenario <- match.arg(scenario)
  att_type <- match.arg(att_type)
  a_method <- match.arg(a_method)

  dat <- generate_data(N = N, TT = TT, scenario = scenario, seed = seed)
  dm <- build_design_matrix(dat, att_type = att_type)

  if (!dm$has_intercept) stop("Invariance check requires intercept in nuisance regressors.")

  ols <- estimate_ols_residuals(dm$y, dm$X)
  Omega_hat <- estimate_Omega_hat(ols$resid, N = N, TT = TT)

  Omega_inv_std <- safe_inverse(Omega_hat, require_spd = FALSE)
  est_std <- fgls_estimate(dm$y, dm$X, N, TT, Omega_inv_std, p_treat = ncol(dm$Z))

  a_obj <- choose_a_transform(
    Omega_hat,
    method = a_method,
    require_spd = FALSE,
    kappa_threshold = kappa_threshold
  )
  est_tr <- fgls_estimate(dm$y, dm$X, N, TT, a_obj$Omega_a_inv, p_treat = ncol(dm$Z))

  max_diff <- max(abs(est_std$coef_treat - est_tr$coef_treat), na.rm = TRUE)
  if (!is.finite(max_diff)) max_diff <- NA_real_
  invariant_violation <- is.finite(max_diff) && (max_diff > tol_invariance)

  debug_obj <- NULL
  if (invariant_violation) {
    debug_obj <- list(
      Omega_hat = Omega_hat,
      Omega_a = a_obj$Omega_a,
      coef_std = est_std$coef_treat,
      coef_trans = est_tr$coef_treat
    )
  }

  list(
    scenario = scenario,
    treat_names = dm$treat_names,
    true_treat = dm$true_treat,
    beta_std = est_std$coef_treat,
    beta_tr = est_tr$coef_treat,
    max_abs_diff = max_diff,
    invariant_violation = invariant_violation,
    kappa_Omega = safe_condition_number(Omega_hat),
    kappa_Omega_a = a_obj$kappa_after,
    min_eig_Omega = min(safe_eigen(Omega_hat)),
    min_eig_Omega_a = a_obj$min_eig_after,
    a_chosen = a_obj$a,
    std_fail = !est_std$success,
    tr_fail = !est_tr$success,
    debug = debug_obj
  )
}

summarize_results <- function(res_list) {
  scenarios <- unique(vapply(res_list, function(x) x$scenario, character(1)))

  summ_coef <- list()
  summ_scen <- list()

  for (sc in scenarios) {
    sub <- res_list[vapply(res_list, function(x) x$scenario == sc, logical(1))]
    treat_names <- sub[[1]]$treat_names
    true_vec <- sub[[1]]$true_treat

    B_std <- do.call(rbind, lapply(sub, function(x) x$beta_std))
    B_tr <- do.call(rbind, lapply(sub, function(x) x$beta_tr))

    coef_tab <- data.frame(
      scenario = sc,
      coef = treat_names,
      true = as.numeric(true_vec),
      mean_std = colMeans(B_std, na.rm = TRUE),
      mean_tr = colMeans(B_tr, na.rm = TRUE),
      bias_std = colMeans(B_std, na.rm = TRUE) - as.numeric(true_vec),
      bias_tr = colMeans(B_tr, na.rm = TRUE) - as.numeric(true_vec),
      sd_std = apply(B_std, 2, sd, na.rm = TRUE),
      sd_tr = apply(B_tr, 2, sd, na.rm = TRUE),
      rmse_std = sqrt(colMeans((t(t(B_std) - as.numeric(true_vec)))^2, na.rm = TRUE)),
      rmse_tr = sqrt(colMeans((t(t(B_tr) - as.numeric(true_vec)))^2, na.rm = TRUE))
    )
    summ_coef[[sc]] <- coef_tab

    k0 <- vapply(sub, function(x) x$kappa_Omega, numeric(1))
    ka <- vapply(sub, function(x) x$kappa_Omega_a, numeric(1))
    md <- vapply(sub, function(x) x$max_abs_diff, numeric(1))

    scen_tab <- data.frame(
      scenario = sc,
      n_rep = length(sub),
      fail_std_rate = mean(vapply(sub, function(x) x$std_fail, logical(1))),
      fail_tr_rate = mean(vapply(sub, function(x) x$tr_fail, logical(1))),
      invariant_violation_rate = mean(vapply(sub, function(x) x$invariant_violation, logical(1))),
      mean_kappa_Omega = mean(k0, na.rm = TRUE),
      median_kappa_Omega = median(k0, na.rm = TRUE),
      mean_kappa_Omega_a = mean(ka, na.rm = TRUE),
      median_kappa_Omega_a = median(ka, na.rm = TRUE),
      mean_max_abs_diff = mean(md, na.rm = TRUE),
      p95_max_abs_diff = as.numeric(quantile(md, 0.95, na.rm = TRUE))
    )
    summ_scen[[sc]] <- scen_tab
  }

  list(
    by_coefficient = do.call(rbind, summ_coef),
    scenario_overview = do.call(rbind, summ_scen)
  )
}

run_mc <- function(
    n_rep = 1000L,
    scenarios = c("A", "B", "C"),
    N = 1000L,
    TT = 10L,
    seed = 123,
    att_type = c("cohort_time"),
    a_method = c("grid", "adaptive"),
    progress_every = 50L,
    tol_invariance = 1e-8,
    kappa_threshold = 1e8) {

  att_type <- match.arg(att_type)
  a_method <- match.arg(a_method)
  set.seed(seed)

  out <- list()
  cc <- 1L
  for (sc in scenarios) {
    for (r in seq_len(n_rep)) {
      rep_seed <- sample.int(.Machine$integer.max, 1)
      out[[cc]] <- run_one_rep(
        scenario = sc,
        N = N,
        TT = TT,
        seed = rep_seed,
        att_type = att_type,
        a_method = a_method,
        tol_invariance = tol_invariance,
        kappa_threshold = kappa_threshold
      )
      if (r %% progress_every == 0L) {
        message(sprintf("Scenario %s: completed %d / %d", sc, r, n_rep))
      }
      cc <- cc + 1L
    }
  }

  summary_tables <- summarize_results(out)
  list(raw = out, summary_tables = summary_tables)
}

invariance_grid_test <- function(
    scenario = "C",
    N = 1000L,
    TT = 10L,
    seed = 999,
    att_type = "cohort_time",
    grid_length = 201L,
    tol = 1e-8) {

  dat <- generate_data(N = N, TT = TT, scenario = scenario, seed = seed)
  dm <- build_design_matrix(dat, att_type = att_type)
  ols <- estimate_ols_residuals(dm$y, dm$X)
  Omega_hat <- estimate_Omega_hat(ols$resid, N = N, TT = TT)

  one <- matrix(1, TT, 1)
  J <- one %*% t(one)
  scale_a <- 2 * mean(diag(Omega_hat)) / TT
  avec <- seq(-scale_a, scale_a, length.out = grid_length)

  base_inv <- safe_inverse(Omega_hat)
  base_est <- fgls_estimate(dm$y, dm$X, N, TT, base_inv, ncol(dm$Z))$coef_treat

  rec <- data.frame(a = avec, valid = FALSE, kappa = NA_real_, max_abs_diff = NA_real_)
  for (k in seq_along(avec)) {
    Oa <- Omega_hat - avec[k] * J
    inv <- safe_inverse(Oa)
    if (!is.null(inv)) {
      b <- fgls_estimate(dm$y, dm$X, N, TT, inv, ncol(dm$Z))$coef_treat
      rec$valid[k] <- TRUE
      rec$kappa[k] <- safe_condition_number(Oa)
      rec$max_abs_diff[k] <- max(abs(b - base_est), na.rm = TRUE)
    }
  }

  valid_diffs <- rec$max_abs_diff[rec$valid]
  cat("Invariance grid test:\n")
  cat(" valid a count:", sum(rec$valid), "of", length(avec), "\n")
  cat(" max coefficient difference across valid a:", max(valid_diffs, na.rm = TRUE), "\n")
  cat(" all within tolerance?", all(valid_diffs <= tol, na.rm = TRUE), "\n")

  print(head(rec[rec$valid, ], 10))
  print(tail(rec[rec$valid, ], 10))

  invisible(rec)
}

# Example quick run (100 replications):
# quick_res <- run_mc(n_rep = 100, seed = 123, a_method = "adaptive")
# quick_res$summary_tables$scenario_overview
#
# Example full run (1000 replications):
# full_res <- run_mc(n_rep = 1000, seed = 123, a_method = "adaptive")
# full_res$summary_tables$scenario_overview
#
# Standalone invariance check on one dataset:
# inv_grid <- invariance_grid_test(scenario = "C", seed = 42)
