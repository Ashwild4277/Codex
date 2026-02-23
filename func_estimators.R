# Estimation routines: POLS, FGLS, 2SLS (pooled/system), two-step optimal GMM, J-test

safe_chol_solve <- function(M, b, ridge = 0) {
  k <- ncol(M)
  M2 <- M + diag(ridge, k)
  out <- tryCatch({
    R <- chol((M2 + t(M2)) / 2)
    backsolve(R, forwardsolve(t(R), b))
  }, error = function(e) NULL)
  out
}

safe_invert <- function(M, ridge = 0) {
  k <- ncol(M)
  M2 <- M + diag(ridge, k)
  out <- tryCatch({
    R <- chol((M2 + t(M2)) / 2)
    chol2inv(R)
  }, error = function(e) {
    tryCatch(solve(M2), error = function(e2) NULL)
  })
  out
}

fit_pols <- function(y, A) {
  fit <- lm.fit(x = A, y = y)
  e <- fit$residuals
  list(beta = fit$coefficients, resid = e, fail = anyNA(fit$coefficients), rank = fit$rank)
}

estimate_omega <- function(resid, N, T) {
  E <- matrix(resid, nrow = N, ncol = T, byrow = TRUE)
  crossprod(E) / N
}

fit_fgls <- function(y, A, N, T, ridge = 1e-8) {
  ols <- fit_pols(y, A)
  Omega <- estimate_omega(ols$resid, N, T)
  Omega_inv <- safe_invert(Omega, ridge = ridge)
  if (is.null(Omega_inv)) {
    return(list(beta = rep(NA_real_, ncol(A)), fail = TRUE, reason = "Omega invert fail", kappa = Inf, Omega = Omega))
  }

  ymat <- matrix(y, nrow = N, ncol = T, byrow = TRUE)
  Avec <- array(A, dim = c(N, T, ncol(A)))
  Sxx <- matrix(0, ncol(A), ncol(A)); Sxy <- rep(0, ncol(A))
  for (i in seq_len(N)) {
    Ai <- Avec[i, , , drop = FALSE][1, , ]
    yi <- ymat[i, ]
    Sxx <- Sxx + t(Ai) %*% Omega_inv %*% Ai
    Sxy <- Sxy + t(Ai) %*% Omega_inv %*% yi
  }
  beta <- safe_chol_solve(Sxx, Sxy, ridge = ridge)
  list(beta = if (is.null(beta)) rep(NA_real_, ncol(A)) else as.numeric(beta), fail = is.null(beta),
       reason = if (is.null(beta)) "FGLS normal eq fail" else "", kappa = kappa(Omega), Omega = Omega)
}

fit_2sls <- function(y, A, Z, ridge = 1e-8, label = "S2SLS") {
  ZtZ_inv <- safe_invert(crossprod(Z), ridge = ridge)
  if (is.null(ZtZ_inv)) {
    return(list(beta = rep(NA_real_, ncol(A)), fail = TRUE, reason = "ZtZ invert fail", label = label, rank_Z = qr(Z)$rank))
  }
  PZA <- Z %*% (ZtZ_inv %*% crossprod(Z, A))
  PZy <- Z %*% (ZtZ_inv %*% crossprod(Z, y))
  beta <- safe_chol_solve(crossprod(A, PZA), crossprod(A, PZy), ridge = ridge)
  list(beta = if (is.null(beta)) rep(NA_real_, ncol(A)) else as.numeric(beta), fail = is.null(beta),
       reason = if (is.null(beta)) "2SLS normal eq fail" else "", label = label, rank_Z = qr(Z)$rank)
}

compute_unit_moments <- function(y, A, Z, beta, N, T) {
  ymat <- matrix(y, N, T, byrow = TRUE)
  Avec <- array(A, c(N, T, ncol(A)))
  Zvec <- array(Z, c(N, T, ncol(Z)))
  g <- matrix(0, N, ncol(Z))
  u_store <- matrix(0, N, T)
  for (i in seq_len(N)) {
    Ai <- Avec[i, , , drop = FALSE][1, , ]
    Zi <- Zvec[i, , , drop = FALSE][1, , ]
    ui <- ymat[i, ] - Ai %*% beta
    g[i, ] <- as.numeric(t(Zi) %*% ui)
    u_store[i, ] <- ui
  }
  list(g = g, u = u_store)
}

fit_gmm_two_step <- function(y, A, Z, N, T, ridge = 1e-8) {
  # Step 1: identity weighting == 2SLS-style start
  start <- fit_2sls(y, A, Z, ridge = ridge, label = "GMM_step1")
  if (start$fail) {
    return(list(beta = start$beta, fail = TRUE, reason = paste("Step1 failed:", start$reason),
                Lambda = matrix(NA_real_, ncol(Z), ncol(Z)), kappa_Lambda = Inf, J = NA_real_, J_df = NA_integer_, J_p = NA_real_))
  }

  moms <- compute_unit_moments(y, A, Z, start$beta, N, T)
  Lambda <- crossprod(moms$g) / N
  W <- safe_invert(Lambda, ridge = ridge)
  if (is.null(W)) {
    return(list(beta = rep(NA_real_, ncol(A)), fail = TRUE, reason = "Lambda invert fail",
                Lambda = Lambda, kappa_Lambda = kappa(Lambda), J = NA_real_, J_df = NA_integer_, J_p = NA_real_))
  }

  AZWZA <- crossprod(A, Z %*% (W %*% crossprod(Z, A)))
  AZWZy <- crossprod(A, Z %*% (W %*% crossprod(Z, y)))
  beta2 <- safe_chol_solve(AZWZA, AZWZy, ridge = ridge)
  if (is.null(beta2)) {
    return(list(beta = rep(NA_real_, ncol(A)), fail = TRUE, reason = "Step2 normal eq fail",
                Lambda = Lambda, kappa_Lambda = kappa(Lambda), J = NA_real_, J_df = NA_integer_, J_p = NA_real_))
  }

  moms2 <- compute_unit_moments(y, A, Z, beta2, N, T)
  gbar <- colMeans(moms2$g)
  J <- N * as.numeric(t(gbar) %*% W %*% gbar)
  df <- ncol(Z) - ncol(A)
  p <- if (df > 0) 1 - pchisq(J, df = df) else NA_real_

  list(beta = as.numeric(beta2), fail = FALSE, reason = "", Lambda = Lambda,
       kappa_Lambda = kappa(Lambda), J = J, J_df = df, J_p = p,
       reject_5pct = ifelse(is.na(p), NA, p < 0.05))
}

j_test_excluded_only <- function(y, A, Z, excluded_terms, beta, N, T, ridge = 1e-8) {
  if (length(excluded_terms) == 0L) return(list(J = NA_real_, df = 0L, p = NA_real_, reject_5pct = NA))
  idx <- match(excluded_terms, colnames(Z))
  moms <- compute_unit_moments(y, A, Z[, idx, drop = FALSE], beta, N, T)
  gbar <- colMeans(moms$g)
  Lambda <- crossprod(moms$g) / N
  W <- safe_invert(Lambda, ridge = ridge)
  if (is.null(W)) return(list(J = NA_real_, df = length(idx), p = NA_real_, reject_5pct = NA))
  J <- N * as.numeric(t(gbar) %*% W %*% gbar)
  p <- 1 - pchisq(J, df = length(idx))
  list(J = J, df = length(idx), p = p, reject_5pct = p < 0.05)
}
