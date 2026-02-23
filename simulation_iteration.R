# Single iteration wrapper

run_one_iteration <- function(iter, params) {
  dat <- simulate_panel_data(params, iter_seed = params$seed + iter)
  des <- build_panel_design(dat$panel, params)
  vcheck <- validate_design(des, params)

  pols <- fit_pols(des$y, des$A)
  fgls <- fit_fgls(des$y, des$A, N = params$N, T = params$T, ridge = params$ridge)
  s2sls <- fit_2sls(des$y, des$A, des$Z, ridge = params$ridge, label = "S2SLS")
  gmm <- fit_gmm_two_step(des$y, des$A, des$Z, N = params$N, T = params$T, ridge = params$ridge)

  jt <- j_test_excluded_only(des$y, des$A, des$Z, des$pre_terms, beta = if (!gmm$fail) gmm$beta else s2sls$beta,
                             N = params$N, T = params$T, ridge = params$ridge)

  att_terms <- des$post_terms
  make_est_df <- function(beta, estimator_name) {
    idx <- match(att_terms, colnames(des$A))
    data.frame(iter = iter, estimator = estimator_name, term = att_terms,
               estimate = beta[idx], se = NA_real_, converge_flag = 1L,
               fail_flag = as.integer(anyNA(beta[idx])), warnings = "", stringsAsFactors = FALSE)
  }

  out <- rbind(
    make_est_df(pols$beta, "POLS"),
    make_est_df(fgls$beta, "FGLS"),
    make_est_df(s2sls$beta, "S2SLS"),
    make_est_df(gmm$beta, "GMM")
  )

  truth <- data.frame(term = names(dat$true_att), truth = as.numeric(dat$true_att), stringsAsFactors = FALSE)

  diagnostics <- data.frame(
    iter = iter,
    fail_pols = pols$fail,
    fail_fgls = fgls$fail,
    fail_s2sls = s2sls$fail,
    fail_gmm = gmm$fail,
    rank_A = qr(des$A)$rank,
    rank_Z = qr(des$Z)$rank,
    kappa_A = kappa(des$A),
    kappa_Z = kappa(des$Z),
    kappa_Omega = ifelse(is.null(fgls$Omega), NA_real_, kappa(fgls$Omega)),
    kappa_Lambda = ifelse(is.null(gmm$Lambda), NA_real_, gmm$kappa_Lambda),
    design_ok = vcheck$ok,
    design_issues = paste(vcheck$issues, collapse = " | "),
    J_overid = gmm$J,
    J_overid_df = gmm$J_df,
    J_overid_p = gmm$J_p,
    J_excl = jt$J,
    J_excl_df = jt$df,
    J_excl_p = jt$p,
    J_excl_reject = jt$reject_5pct
  )

  # Validation checks
  checks <- list(
    dim_A = dim(des$A),
    dim_Z = dim(des$Z),
    post_in_A = all(des$post_terms %in% colnames(des$A)),
    pre_in_Z_not_A = all(des$pre_terms %in% colnames(des$Z)) && !any(des$pre_terms %in% colnames(des$A))
  )

  list(estimates = out, truth = truth, diagnostics = diagnostics, checks = checks)
}
