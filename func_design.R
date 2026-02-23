# Design matrix construction for structural regressors A and instruments Z

build_panel_design <- function(panel, params) {
  N <- params$N; T <- params$T; T0 <- params$T0
  gset <- seq.int(T0 + 1L, T)

  df <- panel
  df$cohort_label <- ifelse(is.finite(df$cohort), sprintf("d%04d", as.integer(df$cohort)), "dInf")
  df$time_label <- sprintf("f%02d", df$t)

  # Structural regressors in A
  intercept <- rep(1, nrow(df))
  time_mm <- model.matrix(~ factor(time_label) - 1, df)
  colnames(time_mm) <- sub("factor\\(time_label\\)", "", colnames(time_mm))
  if (ncol(time_mm) > 0L) time_mm <- time_mm[, -1, drop = FALSE] # baseline f01

  cohort_mm <- model.matrix(~ factor(cohort_label) - 1, df)
  colnames(cohort_mm) <- sub("factor\\(cohort_label\\)", "", colnames(cohort_mm))
  # baseline never-treated dInf when present, otherwise first treated cohort
  if ("dInf" %in% colnames(cohort_mm)) {
    cohort_mm <- cohort_mm[, setdiff(colnames(cohort_mm), "dInf"), drop = FALSE]
  } else if (ncol(cohort_mm) > 0L) {
    cohort_mm <- cohort_mm[, -1, drop = FALSE]
  }

  x_base <- as.matrix(df[, c("x1", "x2")])
  colnames(x_base) <- c("x1", "x2")

  cohort_x <- matrix(numeric(), nrow(df), 0)
  if (ncol(cohort_mm) > 0L) {
    cohort_x <- do.call(cbind, lapply(seq_len(ncol(cohort_mm)), function(j) {
      cbind(cohort_mm[, j] * df$x1, cohort_mm[, j] * df$x2)
    }))
    cx_names <- unlist(lapply(colnames(cohort_mm), function(cn) c(paste0(cn, ":x1"), paste0(cn, ":x2"))))
    colnames(cohort_x) <- cx_names
  }

  time_x <- matrix(numeric(), nrow(df), 0)
  if (ncol(time_mm) > 0L) {
    time_x <- do.call(cbind, lapply(seq_len(ncol(time_mm)), function(j) {
      cbind(time_mm[, j] * df$x1, time_mm[, j] * df$x2)
    }))
    tx_names <- unlist(lapply(colnames(time_mm), function(tn) c(paste0(tn, ":x1"), paste0(tn, ":x2"))))
    colnames(time_x) <- tx_names
  }

  post_cols <- list(); pre_cols <- list(); post_x_cols <- list(); pre_x_cols <- list()
  for (g in gset) {
    gnm <- sprintf("d%04d", g)
    for (tt in seq_len(T)) {
      tnm <- sprintf("f%02d", tt)
      v <- as.numeric(df$cohort == g & df$t == tt)
      nm <- sprintf("int_%s_%s", gnm, tnm)
      if (tt >= g) {
        post_cols[[nm]] <- v
        if (params$include_x_heter_att) {
          post_x_cols[[paste0(nm, ":x1")]] <- v * df$x1
          post_x_cols[[paste0(nm, ":x2")]] <- v * df$x2
        }
      } else {
        pre_cols[[nm]] <- v
        if (params$include_x_heter_att) {
          pre_x_cols[[paste0(nm, ":x1")]] <- v * df$x1
          pre_x_cols[[paste0(nm, ":x2")]] <- v * df$x2
        }
      }
    }
  }

  post_mm <- if (length(post_cols)) as.matrix(as.data.frame(post_cols)) else matrix(numeric(), nrow(df), 0)
  pre_mm <- if (length(pre_cols)) as.matrix(as.data.frame(pre_cols)) else matrix(numeric(), nrow(df), 0)
  post_x_mm <- if (length(post_x_cols)) as.matrix(as.data.frame(post_x_cols)) else matrix(numeric(), nrow(df), 0)
  pre_x_mm <- if (length(pre_x_cols)) as.matrix(as.data.frame(pre_x_cols)) else matrix(numeric(), nrow(df), 0)

  A <- cbind(`(Intercept)` = intercept, time_mm, cohort_mm, x_base)
  if (ncol(cohort_x) > 0L) A <- cbind(A, cohort_x)
  if (ncol(time_x) > 0L) A <- cbind(A, time_x)
  if (ncol(post_mm) > 0L) A <- cbind(A, post_mm)
  if (ncol(post_x_mm) > 0L) A <- cbind(A, post_x_mm)

  Z <- cbind(A)
  if (ncol(pre_mm) > 0L) Z <- cbind(Z, pre_mm)
  if (ncol(pre_x_mm) > 0L) Z <- cbind(Z, pre_x_mm)

  list(
    y = df$y,
    A = A,
    Z = Z,
    post_terms = colnames(post_mm),
    pre_terms = colnames(pre_mm),
    post_x_terms = colnames(post_x_mm),
    pre_x_terms = colnames(pre_x_mm),
    term_order_A = colnames(A),
    term_order_Z = colnames(Z)
  )
}

validate_design <- function(design, params) {
  issues <- character(0)
  if (!all(design$post_terms %in% colnames(design$A))) issues <- c(issues, "Post ATT terms missing from A")
  if (!all(design$pre_terms %in% colnames(design$Z))) issues <- c(issues, "Pre terms missing from Z")
  if (any(design$pre_terms %in% colnames(design$A))) issues <- c(issues, "Pre terms leaked into A")
  if (params$include_x_heter_att) {
    if (!all(design$post_x_terms %in% colnames(design$A))) issues <- c(issues, "Post X ATT terms missing from A")
    if (!all(design$pre_x_terms %in% colnames(design$Z))) issues <- c(issues, "Pre X terms missing from Z")
  }
  list(ok = length(issues) == 0L, issues = issues)
}
