# =============================================================================
# QUESTIONS: outcome-model + small-sample validation (covariate-aware)
# Append to questions.R. Produces the Table 1 columns directly and now:
#   - accepts optional covariates, carried through EVERY component (full model,
#     reduced/permutation model, bootstrap, CV) so the spillover test stays a
#     PARTIAL test conditional on direct + covariates;
#   - consumes the long-format object from prepare_spillover_terms() via
#     fit_spillover_screen(), so no manual per-candidate looping.
# Depends only on stats.
# =============================================================================


#' Fit the direct-plus-spillover outcome model for one module and validate it
#'
#' Fits \eqn{y = \beta_0 + \beta_{dir}\,\mathrm{direct} + \beta_{spill}\,\mathrm{spill}
#' + \gamma^\top Z + \varepsilon} with direct and spillover standardized, and
#' evaluates the partial spillover coefficient with the small-sample battery:
#' a Freedman--Lane permutation test (permuting residuals of the reduced model,
#' which contains direct AND any covariates, so the test conditions on both),
#' a case bootstrap for sign-stability and a CI, repeated k-fold CV for the added
#' out-of-sample R^2 over a direct+covariates baseline, and leave-one-out.
#'
#' @param y Numeric outcome (one value per unit/row).
#' @param direct,spill Numeric vectors, same length/order as y.
#' @param covariates Optional numeric matrix/data.frame of additional predictors
#'   (rows = units, same order as y). Included in the full model, the reduced
#'   (permutation/CV baseline) model, the bootstrap, and CV. Default NULL.
#' @param family GLM family (default gaussian()); permutation/CV assume identity-link LM.
#' @param n_perm,n_boot Permutation and bootstrap counts (defaults 5000).
#' @param cv_k,cv_reps CV folds and repeats (defaults 5, 50).
#' @param conf Bootstrap CI level (default 0.95).
#' @param seed Optional integer seed.
#' @return One-row data.frame: n, beta_spill, ci_lo, ci_hi, sign_stab, p_perm,
#'   cv_dR2, p_direct, loo_max_p.
#' @seealso \code{\link{fit_spillover_screen}}, \code{\link{prepare_spillover_terms}}
#' @export
fit_spillover_outcome <- function(y, direct, spill, covariates = NULL,
                                  family = stats::gaussian(),
                                  n_perm = 5000, n_boot = 5000,
                                  cv_k = 5, cv_reps = 50,
                                  conf = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z <- if (is.null(covariates)) NULL else as.matrix(covariates)
  ok <- if (is.null(Z)) stats::complete.cases(y, direct, spill) else stats::complete.cases(y, direct, spill, Z)
  y <- as.numeric(y)[ok]; direct <- as.numeric(direct)[ok]; spill <- as.numeric(spill)[ok]
  if (!is.null(Z)) Z <- Z[ok, , drop = FALSE]
  n <- length(y); if (n < 5) stop("Need at least 5 complete units.")
  is_lm <- family$family == "gaussian" && family$link == "identity"

  zs <- function(v) as.numeric(scale(v))
  d <- zs(direct); s <- zs(spill)
  dat <- data.frame(y = y, d = d, s = s)
  if (!is.null(Z)) { colnames(Z) <- paste0("Z", seq_len(ncol(Z))); dat <- cbind(dat, Z) }
  cov_terms <- if (is.null(Z)) "" else paste0(" + ", paste(colnames(Z), collapse = " + "))

  f_full <- stats::as.formula(paste0("y ~ d + s", cov_terms))
  f_red  <- stats::as.formula(paste0("y ~ d",     cov_terms))   # reduced = direct (+ covariates)
  beta_of <- function(dd) stats::coef(stats::glm(f_full, family = family, data = dd))[["s"]]

  full     <- stats::glm(f_full, family = family, data = dat)
  beta_spill <- stats::coef(full)[["s"]]
  p_direct <- summary(full)$coefficients["d", 4]

  # --- Freedman-Lane: permute residuals of the reduced model (direct + covariates) ---
  reduced <- stats::glm(f_red, family = family, data = dat)
  fit_r <- stats::predict(reduced, type = "response"); res_r <- y - fit_r
  perm <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    dd <- dat; dd$y <- fit_r + res_r[sample.int(n)]
    perm[b] <- beta_of(dd)
  }
  p_perm <- (1 + sum(abs(perm) >= abs(beta_spill))) / (n_perm + 1)

  # --- Case bootstrap: sign-stability + CI ---
  boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    ix <- sample.int(n, replace = TRUE)
    boot[b] <- tryCatch(beta_of(dat[ix, , drop = FALSE]), error = function(e) NA_real_)
  }
  boot <- boot[is.finite(boot)]; a <- (1 - conf) / 2
  ci <- stats::quantile(boot, c(a, 1 - a), names = FALSE)
  sign_stab <- mean(sign(boot) == sign(beta_spill))

  # --- Repeated k-fold CV: added out-of-sample R^2 over the direct+covariates baseline ---
  cv_dR2 <- NA_real_
  if (is_lm) {
    gains <- numeric(cv_reps)
    for (r in seq_len(cv_reps)) {
      fold <- sample(rep_len(seq_len(cv_k), n)); pf <- pr <- rep(NA_real_, n)
      for (f in seq_len(cv_k)) {
        tr <- fold != f; te <- !tr
        if (sum(te) == 0 || sum(tr) < (3 + (if (is.null(Z)) 0 else ncol(Z)))) next
        mf <- stats::lm(f_full, data = dat[tr, , drop = FALSE])
        mr <- stats::lm(f_red,  data = dat[tr, , drop = FALSE])
        pf[te] <- stats::predict(mf, dat[te, , drop = FALSE])
        pr[te] <- stats::predict(mr, dat[te, , drop = FALSE])
      }
      good <- is.finite(pf) & is.finite(pr); sst <- sum((y[good] - mean(y[good]))^2)
      gains[r] <- (1 - sum((y[good] - pf[good])^2) / sst) - (1 - sum((y[good] - pr[good])^2) / sst)
    }
    cv_dR2 <- mean(gains, na.rm = TRUE)
  }

  # --- Leave-one-out influence ---
  loo <- numeric(n)
  for (i in seq_len(n)) loo[i] <- summary(stats::glm(f_full, family = family,
                                    data = dat[-i, , drop = FALSE]))$coefficients["s", 4]

  data.frame(n = n, beta_spill = beta_spill, ci_lo = ci[1], ci_hi = ci[2],
             sign_stab = sign_stab, p_perm = p_perm, cv_dR2 = cv_dR2,
             p_direct = p_direct, loo_max_p = max(loo), row.names = NULL)
}


#' Screen candidate module x timepoint combinations from a long-format object
#'
#' Consumes the output of prepare_spillover_terms() directly. For each candidate
#' (module, timepoint) it selects that timepoint's rows for the module, fits the
#' covariate-aware outcome model, and Benjamini--Hochberg adjusts the permutation
#' p-values across the WHOLE candidate family (the reported FDR denominator).
#'
#' @param long Long data.frame from prepare_spillover_terms() (needs columns
#'   module, timepoint, y, direct, spill, and any covariate columns).
#' @param candidates data.frame(module, timepoint) from select_candidates().
#' @param covariate_cols Optional character vector of column names in `long` to
#'   pass as covariates to the outcome model. Default NULL (none).
#' @param family,n_perm,n_boot,cv_k,cv_reps,conf Passed to fit_spillover_outcome().
#' @param seed Optional integer seed applied once before the screen.
#' @param n_family Optional integer overriding the BH family size (default = number
#'   of candidates actually fit). Set to the pre-specified family size if some
#'   candidates are dropped for insufficient data but should still count.
#' @return data.frame sorted by p_BH: module, timepoint, n, beta_spill, ci_lo,
#'   ci_hi, sign_stab, p_perm, p_BH, cv_dR2, p_direct, loo_max_p.
#' @seealso \code{\link{prepare_spillover_terms}}, \code{\link{select_candidates}}
#' @export
fit_spillover_screen <- function(long, candidates, covariate_cols = NULL,
                                 family = stats::gaussian(),
                                 n_perm = 5000, n_boot = 5000, cv_k = 5, cv_reps = 50,
                                 conf = 0.95, seed = NULL, n_family = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rows <- list()
  for (r in seq_len(nrow(candidates))) {
    m <- as.character(candidates$module[r]); t <- as.character(candidates$timepoint[r])
    sub <- long[long$module == m & long$timepoint == t, , drop = FALSE]
    sub <- sub[stats::complete.cases(sub[, c("y", "direct", "spill")]), , drop = FALSE]
    if (nrow(sub) < 8) next
    Z <- if (is.null(covariate_cols)) NULL else sub[, covariate_cols, drop = FALSE]
    fit <- tryCatch(
      fit_spillover_outcome(sub$y, sub$direct, sub$spill, covariates = Z, family = family,
                            n_perm = n_perm, n_boot = n_boot, cv_k = cv_k, cv_reps = cv_reps, conf = conf),
      error = function(e) NULL)
    if (is.null(fit)) next
    rows[[length(rows) + 1]] <- data.frame(module = m, timepoint = t, fit,
                                            stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out))
    stop("No candidates could be fit (check sample counts and long-format column names).")
  out$p_perm <- as.numeric(out$p_perm)
  fam <- if (is.null(n_family)) nrow(out) else n_family
  out$p_BH <- stats::p.adjust(out$p_perm, method = "BH", n = fam)
  out <- out[order(as.numeric(out$p_BH), as.numeric(out$p_perm)), , drop = FALSE]
  rownames(out) <- NULL
  out[, c("module","timepoint","n","beta_spill","ci_lo","ci_hi","sign_stab",
          "p_perm","p_BH","cv_dR2","p_direct","loo_max_p")]
}
