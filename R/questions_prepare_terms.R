# =============================================================================
# QUESTIONS: substrate construction (data -> long-format direct/spillover terms)
# Append to questions.R alongside the outcome-model functions. Turns raw inputs
# (network, expression, sample labels, baseline definition) into the per
# animal x timepoint x module direct and spillover terms the outcome model
# consumes. Matches the manuscript construction exactly:
#   - between matrix: row-normalize the FULL graph, THEN drop within-module
#     entries, NO renormalization (row sum = fraction of influence leaving module)
#   - departure D = expr - each unit's own mean over baseline timepoints
#   - direct = sign-aligned 1st PC of standardized module expression, minus each
#     unit's baseline mean of that component
#   - spillover = mean over module genes of (W_between %*% D)
# Reuses normalize_rows(), coerce_to_dgCMatrix(), extract_network_matrix().
# =============================================================================


#' Between-module transition weights (cross-module influence only)
#'
#' Row-normalizes the full adjacency so each row is a gene's distribution of
#' influence over all neighbors, then retains only between-module entries WITHOUT
#' renormalizing. The row sum of the result is therefore the fraction of a gene's
#' influence that leaves its own module, so genes wired mostly within-module
#' contribute little cross-module exposure. This ordering matches the QUESTIONS
#' spillover definition; renormalizing after dropping would inflate peripheral
#' genes back to unit outgoing weight and is deliberately NOT done.
#'
#' @param network Adjacency/TOM matrix or WGCNA-style object (passed through extract_network_matrix()).
#' @param membership Named vector mapping gene -> module (names must match network genes).
#' @return A sparse dgCMatrix of between-module weights (rows = source gene).
#' @seealso \code{\link{prepare_spillover_terms}}
#' @export
build_between_matrix <- function(network, membership) {
  W <- coerce_to_dgCMatrix(abs(extract_network_matrix(network)))
  g <- rownames(W)
  if (is.null(g)) stop("network must have gene names on rows/columns.")
  membership <- membership[g]
  if (anyNA(membership)) stop("membership must cover every gene in the network.")
  Wn <- normalize_rows(W, verbose = FALSE)                 # normalize the FULL graph first
  sm <- Matrix::summary(Wn)
  keep <- membership[sm$i] != membership[sm$j]             # then drop within-module entries
  Matrix::sparseMatrix(i = sm$i[keep], j = sm$j[keep], x = sm$x[keep],
                       dims = dim(Wn), dimnames = dimnames(Wn))   # NOT renormalized
}


#' Per-sample departures from each unit's own baseline
#'
#' Computes `D = expr - baseline[, unit]`, where each unit's baseline is its mean
#' expression over the samples flagged as baseline (e.g. pre-infection timepoints).
#'
#' @param expr Gene x sample expression matrix (rownames = genes, colnames = sample IDs).
#' @param sample_unit Character vector (length = ncol(expr)) giving each sample's unit/animal ID.
#' @param baseline Logical vector (length = ncol(expr)) TRUE for baseline samples.
#' @return A gene x sample matrix of departures (same dim/dimnames as expr).
#' @seealso \code{\link{prepare_spillover_terms}}
#' @export
baseline_departures <- function(expr, sample_unit, baseline) {
  stopifnot(ncol(expr) == length(sample_unit), length(sample_unit) == length(baseline))
  units <- unique(sample_unit)
  bl <- vapply(units, function(u) {
    ci <- which(sample_unit == u & baseline)
    if (!length(ci)) rep(NA_real_, nrow(expr)) else Matrix::rowMeans(expr[, ci, drop = FALSE])
  }, numeric(nrow(expr)))
  colnames(bl) <- units
  as.matrix(expr) - bl[, sample_unit]
}


#' Sign-aligned module eigengene expressed as a departure
#'
#' First principal component of the module's standardized gene expression,
#' sign-aligned to correlate positively with mean module expression, then
#' expressed as a departure by subtracting each unit's mean of that component
#' over baseline samples. This is the DIRECT term in the outcome model.
#'
#' @param expr Gene x sample expression matrix.
#' @param module_genes Genes belonging to the module.
#' @param sample_unit Per-sample unit IDs (length ncol(expr)).
#' @param baseline Logical per-sample baseline flag (length ncol(expr)).
#' @return Named numeric vector (per sample) of the eigengene departure.
#' @export
module_eigengene_departure <- function(expr, module_genes, sample_unit, baseline) {
  gm <- intersect(module_genes, rownames(expr))
  X  <- t(scale(t(as.matrix(expr[gm, , drop = FALSE]))))          # standardize genes across samples
  pc <- stats::prcomp(t(X), center = FALSE, scale. = FALSE)$x[, 1]
  if (stats::cor(pc, colMeans(as.matrix(expr[gm, , drop = FALSE]))) < 0) pc <- -pc  # sign-align
  bl <- tapply(pc[baseline], sample_unit[baseline], mean)         # each unit's baseline level
  stats::setNames(pc - bl[sample_unit], colnames(expr))
}


#' Assemble long-format direct + spillover terms for all modules and timepoints
#'
#' The one-call substrate builder. Returns a long data.frame with one row per
#' sample x module, carrying the outcome, the direct and spillover terms, and the
#' sample/unit/timepoint/module keys (plus any covariates joined through), ready
#' for select_candidates() and fit_spillover_screen(). Constructs the between
#' matrix, per-sample departures, gene-level spillover (W_between %*% D), and the
#' module eigengene departures exactly as specified in the QUESTIONS methods.
#'
#' @param network Adjacency/TOM/WGCNA object.
#' @param membership Named gene -> module vector.
#' @param expr Gene x sample expression matrix.
#' @param sample_unit Per-sample unit/animal IDs (length ncol(expr)).
#' @param sample_time Per-sample timepoint labels (length ncol(expr)).
#' @param outcome Named numeric outcome, names = unit IDs (constant per unit). Apply
#'   any transform (e.g. log1p) before passing. NULL to omit (adds no y column).
#' @param baseline_times Character vector of timepoint labels that define baseline.
#' @param covariates Optional data.frame with rownames = unit IDs, joined per unit
#'   (constant across timepoints/modules). Passed through to the model untouched.
#' @param modules Optional subset of module labels to build (default: all modules).
#' @param verbose Logical.
#' @return A long data.frame: sample, unit, timepoint, module, y (if outcome given),
#'   direct, spill, and any covariate columns. One row per sample x module.
#' @seealso \code{\link{select_candidates}}, \code{\link{fit_spillover_screen}}
#' @export
prepare_spillover_terms <- function(network, membership, expr,
                                    sample_unit, sample_time, outcome = NULL,
                                    baseline_times = c("Pre1", "Pre2"),
                                    covariates = NULL, modules = NULL, verbose = TRUE) {
  g <- rownames(expr)
  membership <- membership[g]
  baseline <- sample_time %in% baseline_times

  if (verbose) message("Building between-module weights ...")
  Wbtw <- build_between_matrix(network, membership)
  Wbtw <- Wbtw[g, g]

  if (verbose) message("Computing departures and gene-level spillover ...")
  D  <- baseline_departures(expr, sample_unit, baseline)            # genes x samples
  CS <- as.matrix(Wbtw %*% D)                                       # gene-level spillover, genes x samples

  mods <- if (is.null(modules)) sort(unique(membership)) else modules
  genes_by_mod <- split(names(membership), membership)

  rows <- vector("list", length(mods))
  for (k in seq_along(mods)) {
    m  <- mods[k]; gm <- intersect(genes_by_mod[[as.character(m)]], g)
    if (length(gm) < 2) next
    direct <- module_eigengene_departure(expr, gm, sample_unit, baseline)
    spill  <- colMeans(CS[gm, , drop = FALSE])
    rows[[k]] <- data.frame(sample = colnames(expr), unit = sample_unit,
                            timepoint = sample_time, module = as.character(m),
                            direct = as.numeric(direct), spill = as.numeric(spill),
                            stringsAsFactors = FALSE)
  }
  long <- do.call(rbind, rows)

  if (!is.null(outcome)) long$y <- unname(outcome[long$unit])
  if (!is.null(covariates)) {
    cov <- covariates[long$unit, , drop = FALSE]; rownames(cov) <- NULL
    long <- cbind(long, cov)
  }
  if (verbose) message(sprintf("Prepared %d rows: %d modules x %d samples.",
                               nrow(long), length(unique(long$module)), ncol(expr)))
  long
}


#' Select candidate module x timepoint combinations, independent of outcome
#'
#' Encodes the outcome-INDEPENDENT selection rule: rank modules by the peak
#' absolute descriptive spillover across timepoints and keep the top n, then form
#' the family as those modules crossed with all non-baseline timepoints. Selection
#' never sees the outcome, which is what makes the multiple-testing family honest.
#'
#' The ranking statistic can be supplied two ways. Preferred and reproducible:
#' pass \code{module_spillover}, the modules x timepoints descriptive spillover
#' matrix from your pipeline (e.g. \code{spillover_results_baseline.rds$module_spillover});
#' modules are then ranked by \code{max(abs(.))} across its columns, exactly as in
#' the source analysis. If omitted, the ranking is derived from \code{long} by the
#' mean per-unit spillover per module x timepoint, which is convenient but is NOT
#' guaranteed to match a precomputed cohort-level matrix.
#'
#' @param long Output of prepare_spillover_terms().
#' @param n Number of modules to keep (default 15).
#' @param module_spillover Optional modules x timepoints matrix (rownames = module
#'   labels matching \code{long$module}) used as the ranking statistic. Preferred.
#' @param exclude_times Timepoints to drop from the family (e.g. baseline labels).
#' @return A data.frame(module, timepoint) listing the candidate family.
#' @seealso \code{\link{prepare_spillover_terms}}, \code{\link{fit_spillover_screen}}
#' @export
select_candidates <- function(long, n = 15, module_spillover = NULL,
                              exclude_times = character(0)) {
  if (!is.null(module_spillover)) {
    peak_val <- apply(abs(as.matrix(module_spillover)), 1, max, na.rm = TRUE)
    keep_mod <- names(sort(peak_val, decreasing = TRUE))[seq_len(min(n, length(peak_val)))]
  } else {
    agg  <- stats::aggregate(spill ~ module + timepoint, long,
                             FUN = function(z) mean(z, na.rm = TRUE))
    peak <- stats::aggregate(spill ~ module, agg, FUN = function(z) max(abs(z)))
    keep_mod <- peak$module[order(-peak$spill)][seq_len(min(n, nrow(peak)))]
  }
  tps <- setdiff(unique(long$timepoint), exclude_times)
  expand.grid(module = keep_mod, timepoint = tps, stringsAsFactors = FALSE)
}


#' Baseline-referenced per-timepoint effect vectors (limma contrasts)
#'
#' Builds one gene-level effect vector per post-baseline timepoint: each
#' timepoint contrasted against the POOLED baseline mean, via a limma model with
#' repeated-measures blocking. These are the descriptive effects propagated to
#' rank candidate modules. Reproduces spillover_baseline_variant.R.
#'
#' @param expr Gene x sample expression matrix (responsive-gene scoped).
#' @param group Factor of timepoints (length ncol(expr)).
#' @param block Optional repeated-measures factor (e.g. animal).
#' @param baseline_levels Levels of `group` that define baseline (default Pre1, Pre2).
#' @return Named list of effect vectors, one per post-baseline timepoint.
#' @seealso \code{\link{descriptive_module_spillover}}
#' @export
baseline_effect_vectors <- function(expr, group, block = NULL,
                                    baseline_levels = c("Pre1", "Pre2")) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  tp <- droplevels(factor(group)); lv <- levels(tp)
  design <- stats::model.matrix(~ 0 + tp); colnames(design) <- lv
  if (!is.null(block)) {
    block <- factor(block)
    cc  <- limma::duplicateCorrelation(expr, design, block = block)$consensus.correlation
    fit <- limma::lmFit(expr, design, block = block, correlation = cc)
  } else fit <- limma::lmFit(expr, design)
  pre  <- intersect(baseline_levels, lv); post <- setdiff(lv, pre)
  pre_term <- paste0("(", paste(pre, collapse = " + "), ")/", length(pre))
  cn <- paste0(post, " - ", pre_term)
  contr <- limma::makeContrasts(contrasts = cn, levels = design)
  fit2  <- limma::eBayes(limma::contrasts.fit(fit, contr))
  eff <- lapply(seq_along(cn), function(j) stats::setNames(fit2$coefficients[, j], rownames(fit2)))
  stats::setNames(eff, post)
}


#' Descriptive module-level spillover matrix (candidate-ranking statistic)
#'
#' Propagates each per-timepoint effect vector over the FIXED consensus network
#' with calculate_spillover() and averages gene-level spillover within each module
#' (restricted to modules of at least `min_size` genes). Returns the modules x
#' timepoints matrix that select_candidates() ranks on. This replaces the earlier
#' degenerate single-network use of spillover_timecourse() with a direct,
#' tested calculate_spillover() call; it reproduces spillover_results_baseline.rds
#' $module_spillover.
#'
#' @param network Adjacency matrix or WGCNA-style object (row-normalized internally).
#' @param membership Named gene -> module vector. May be integer ids or "M"-style
#'   labels; rownames of the returned matrix are always "M<k>".
#' @param effects Named list of per-timepoint gene effect vectors
#'   (e.g. from baseline_effect_vectors()).
#' @param min_size Minimum module size to summarize (default 10).
#' @return A modules x timepoints matrix (rownames "M<k>") of mean module spillover.
#' @seealso \code{\link{baseline_effect_vectors}}, \code{\link{select_candidates}}
#' @export
descriptive_module_spillover <- function(network, membership, effects, min_size = 10) {
  W <- extract_network_matrix(network)
  mem <- membership[rownames(W)]
  tab  <- table(mem)
  keep <- names(tab)[tab >= min_size]                 # module ids as-is (integer or "M8")
  num  <- suppressWarnings(as.integer(sub("^M", "", keep)))   # numeric part if present
  keep <- keep[order(ifelse(is.na(num), Inf, num), keep)]     # stable ascending order
  rn   <- ifelse(grepl("^M", keep), keep, paste0("M", keep))  # display rownames always "M<k>"
  tps  <- names(effects)
  M <- matrix(0, length(keep), length(tps), dimnames = list(rn, tps))
  for (t in tps) {
    sv <- calculate_spillover(W, initiator_genes = rownames(W),
                              gene_effects = effects[[t]], row_normalize = TRUE)
    for (mi in seq_along(keep)) M[mi, t] <- mean(sv[names(mem)[mem == keep[mi]]], na.rm = TRUE)
  }
  M
}
