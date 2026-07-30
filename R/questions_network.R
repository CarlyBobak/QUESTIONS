# =============================================================================
# QUESTIONS: network construction (tested core)
# Responsive-gene scoping + signed mutual-kNN consensus network + Infomap modules.
# Extracted from the locked analysis scripts (de_select_and_sweep.R,
# build_final_network.R) so the network that spillover propagates over is built by
# documented package functions. These reproduce spillover_network.rds exactly at
# the published settings (FDR 0.01, signed, k = 50, Infomap).
# =============================================================================


#' Select infection- (or condition-) responsive genes by a moderated F-test
#'
#' Fits a limma linear model for change across a grouping factor (e.g. timepoint),
#' accommodating repeated measures within a unit via duplicateCorrelation blocking,
#' and returns the genes passing an FDR threshold on the moderated F-test. Scoping
#' the network to responsive genes is what lets the co-expression graph modularize;
#' the full transcriptome does not.
#'
#' @param expr Gene x sample expression matrix (rows = genes).
#' @param group Factor (length ncol(expr)) whose levels define the conditions to
#'   test change across (e.g. timepoint). Coefficients are each level vs the first.
#' @param block Optional factor (length ncol(expr)) for repeated-measures blocking
#'   (e.g. animal); passed to limma::duplicateCorrelation. NULL for none.
#' @param fdr FDR threshold on the moderated F-test adj.P.Val (default 0.01).
#' @param verbose Logical.
#' @return A list: de_genes (character), expr_de (responsive-gene submatrix),
#'   top (full limma topTable), fdr, consensus_correlation (or NA).
#' @seealso \code{\link{build_knn_network}}
#' @export
select_responsive_genes <- function(expr, group, block = NULL, fdr = 0.01, verbose = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  v <- apply(expr, 1, stats::var); expr <- expr[v > 0, , drop = FALSE]      # drop zero-variance genes
  group <- droplevels(factor(group))
  design <- stats::model.matrix(~ group)
  cc <- NA_real_
  if (!is.null(block)) {
    block <- factor(block)
    corfit <- limma::duplicateCorrelation(expr, design, block = block)
    cc <- corfit$consensus.correlation
    if (verbose) message(sprintf("consensus within-block correlation: %.3f", cc))
    fit <- limma::lmFit(expr, design, block = block, correlation = cc)
  } else {
    fit <- limma::lmFit(expr, design)
  }
  fit <- limma::eBayes(fit)
  tt  <- limma::topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = "none")  # F over levels
  de_genes <- rownames(tt)[tt$adj.P.Val < fdr]
  if (verbose) message(sprintf("Selected %d responsive genes (FDR < %.2f)", length(de_genes), fdr))
  list(de_genes = de_genes, expr_de = expr[de_genes, , drop = FALSE],
       top = tt, fdr = fdr, consensus_correlation = cc)
}


#' Signed mutual-kNN consensus co-expression network with Infomap modules
#'
#' Builds the propagation substrate: correlation similarity (signed or absolute),
#' a mutual (reciprocal) k-nearest-neighbor graph, and Infomap community detection.
#' The reciprocal rule keeps an edge only when each gene is in the other's
#' k-neighborhood, which prunes one-sided ties and yields a sparse, modular graph
#' that preserves between-module bridges. Reproduces build_final_network.R.
#'
#' @param expr Gene x sample expression matrix (already scoped to responsive genes).
#' @param k Neighbors per gene before the mutual filter (default 50).
#' @param sim "signed" (positive correlations only) or "abs" (default "signed").
#' @param block_size Column-block size for the correlation scan (default 500).
#' @param n_cores Cores for the block scan via parallel::mclapply (default 1).
#' @param seed Optional seed (Infomap is randomized).
#' @param verbose Logical.
#' @return A list matching spillover_network.rds: W_graph (sparse adjacency),
#'   membership (named integer module vector), genes, k, sim, n_samples,
#'   module_sizes, plus diagnostics (modularity, within_fraction, singleton_fraction).
#' @seealso \code{\link{select_responsive_genes}}, \code{\link{prepare_spillover_terms}}
#' @export
build_knn_network <- function(expr, k = 50, sim = c("signed", "abs"),
                              block_size = 500, n_cores = 1, seed = NULL, verbose = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")
  sim <- match.arg(sim)
  X <- as.matrix(expr); p <- nrow(X); n <- ncol(X); gname <- rownames(X)
  Zs <- t(scale(t(X))); ZsT <- t(Zs)
  blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))
  top_idx <- function(vv, kk) {
    th <- sort(vv, partial = length(vv) - kk)[length(vv) - kk]
    idx <- which(vv > th); idx[order(vv[idx], decreasing = TRUE)][seq_len(min(kk, length(idx)))]
  }
  proc <- function(b) {
    C <- (Zs[b, , drop = FALSE] %*% ZsT) / (n - 1)
    rv <- if (sim == "signed") C else abs(C)
    for (r in seq_along(b)) rv[r, b[r]] <- -Inf
    if (sim == "signed") rv[rv <= 0] <- -Inf
    ii <- jj <- ww <- vector("list", length(b))
    for (r in seq_along(b)) {
      ix <- top_idx(rv[r, ], k); ix <- ix[is.finite(rv[r, ix])]
      ii[[r]] <- rep.int(b[r], length(ix)); jj[[r]] <- ix; ww[[r]] <- abs(C[r, ix])
    }
    list(i = unlist(ii), j = unlist(jj), w = unlist(ww))
  }
  res <- if (n_cores > 1) parallel::mclapply(blocks, proc, mc.cores = n_cores) else lapply(blocks, proc)
  di <- unlist(lapply(res, `[[`, "i")); dj <- unlist(lapply(res, `[[`, "j")); dw <- unlist(lapply(res, `[[`, "w"))

  mutual <- (dj * (p + 1) + di) %in% (di * (p + 1) + dj)          # reciprocal-edge filter
  mi <- di[mutual]; mj <- dj[mutual]; mw <- dw[mutual]
  a <- pmin(mi, mj); b <- pmax(mi, mj); u <- !duplicated(a * (p + 1) + b)
  a <- a[u]; b <- b[u]; w <- mw[u]
  M <- Matrix::sparseMatrix(i = c(a, b), j = c(b, a), x = c(w, w), dims = c(p, p),
                            dimnames = list(gname, gname))

  g  <- igraph::graph_from_adjacency_matrix(M, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (!is.null(seed)) set.seed(seed)
  cm  <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
  mem <- igraph::membership(cm); deg <- igraph::degree(g)
  within <- sum(mem[a] == mem[b]); between <- sum(mem[a] != mem[b])
  csize  <- sort(as.integer(table(mem)), decreasing = TRUE)
  modq   <- igraph::modularity(g, mem, weights = igraph::E(g)$weight)
  wfrac  <- within / (within + between); sfrac <- mean(deg == 0)
  if (verbose) message(sprintf("edges %d | within %.1f%% | isolated %.1f%% | modularity %.3f",
                               length(a), 100 * wfrac, 100 * sfrac, modq))
  list(W_graph = M, membership = mem, genes = gname, k = k, sim = sim, n_samples = n,
       module_sizes = csize, modularity = modq,
       within_fraction = wfrac, singleton_fraction = sfrac)
}
