# =============================================================================
# QUESTIONS: gene-network spillover (consolidated)
#
# Changes from the previous version:
#   - Removed the alpha damping parameter everywhere. Spillover is now the
#     parameter-free neighbor-weighted exposure: spillover = W %*% effects.
#   - Replaced global-threshold prune_matrix() with a disparity-filter backbone
#     (Serrano, Boguna & Vespignani, PNAS 2009). Its significance parameter is
#     named signif_level (not alpha) to avoid confusion with the removed knob.
#   - Development functions (spillover_timecourse, build_ggm_network, and the
#     multi-timepoint summaries) moved to questions_development.R.
#   - Fixed a trailing-comma call and the net$TOM access that errored on bare
#     matrices (now extracts the matrix before reading gene names).
# =============================================================================


# -----------------------------------------------------------------------------
# Internal: disparity-filter backbone (Serrano et al. 2009)
# -----------------------------------------------------------------------------
# For each node i of degree k_i and strength s_i, the normalized weight of edge
# (i,j) is p_ij = w_ij / s_i. Under the null that a node's normalized weights are
# a uniform random partition of [0,1], the edge p-value from i's side is
# (1 - p_ij)^(k_i - 1). An undirected edge is kept if it is significant from
# EITHER endpoint (the standard multiscale "OR" rule), which preserves
# significant edges at both hubs and the periphery.
.disparity_backbone <- function(W, signif_level = 0.05) {
  sm <- Matrix::summary(W)                       # i, j, x triplets (both triangles)
  s  <- Matrix::rowSums(W)                        # node strengths
  k  <- Matrix::rowSums(W > 0)                    # node degrees

  # directed p-value per stored entry; degree-1 nodes cannot declare significance
  pv <- ifelse(k[sm$i] > 1, (1 - sm$x / s[sm$i])^(k[sm$i] - 1), 1)

  # undirected keep decision: significant from either side (integer edge key)
  n   <- nrow(W)
  key <- (pmin(sm$i, sm$j) - 1) * n + pmax(sm$i, sm$j)
  min_pv   <- tapply(pv, key, min)
  keep_key <- as.numeric(names(min_pv))[min_pv < signif_level]
  keep     <- key %in% keep_key                   # symmetric: both directed entries kept together

  Matrix::sparseMatrix(i = sm$i[keep], j = sm$j[keep], x = sm$x[keep],
                       dims = dim(W), dimnames = dimnames(W))
}


#' Prune a similarity matrix to its statistically significant backbone
#'
#' Extracts the multiscale disparity-filter backbone (Serrano et al. 2009) of a
#' weighted, symmetric, non-negative similarity matrix (e.g. an absolute-correlation
#' adjacency or TOM). Unlike a global threshold, the disparity filter tests each edge
#' locally against the node's own weight distribution, so it retains significant edges
#' at both hubs and the periphery. A genuine periphery is what makes neighbor-weighted
#' spillover measurable.
#'
#' @param mat A square, symmetric, non-negative similarity matrix (dense or sparse).
#' @param signif_level Backbone significance level (default 0.05). Smaller keeps fewer edges.
#' @param verbose If TRUE, prints edge-count and density diagnostics.
#'
#' @return A pruned sparse dgCMatrix containing only backbone edges.
#' @seealso \code{\link{calculate_spillover}}
#' @export
prune_matrix <- function(mat, signif_level = 0.05, verbose = TRUE) {
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) stop("Input must be a matrix.")
  if (!identical(rownames(mat), colnames(mat)))
    stop("Matrix must have matching row and column names.")

  W <- coerce_to_dgCMatrix(mat)                  # dgCMatrix is general, so both triangles are stored
  if (any(W@x < 0)) stop("Disparity filter requires non-negative weights; take abs() upstream.")
  Matrix::diag(W) <- 0
  W <- Matrix::drop0(W)

  n_before <- length(W@x) / 2
  out <- .disparity_backbone(W, signif_level = signif_level)
  n_after <- length(out@x) / 2

  if (verbose) {
    dens <- n_after / (nrow(W) * (nrow(W) - 1) / 2)
    message(sprintf("Disparity backbone: %d -> %d edges (density %.4f) at signif_level %.3f",
                    n_before, n_after, dens, signif_level))
  }
  out
}


#' Calculate gene-level spillover (neighbor-weighted exposure)
#'
#' Propagates a set of initiator effects one step over a gene similarity network as a
#' neighbor-weighted exposure: spillover = W %*% effects, with W optionally row-normalized
#' to a row-stochastic transition matrix. This is the exposure-mapping quantity shared with
#' the ANSWERS peer-exposure construction. There is no damping parameter: sparsification,
#' not a tuned alpha, controls locality.
#'
#' @param network An adjacency matrix, TOM matrix, or WGCNA-style result object.
#' @param initiator_genes Genes to initialize spillover from (used if gene_effects is NULL).
#' @param gene_effects Optional named numeric vector of initiator effects (e.g. logFC).
#' @param row_normalize Logical. Row-normalize W to row-stochastic before propagating (default TRUE).
#' @param prune Logical. Apply the disparity backbone before propagating (default FALSE).
#' @param prune_args Named list passed to prune_matrix(), e.g. list(signif_level = 0.05).
#'
#' @return A named numeric vector of spillover scores over all genes in the network.
#' @seealso \code{\link{prune_matrix}}
#' @export
calculate_spillover <- function(network,
                                initiator_genes,
                                gene_effects = NULL,
                                row_normalize = TRUE,
                                prune = FALSE,
                                prune_args = list(signif_level = 0.05)) {
  # --- Step 1: Extract adjacency matrix ---
  W <- extract_network_matrix(network)

  # --- Convert to sparse dgCMatrix if not already ---
  W <- coerce_to_dgCMatrix(W)

  if (prune) {
    W <- do.call(prune_matrix, c(list(mat = W), prune_args))
  }

  # --- Step 2: Ensure square matrix and named ---
  gene_names <- rownames(W)
  if (is.null(gene_names) || any(gene_names != colnames(W))) {
    stop("The network matrix must have matching row and column names.")
  }

  # --- Step 3: Initialize effect vector ---
  init_vec <- rep(0, length(gene_names))
  names(init_vec) <- gene_names

  if (!is.null(gene_effects)) {
    gene_effects <- gene_effects[names(gene_effects) %in% gene_names]
    init_vec[names(gene_effects)] <- gene_effects
  } else {
    init_vec[initiator_genes[initiator_genes %in% gene_names]] <- 1
  }

  # --- Step 4: Row-normalize W to make stochastic ---
  if (row_normalize) {
    W <- normalize_rows(W)
  }

  # --- Step 5: Spillover calculation (pure neighbor-weighted exposure) ---
  spillover <- W %*% init_vec

  spill_vec <- as.numeric(spillover)
  names(spill_vec) <- rownames(W)
  return(spill_vec)
}

#' Row-normalize a sparse matrix to a row-stochastic transition matrix
#'
#' Divides each row by its sum so rows sum to one. Rows that sum to zero are given
#' a self-loop on the diagonal to avoid division by zero.
#'
#' @param W A square sparse matrix (genes x genes).
#' @param verbose Logical.
#' @return The row-normalized matrix.
#' @seealso \code{\link{is_row_normalized}}, \code{\link{calculate_spillover}}
#' @export
normalize_rows <- function(W, verbose = TRUE) {

  rs <- Matrix::rowSums(W)

  if (any(rs == 0)) {
    ix <- which(rs == 0)
    W[cbind(ix, ix)] <- 1  # self-loops on the DIAGONAL only (W[ix, ix] would fill the whole ix-by-ix block)
    rs[ix] <- 1            # avoid division by zero
  }

  W_norm <- W / rs  # divides row i by rs[i] (column-major recycling with length(rs) == nrow)

  return(W_norm)
}


#' Check if a sparse matrix is row-normalized
#'
#' This function checks whether each row of a sparse matrix sums to 1,
#' within a specified numeric tolerance. This is useful for confirming
#' that a matrix is suitable for diffusion or spillover operations.
#'
#' @param W A sparse matrix of class `dgCMatrix`.
#' @param tolerance Numeric. Acceptable deviation from 1 for each row sum (default = 1e-6).
#' @param verbose Logical. If TRUE, prints the number of rows failing the check.
#'
#' @return Logical. TRUE if all rows sum to 1 within the specified tolerance, FALSE otherwise.
#'
#' @examples
#' W <- Matrix::Matrix(c(0, 1, 1, 0), 2, 2, sparse = TRUE)
#' is_row_normalized(coerce_to_dgCMatrix(normalize_rows(W)))
#'
#' @export
is_row_normalized <- function(W, tolerance = 1e-6, verbose = TRUE) {
  if (!inherits(W, "dgCMatrix")) {
    stop("Input matrix W must be of class 'dgCMatrix'.")
  }

  rs <- Matrix::rowSums(W)
  diffs <- abs(rs - 1)
  not_ok <- which(diffs > tolerance)

  if (verbose) {
    if (length(not_ok) == 0) {
      message("Matrix appears to be row-normalized.")
    } else {
      message(sprintf("%d rows are not row-normalized (tolerance = %.1e).", length(not_ok), tolerance))
    }
  }

  return(length(not_ok) == 0)
}


#' Extract a matrix from a network input
#'
#' This function handles a network object that may be a raw matrix, sparse matrix,
#' or a WGCNA result list containing `TOM` or `adjacency`. It extracts the appropriate
#' matrix for downstream use.
#'
#' @param network A matrix, sparse Matrix, or a list from WGCNA with `TOM` or `adjacency`.
#'
#' @return A base or sparse matrix suitable for analysis.
#' @export
extract_network_matrix <- function(network) {
  if (is.list(network)) {
    if (!is.null(network$TOM)) {
      return(network$TOM)
    } else if (!is.null(network$adjacency)) {
      return(network$adjacency)
    } else if (!is.null(network$datExpr)) {
      stop("Please provide a TOM or adjacency matrix in the WGCNA object.")
    } else {
      stop("Unrecognized network format. Provide TOM, adjacency matrix, or WGCNA-style list.")
    }
  } else if (inherits(network, "matrix") || inherits(network, "Matrix")) {
    return(network)
  } else {
    stop("Invalid input for 'network'. Must be matrix or list with TOM/adjacency.")
  }
}


#' Coerce a matrix to sparse dgCMatrix format
#'
#' Converts input matrix to a general compressed sparse column matrix (dgCMatrix),
#' with special handling for symmetric sparse matrices.
#'
#' @param W A base matrix, Matrix, or symmetric sparse matrix (e.g., dsCMatrix).
#'
#' @return A `dgCMatrix` object.
#' @importFrom methods as
#' @export
coerce_to_dgCMatrix <- function(W) {
  if (!inherits(W, "dgCMatrix")) {
    # Modern Matrix idiom: route through the virtual superclasses rather than
    # the deprecated direct as(<dsCMatrix>, "dgCMatrix"). Coercing to
    # "CsparseMatrix" gives a column-compressed sparse form, and "generalMatrix"
    # expands any symmetric storage to both triangles. The result is always a
    # general double column-compressed matrix, i.e. a dgCMatrix.
    W <- as(as(W, "CsparseMatrix"), "generalMatrix")
  }
  return(W)
}