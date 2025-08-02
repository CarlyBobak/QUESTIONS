#' Prune a similarity matrix (TOM or adjacency) to enhance sparsity
#'
#' Supports thresholding or top-k nearest neighbors to reduce matrix density,
#' improving interpretability and biological signal in diffusion tasks.
#'
#' @param mat A square similarity matrix (e.g., TOM or adjacency)
#' @param value A numeric threshold (if method = "threshold") or number of neighbors (if method = "topk")
#' @param auto If TRUE, determine threshold automatically to achieve target sparsity
#' @param target_sparsity Desired sparsity level (only used if auto = TRUE and method = "threshold")
#' @param verbose If TRUE, prints pruning diagnostics
#'
#' @return A pruned matrix with zeroed weak connections
#' @export
prune_matrix <- function(mat,
                         value = NULL,
                         auto = FALSE,
                         target_sparsity = 0.05,
                         verbose = TRUE) {
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) stop("Input must be a matrix.")
  if (!identical(rownames(mat), colnames(mat))) stop("Matrix must have matching row and column names.")

  pruned <- mat

  if((is.null(value)&(auto==F))){
    stop("Please provide a threshold value, or set auto=TRUE")
  }

  if (auto) {
    # Total number of elements in full matrix
    total_elements <- nrow(mat)^2  # square matrix

    # Number of non-NA, non-zero entries you want to retain
    target_nonzero_count <- round(target_sparsity * total_elements,0)

    # Get all values (including zeros, if dense) or just non-zero if sparse
    all_vals <- as.vector(mat@x)
    all_vals <- all_vals[!is.na(all_vals)]  # remove NA if any

    # Get cutoff corresponding to retaining `target_nonzero_count` largest values
    sorted_vals <- sort(all_vals, decreasing = TRUE)
    value <- sorted_vals[target_nonzero_count]
    if (verbose) message(sprintf("Auto-threshold selected: %.4f for target sparsity %.2f", value, target_sparsity))
  }

  pruned@x[pruned@x <= value] <- 0  # this zeroes out small values
  pruned <- Matrix::drop0(pruned)

  return(pruned)
}


#' Calculate gene-level spillover from time 1 to time 2
#'
#' This function estimates spillover effects from a set of initiator genes (with optional initial effects)
#' using a transition matrix derived from a TOM matrix, adjacency matrix, or WGCNA object. Optionally, the
#' network matrix can be pruned to remove weak or noisy connections using either a threshold or top-k filtering strategy.
#'
#' @param network An adjacency matrix, TOM matrix, or a WGCNA result object (list with `TOM`, `adjacency`, or `datExpr`).
#' @param initiator_genes A vector of gene names to initialize spillover from. Can represent hub genes, biomarkers, etc.
#' @param gene_effects Optional named numeric vector with effect sizes (e.g., logFC or median shift) for initiators. Default = 1.
#' @param alpha Damping parameter controlling conservation of original signal (default = 1). Set between 0 (no propagation) and 1 (full propagation).
#' @param row_normalize Logical. Whether to row-wise normalize the matrix before computing spillover values. (default = TRUE).
#' @param prune Logical. Whether to prune the network matrix before computing spillover (default = FALSE).
#' @param prune_args A named list of arguments passed to `prune_matrix()`, such as `value`, `auto`, and `target_sparsity`.
#'
#'#' @details
#' **Spillover Calculation and Interpretation**
#'
#' This function estimates gene-level *spillover*, representing the indirect influence of upstream initiator genes
#' as propagated through a gene similarity network (e.g., adjacency or TOM matrix). The propagation follows
#' a dampened diffusion model:
#'
#' \deqn{
#'   \text{spillover} = \alpha \cdot W \cdot \text{initial\_effects} + (1 - \alpha) \cdot \text{initial\_effects}
#' }
#'
#' where:
#' - \eqn{W} is the (optionally normalized) similarity or transition matrix,
#' - \eqn{\text{initial\_effects}} is a vector of user-specified gene-level initiator values (default = 1),
#' - \eqn{\alpha} is a damping parameter that controls how much signal is allowed to propagate versus being retained.
#'
#' The output is a named numeric vector of spillover scores for each gene, reflecting the degree to which
#' each gene was influenced by the initiators, either directly (if it is an initiator itself) or indirectly
#' through connected neighbors in the network.
#'
#' - Higher spillover values indicate greater expected influence from initiator genes.
#' - When `row_normalize = TRUE`, each row of the matrix \eqn{W} is scaled to sum to 1, making the spillover
#'   interpretable as a probabilistic distribution of influence across neighbors.
#' - When `alpha = 1`, only indirect signal is retained; lower alpha retains more direct signal.
#'
#' **Interpretation Notes:**
#' - If pruning is enabled, weak connections are removed before computing spillover. This can sharpen signal
#'   and improve interpretability but may reduce sensitivity to diffuse propagation.
#' - Raw spillover values are continuous, but follow-up steps can apply thresholds or binarization to select
#'   genes with biologically meaningful levels of inferred influence.
#'
#' This function serves as the core engine for timecourse analyses (see `spillover_timecourse()`) and
#' downstream integration with expression and phenotype data.
#'
#' @return A named numeric vector of spillover scores at time 2 for all genes in the network.
#'
#' @seealso \code{\link{prune_matrix}} for preprocessing options to filter weak edges before spillover.
#'
#' @examples
#' \dontrun{
#' spillover <- calculate_spillover(
#'   network = tom_matrix,
#'   initiator_genes = hub_genes,
#'   alpha = 0.85,
#'   prune = TRUE,
#'   prune_args = list(method = "threshold", auto = TRUE, target_sparsity = 0.05)
#' )
#' }
#'
#'@export

calculate_spillover <- function(network,
                                initiator_genes,
                                gene_effects = NULL,
                                alpha = 1, row_normalize=FALSE, prune = FALSE,
                                prune_args = list(value=NULL, auto = TRUE, target_sparsity = 0.05)) {
  # --- Step 1: Extract adjacency matrix ---
  W<-extract_network_matrix(network)

  # --- Convert to sparse dgCMatrix if not already ---
  W<-coerce_to_dgCMatrix(W)

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

  # Use gene_effects if provided, else default to 1
  if (!is.null(gene_effects)) {
    gene_effects <- gene_effects[names(gene_effects) %in% gene_names]
    init_vec[names(gene_effects)] <- gene_effects
  } else {
    init_vec[initiator_genes[initiator_genes %in% gene_names]] <- 1
  }

  # --- Step 4: Row-normalize W to make stochastic ---
  if(row_normalize){
    W<-normalize_rows(W)
  }

  # --- Step 5: Spillover calculation ---
  spillover <- alpha * (W %*% init_vec) + (1 - alpha) * init_vec

  spill_vec <- as.numeric(spillover)
  names(spill_vec) <- rownames(W)
  return(spill_vec)  # return named vector
}

#' Calculate timecourse spillover across multiple networks
#'
#' This function iteratively computes gene-level spillover across a series of networks (e.g., WGCNA results or TOM/adjacency matrices),
#' using a propagation model. Spillover can begin from a shared initiator list or from timepoint-specific gene sets and optionally
#' incorporate gene-level effects such as fold changes or module membership. Propagation can be binary or continuous, and supports
#' persistent initiators and feedback from previous spillover states. The function also supports optional pruning of dense networks
#' using a similarity threshold or top-k strongest connections per row.
#'
#' @param network_list A named list of networks (matrices or WGCNA result objects), one per timepoint (T2 to Tn).
#' @param initiator_genes Either a character vector of shared initiators or a list of vectors, one per timepoint (T1 to Tn-1).
#' @param gene_effects Optional: either a named numeric vector of gene effects (e.g., logFCs) or a list of such vectors per timepoint.
#'                     If not provided, defaults to 1 for initiators and 0 otherwise.
#' @param alpha Numeric. Damping parameter controlling signal retention (default = 1). Set between 0 (no propagation) and 1 (full propagation).
#' @param row_normalize Logical. Whether to row-wise normalize the matrix before computing spillover values. (default = TRUE).
#' @param binary Logical. If TRUE, converts output at each timepoint to binary using a threshold (default = TRUE).
#' @param threshold_spillover Logical. If TRUE only spillove values above threshold will be retained (defualt = FALSE).
#' @param threshold Numeric. Used if `binary = TRUE` or `threshold_spillover = TRUE` to define the activation cutoff (default = 0).
#' @param calc_threshold Logical. If TRUE, uses spillover_threshold_from_degree to calculate a topology informed threshold for significant spillover.
#' @param quantile Numeric. Value between 0 and 1 (default = 0.95). The spilllover threshold is set at this quantile of the expected background distribution.
#' @param persistent_initiators Logical. If TRUE, initiators accumulate over time and remain active after first activation (default = FALSE).
#' @param use_previous_spillover Logical. If TRUE, uses spillover output from the previous timepoint as the gene effects at the next timepoint (default = FALSE).
#' @param prune Logical. Whether to prune the network matrix before computing spillover (default = FALSE).
#' @param prune_args A named list of arguments passed to `prune_matrix()`, such as `method`, `value`, `auto`, and `target_sparsity`.
#' @param verbose Logical. If TRUE, prints diagnostic and progress messages (default = TRUE).
#'
#'#' @details
#' **Overview and Interpretation**
#'
#' This function generalizes the single-step spillover propagation to an entire *timecourse* of network states,
#' allowing users to track how influence from an initial set of genes evolves across multiple biological snapshots (e.g., timepoints, conditions).
#'
#' For each timepoint, spillover is calculated using the following damped propagation model:
#'
#' \deqn{
#'   \text{spillover}_{t} = \alpha \cdot W_{t} \cdot \text{effects}_{t} + (1 - \alpha) \cdot \text{effects}_{t}
#' }
#'
#' where:
#' - \eqn{W_{t}} is the network matrix at timepoint \eqn{t}, optionally pruned and row-normalized,
#' - \eqn{\text{effects}_{t}} is a gene-level effect vector derived from initiators or the previous spillover state,
#' - \eqn{\alpha} controls how much signal propagates vs. remains locally conserved.
#'
#' **Options and Flexibility**
#'
#' - **Initiator Genes**: The function supports either a single set of initiators or a list for per-timepoint variation.
#' - **Persistent Initiators**: If enabled, genes that were active in any prior step remain active in subsequent ones.
#' - **Previous Spillover as Input**: Allows dynamic feedback, where the spillover at \eqn{t-1} seeds the effects at \eqn{t}.
#' - **Thresholding**: Optional filtering of weak spillover using either a fixed numeric cutoff or a topology-aware quantile threshold from `spillover_threshold_from_degree()`.
#' - **Binary Conversion**: If enabled, spillover values are binarized (1/0) to reflect activation status, useful for logic-based modeling or visualization.
#' - **Pruning**: Supports automatic or manual reduction of network density to enhance signal-to-noise and speed.
#'
#' **Interpretation Guidance**
#'
#' - Raw spillover scores reflect the accumulated influence of initiators on each gene at each timepoint.
#' - These scores can be used to identify genes that consistently accumulate signal, detect pathway shifts, or correlate with phenotypic changes.
#' - Thresholding or binarization is helpful when inferring discrete activation or tracing modules step-by-step.
#' - Use `spillover_threshold_from_degree()` for principled cutoffs based on network structure, especially in row-normalized matrices.
#'
#' This function is particularly well-suited for longitudinal studies (e.g., developmental timecourses, intervention responses),
#' where spillover-based diffusion reveals shifting patterns of gene influence across the network's evolving topology.
#'
#' @return A named list of spillover vectors (named numeric vectors) for each timepoint in the input list.
#'
#' @seealso \code{\link{calculate_spillover}}, \code{\link{prune_matrix}}
#'
#'@export
spillover_timecourse <- function(
    network_list,
    initiator_genes,
    gene_effects = NULL,
    alpha = 1,
    row_normalize = TRUE,
    binary = FALSE,
    threshold_spillover = FALSE,
    threshold = NULL,
    calc_threshold = FALSE,
    quantile = 0.95,
    persistent_initiators = FALSE,
    use_previous_spillover = FALSE,
    prune = FALSE,
    prune_args = list(auto = FALSE, target_sparsity = 0.05),
    verbose = TRUE
) {

  # --- Input checks ---
  if (!is.list(network_list) || is.null(names(network_list))) {
    stop("`network_list` must be a named list of networks/WGCNA results, one per timepoint.")
  }

  if (is.list(initiator_genes)) {
    if (length(initiator_genes) != length(network_list)) {
      stop("If `initiator_genes` is a list, it must have the same length as `network_list`.")
    }
  }

  if (!is.list(initiator_genes) & is.null(gene_effects) & !use_previous_spillover) {
    warning("Using the same initiator genes at every timepoint with binary effects and no prior spillover. If this is unintended, consider passing `gene_effects` or enabling `use_previous_spillover = TRUE`.")
  }

  if (!is.null(gene_effects)) {
    if (is.list(gene_effects)) {
      if (length(gene_effects) != length(network_list)) {
        stop("If `gene_effects` is a list, it must have the same length as `network_list`.")
      }
    } else {
      if (is.list(initiator_genes) && length(gene_effects) != 1) {
        stop("If `initiator_genes` is a list, `gene_effects` must be a list of the same length or a single shared vector.")
      }
      if (!use_previous_spillover) {
        stop("If `gene_effects` is a single vector and `use_previous_spillover = FALSE`, it cannot be reused across timepoints.")
      }
    }
  }

  spillover_results <- list()
  n_timepoints <- length(network_list)

  # --- Track persistent initiators if enabled ---
  all_active_genes <- character(0)

  for (i in seq_len(n_timepoints)) {
    if (verbose) message("Processing timepoint ", i, "/", n_timepoints)
    net <- network_list[[i]]

    # --- Determine initiators and effects ---
    if (use_previous_spillover && i > 1) {
      # Use previous spillover directly
      if (binary) {
        genes_on <- names(prev_spill)[prev_spill > threshold]
        effects <- rep(1, length(genes_on))
        names(effects) <- genes_on
      } else {
        effects <- prev_spill
        genes_on <- names(effects)[effects != 0]
      }
    } else {
      # Use user-provided initiators/effects
      if (is.list(initiator_genes)) {
        genes_on <- initiator_genes[[i]]
      } else {
        genes_on <- initiator_genes
      }

      if (persistent_initiators) {
        all_active_genes <- union(all_active_genes, genes_on)
        genes_on <- all_active_genes
      }

      if (is.null(gene_effects)) {
        all_genes <- rownames(net$TOM %||% net$adjacency %||% net)
        effects <- rep(0, length(all_genes))
        names(effects) <- all_genes
        effects[intersect(genes_on, all_genes)] <- 1
      } else if (is.list(gene_effects)) {
        effects <- gene_effects[[i]]
      } else {
        effects <- gene_effects
      }
    }

    net<-extract_network_matrix(net)

    if(prune){
      net <- do.call(prune_matrix, c(list(mat = net, verbose = verbose), prune_args))
    }

    if(row_normalize){
      net<-normalize_rows(net)
    }

    if((row_normalize==F)&calc_threshold){
      warning("Normalizing rows of adjacency matrix for spillover threshold calculation")
      net<-normalize_rows(net)
    }

    # --- Run spillover ---
    spill <- calculate_spillover(
      network = net,
      initiator_genes = genes_on,
      gene_effects = effects,
      alpha = alpha,
      row_normalize = FALSE,
      prune = FALSE,
    )

    # --- Optional thresholding meaningful spillover---
    if(threshold_spillover){
      if(calc_threshold){
        threshold_args <- list(W = net, quantile = quantile, verbose = verbose)
        threshold <- do.call(spillover_threshold_from_degree, threshold_args)
      } else if(is.null(threshold)){
        stop("Please provide threshold for spillover activation cutoff.")
      }
      spill[spill <= threshold]<-0
    }

    spillover_results[[i]] <- spill

    # --- Optional binary thresholding for propagation---
    if (binary) {
      spill <- ifelse(spill > threshold, 1, 0)
    }

    prev_spill <- spill
  }

  names(spillover_results) <- names(network_list)
  return(spillover_results)
}

#' Calculate a degree-aware spillover threshold for row-normalized networks
#'
#' This function estimates a network-specific threshold for spillover activation,
#' based on the expected signal that each node would receive from random initiator sets.
#' It uses the assumption that, in a row-normalized matrix, initiators distribute
#' 1 unit of signal across their neighbors, and thus the expected spillover
#' is proportional to the degree (number of connections) of each node.
#'
#' This threshold can be used to distinguish spillover values that exceed
#' what is expected from random signal diffusion alone, and is intended for
#' use with row-normalized sparse adjacency or TOM matrices.
#'
#' @param W A row-normalized sparse network matrix of class `dgCMatrix`, where rows and columns have matching gene names.
#' @param quantile Numeric value between 0 and 1 (default = 0.95). The threshold is set at this quantile of the expected background distribution.
#' @param verbose Logical. If TRUE, prints the computed threshold and quantile information.
#' @param tolerance Numeric tolerance for row-sum deviation (default = 1e-6).
#'
#' @return A single numeric value representing the recommended spillover threshold.
#'
#' @examples
#' \dontrun{
#' threshold <- spillover_threshold_from_degree(W = my_network, quantile = 0.95)
#' spill_binary <- ifelse(spillover > threshold, 1, 0)
#' }
#'
#' @export
spillover_threshold_from_degree <- function(W, quantile = 0.95, verbose = TRUE, tolerance = 1e-6) {
  if (!inherits(W, "dgCMatrix")) {
    stop("Input matrix W must be of class 'dgCMatrix'. Please convert before running.")
  }

  # --- Check row-normalization ---
  rs <- Matrix::rowSums(W)
  if (any(abs(rs - 1) > tolerance)) {
    warning("Normalizing network matrix")
    W<-normalize_rows(W)
  }

  # --- Compute degree-aware expected signal ---
  degrees <- Matrix::rowSums(W != 0)
  expected_spillover <- degrees / sum(degrees)

  # --- Compute threshold ---
  threshold <- stats::quantile(expected_spillover, probs = quantile, names = FALSE)

  if (verbose) {
    message(sprintf("Spillover threshold set at %.6f (%.0fth percentile of degree-based expectation)",
                    threshold, quantile * 100))
  }

  return(threshold)
}

#' Row-wise normalize a sparse matrix
#'
#' This function normalizes each row of a sparse matrix so that the row sums to 1.
#' This is often required before diffusion or spillover modeling to create a row-stochastic matrix.
#'
#' If any row has a sum of zero (i.e., disconnected node), the function replaces the diagonal entry
#' for that row with 1, effectively preserving its identity and avoiding division by zero.
#'
#' @param W A sparse matrix of class `dgCMatrix`, typically an adjacency or TOM matrix.
#' @param verbose Logical. If TRUE, prints the number of rows that were corrected (i.e., had a sum of zero).
#'
#' @return A row-normalized sparse matrix of the same dimensions and class.
#'
#' @examples
#' \dontrun{
#' W_norm <- normalize_rows(W)
#'}
#' @export
normalize_rows <- function(W, verbose = TRUE) {

  rs <- Matrix::rowSums(W)

  if (any(rs == 0)) {
    ix <- which(rs == 0)
    W[ix, ix] <- 1  # insert self-edges
    rs[ix] <- 1     # avoid division by zero
  }

  W_norm <- W / rs  # performs row-wise broadcasting

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
#' \dontrun{
#' is_row_normalized(W)
#'}
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
#' @export
coerce_to_dgCMatrix <- function(W) {
  if (!inherits(W, "dgCMatrix")) {
    if (inherits(W, "dsCMatrix")) {
      # Expand symmetric to general dense matrix first
      W <- as(Matrix::forceSymmetric(W, uplo = "U"), "generalMatrix")
      W <- Matrix::Matrix(W, sparse = TRUE)
    } else {
      W <- Matrix::Matrix(W, sparse = TRUE)
    }
    W <- methods::as(W, "dgCMatrix")  # Coerce explicitly from general sparse matrix
  }
  return(W)
}

#' Compute sample-level spillover summary scores
#'
#' Aggregates expression data using gene-level spillover weights. Accepts either a single matrix/vector
#' or lists of both for timecourse-style analysis (e.g., paired with `spillover_timecourse()`).
#'
#' @param expr Either a matrix/data.frame (genes × samples) or a list of such matrices.
#' @param spillover Either a named numeric vector of spillover scores or a list of such vectors.
#' @param summary_fn Optional custom function: takes two arguments (expr, spillover) and returns a scalar.
#'                   Default is a weighted sum of expr × spillover.
#' @param normalize_scores Logical. Whether to normalize spillover scores (ie return spillover Z scores) (default = FALSE).
#' @param genes_in_rows Logical. TRUE if rows are genes and columns are samples (default = FALSE).
#' @param abs_expression Logical. If TRUE, take absolute value of expression before summarizing (default = FALSE).
#' @param verbose Logical. Whether to print progress and diagnostics (default = TRUE).
#'
#'#' @details
#' **Interpreting Per-Sample Spillover Scores**
#'
#' This function computes a summary score for each sample, reflecting how strongly the sample expresses genes
#' that were influenced by network-based spillover from upstream signals (e.g., hub genes or initiators).
#'
#' The default scoring function calculates a weighted sum:
#'
#' \deqn{
#'   \text{score}_s = \sum_{g \in G} \text{expression}_{gs} \times \text{spillover}_g
#' }
#'
#' where:
#' - \eqn{G} is the set of genes with non-zero spillover
#' - \eqn{\text{expression}_{gs}} is the expression of gene *g* in sample *s*
#' - \eqn{\text{spillover}_g} is the spillover score assigned to gene *g* at that timepoint
#'
#' This score can be interpreted as the degree to which a sample expresses the genes that received
#' signal via network propagation. Higher values suggest a stronger manifestation of upstream network activity.
#'
#' **Z-score Option (`normalize_samples = TRUE`):**
#' When this is enabled, the resulting spillover scores are standardized across samples:
#'
#' \deqn{
#'   Z_s = \frac{\text{score}_s - \mu}{\sigma}
#' }
#'
#' where \eqn{\mu} and \eqn{\sigma} are the mean and standard deviation of spillover scores across all samples.
#' Z-scores are useful for comparing samples across different datasets or timepoints, or when integrating with other predictors.
#'
#' **Absolute Expression Option (`abs_expression = TRUE`):**
#' This flag takes the absolute value of expression values before applying the spillover weights. This helps avoid
#' signal cancellation if the expression data has been centered, scaled, or otherwise contains negative values.
#' It is especially useful when the presence or magnitude of expression is of interest, rather than its direction.
#'
#' @return A named numeric vector (if single input) or named list of per-timepoint score vectors.
#'
#' @examples
#' \dontrun{
#' # Single timepoint
#' scores <- summarize_spillover_per_sample(expr, spillover_vec)
#'
#' # Timecourse
#' scores_list <- summarize_spillover_per_sample(expr_list, spillover_list)
#'}
#' @export
summarize_spillover_per_sample <- function(expr,
                                           spillover,
                                           summary_fn = NULL,
                                           normalize_scores = FALSE,
                                           genes_in_rows = FALSE,
                                           abs_expression = FALSE,
                                           verbose = TRUE) {
  # --- Handle single vs list input ---
  is_timecourse <- is.list(expr) && is.list(spillover)

  if (is_timecourse) {
    if (length(expr) != length(spillover)) {
      stop("If both `expr` and `spillover` are lists, they must be the same length.")
    }

    result_list <- vector("list", length(expr))
    names(result_list) <- names(expr)

    for (i in seq_along(expr)) {
      if (verbose) message(sprintf("Scoring timepoint %d/%d: %s", i, length(expr), names(expr)[i]))
      result_list[[i]] <- summarize_spillover_per_sample(expr[[i]],
                                                         spillover[[i]],
                                                         summary_fn = summary_fn,
                                                         normalize_scores = normalize_scores,
                                                         genes_in_rows = genes_in_rows,
                                                         abs_expression = abs_expression,
                                                         verbose = verbose)
    }

    return(result_list)
  }

  # --- Single matrix/vector case ---
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("`expr` must be a matrix/data frame or list thereof.")
  }
  if (!is.numeric(spillover)) stop("`spillover` must be a numeric named vector or list thereof.")

  if (genes_in_rows) expr <- t(expr)

  gene_overlap <- intersect(colnames(expr), names(spillover))

  if (length(gene_overlap) < 5) {
    warning("Fewer than 5 overlapping genes between expression data and spillover vector.")
  }

  expr_sub <- expr[, gene_overlap, drop = FALSE]
  spillover_sub <- spillover[gene_overlap]

  if (abs_expression) {
    expr_sub <- abs(expr_sub)
  }

  if (is.null(summary_fn)) {
    scores <- as.numeric(as.matrix(expr_sub) %*% spillover_sub)
  } else {
    scores <- apply(t(expr_sub), 2, function(x) summary_fn(x, spillover_sub))
  }

  if(normalize_scores){
    scores <- scale(scores)[,1]
  }

  names(scores)<-rownames(expr_sub)

  return(scores)
}

#' Summarize sample-level spillover per module across a timecourse
#'
#' This function calculates a spillover summary score for each sample, grouped by gene module,
#' across a timecourse of expression matrices and spillover vectors. It supports optional
#' absolute-value transformation of expression and z-score normalization of the final result.
#'
#' This enables downstream analyses such as:
#' - Tracing how module-specific signal propagates through time
#' - Associating spillover magnitude with phenotypes
#' - Contrasting signal distribution across modules and samples
#'
#' @param expr_list A list of gene expression matrices, one per timepoint (samples × genes).
#'                  Must be in the same order as `spillover_list`.
#' @param spillover_list A list of named numeric vectors of spillover scores, one per timepoint.
#'                       Names must match gene names in expression and module inputs.
#' @param module_membership_list A list of named vectors mapping genes to modules (character or factor), one per timepoint.
#' @param summary_func Optional. A function to compute the summary statistic (default = sum(expr × spillover)).
#'                     The function must accept two numeric vectors: expression and spillover.
#' @param genes_in_rows Logical. TRUE if rows are genes and columns are samples (default = FALSE).
#' @param abs_expression Logical. If TRUE, take absolute value of expression before applying `summary_func` (default = FALSE).
#' @param normalize_scores Logical. If TRUE, z-score standardize the spillover scores across samples (per module) (default = FALSE).
#' @param verbose Logical. If TRUE, print progress messages (default = TRUE).
#'
#' @details
#' **Interpreting Spillover Summary Values**
#'
#' This function calculates a summary score for each sample that reflects the *weighted influence*
#' of genes receiving spillover signal, aggregated within each module. The default score is computed as:
#'
#' \deqn{
#'   \text{score}_{is} = \sum_{g \in M_i} \text{expression}_{sg} \times \text{spillover}_g
#' }
#'
#' where:
#' - \eqn{M_i} is the set of genes in module *i*
#' - \eqn{\text{expression}_{sg}} is the expression of gene *g* in sample *s*
#' - \eqn{\text{spillover}_g} is the spillover score assigned to gene *g* at that timepoint
#'
#' The result reflects how strongly a sample expresses genes that were targeted by spillover, and thus may
#' indicate whether the biological effects of an upstream perturbation are active in that sample.
#'
#' **Z-score Option (`normalize_scores = TRUE`):**
#' When enabled, scores are standardized *per module* (i.e., column-wise) across all samples:
#'
#' \deqn{
#'   Z_{is} = \frac{\text{score}_{is} - \mu_i}{\sigma_i}
#' }
#'
#' where \eqn{\mu_i} and \eqn{\sigma_i} are the mean and standard deviation of scores for module *i* across samples.
#' Z-scoring is useful when comparing modules with different dynamic ranges, or when integrating scores across timepoints.
#'
#' **Absolute Expression Option (`abs_expression = TRUE`):**
#' When this option is set, gene expression values are converted to absolute values before weighting. This can help
#' avoid cancellation effects if the expression matrix has been centered or includes negative values (e.g., z-scores or residuals).
#' Use this option when the biological interpretation of expression is directional (e.g., presence/absence) rather than fold-change.
#'
#' @return A list of data frames, one per timepoint. Each data frame is samples × modules, with spillover summary scores.
#'
#' @examples
#' \dontrun{
#' # Basic usage with z-score
#' module_scores <- summarize_spillover_by_sample_and_module(expr_list, spillover_list, module_list,
#'                    normalize_scores = TRUE)
#'
#' # Use absolute expression with custom summary
#' module_scores <- summarize_spillover_by_sample_and_module(expr_list, spillover_list, module_list,
#'                    abs_expression = TRUE,
#'                    summary_func = function(expr, spill) mean(expr * spill))
#'}
#' @export
summarize_spillover_by_sample_and_module <- function(expr_list,
                                                     spillover_list,
                                                     module_membership_list,
                                                     summary_func = NULL,
                                                     genes_in_rows = FALSE,
                                                     abs_expression = FALSE,
                                                     normalize_scores = FALSE,
                                                     verbose = TRUE) {
  stopifnot(length(expr_list) == length(spillover_list),
            length(expr_list) == length(module_membership_list))

  if (is.null(summary_func)) {
    summary_func <- function(expr, spill) sum(expr * spill, na.rm = TRUE)
  }

  n_timepoints <- length(expr_list)
  result_list <- vector("list", n_timepoints)
  names(result_list) <- names(expr_list)

  for (i in seq_len(n_timepoints)) {
    if (verbose) message("Processing timepoint ", i, "/", n_timepoints)

    expr <- expr_list[[i]]
    spill <- spillover_list[[i]]
    module_map <- module_membership_list[[i]]

    if (genes_in_rows) {
      expr <- t(expr)
    }

    common_genes <- Reduce(intersect, list(colnames(expr), names(spill), names(module_map)))
    if (length(common_genes) < 5) {
      warning(sprintf("Only %d overlapping genes at timepoint %s", length(common_genes), names(expr_list)[i]))
    }

    expr <- expr[, common_genes, drop = FALSE]
    spill <- spill[common_genes]
    module_map <- module_map[common_genes]

    if (abs_expression) {
      expr <- abs(expr)
    }

    modules <- unique(module_map)
    samples <- rownames(expr)

    summary_mat <- matrix(NA_real_, nrow = length(samples), ncol = length(modules),
                          dimnames = list(samples, modules))

    for (m in modules) {
      genes_m <- names(module_map[module_map == m])
      expr_m <- expr[, genes_m, drop = FALSE]
      spill_m <- spill[genes_m]

      for (s in samples) {
        val <- summary_func(expr_m[s, ], spill_m)
        if (length(val) != 1 || !is.numeric(val)) {
          stop("`summary_func` must return a single numeric value per sample-module combination.")
        }
        summary_mat[s, m] <- val
      }
    }

    df <- as.data.frame(summary_mat)
    if (normalize_scores) {
      df <- as.data.frame(scale(df))
    }

    result_list[[i]] <- df
  }

  return(result_list)
}



