
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QUESTIONS

**QUESTIONS** (*Quantifying Expression Spillover Trajectories in Omic
Network Series*) is an R package for modeling how upstream gene signals
propagate through time-ordered biological networks — especially those
derived from WGCNA or other coexpression analyses.

QUESTIONS extends the spillover propagation framework to help
researchers quantify how biological influence shifts over time, identify
key signal-carrying modules, and interpret dynamic omics data with
network-level context.

------------------------------------------------------------------------

## Features

- **Spillover modeling** between initiator genes and their neighbors
- **Timecourse analysis** across a list of network snapshots
- **Network pruning** via thresholds or top-k neighbors
- **Spillover thresholding** informed by network structure
- **Sample-level scoring** based on spillover-weighted expression
- **Module-level summaries** across timepoints
- Native support for sparse matrices and WGCNA objects

------------------------------------------------------------------------

## Installation

Install the development version of `QUESTIONS` from GitHub with
[pak](https://pak.r-lib.org):

``` r
# install.packages("pak")
pak::pak("CarlyBobak/QUESTIONS")
```

## Basic Example

To compute spillover:

1.  Create a similarity matrix (e.g., TOM or adjacency matrix).
2.  Choose one or more initiator genes to seed the propagation.
3.  Compute spillover using `calculate_spillover()` with optional
    pruning or normalization.
4.  Interpret or summarize results at the gene, sample, or module level.

For example, you might start with a 2×2 matrix of gene similarity
values, initialize spillover from one gene, and observe how signal
distributes across the network. Larger networks and timecourse inputs
can be handled using `spillover_timecourse()`.

Sample- and module-level summaries can be generated using
`summarize_spillover_per_sample()` or
`summarize_spillover_by_sample_and_module()`.

``` r
library(QUESTIONS)
library(Matrix)

mat <- Matrix(c(1, 0.3, 0.3, 1), nrow = 2, sparse = TRUE)
rownames(mat) <- colnames(mat) <- c("Gene1", "Gene2")

calculate_spillover(
  network = mat,
  initiator_genes = "Gene1",
  alpha = 0.9
)
```

## 🖋 Citation

If you use this package, please cite:

Carly Bobak (2025).  
*QUESTIONS: Quantifying Expression Spillover Trajectories in Omic
Network Series*.  
R package version 0.1.0.

------------------------------------------------------------------------

## Author & Development

Developed by Carly Bobak in Research Computing and Data, Jaini Shah, and
James O’Malley in The Dartmouth Institute at Dartmouth College to
support time-resolved and network-informed omics analysis.

Questions and contributions welcome via [GitHub
Issues](https://github.com/CarlyBobak/QUESTIONS/issues).

------------------------------------------------------------------------

## License

MIT License © 2025 Carly Bobak
