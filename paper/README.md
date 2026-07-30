# Reproducing the QUESTIONS tuberculosis case study

This directory reproduces the applied results in the QUESTIONS manuscript: the
spillover-outcome analysis of the cynomolgus macaque *Mycobacterium tuberculosis*
model (GEO accession GSE84152), including Table 1 and the Figure 1 heatmaps. It
drives the analysis through the installed QUESTIONS package, so it also serves as
a worked example of the package on real data.

The associations reported are correlational, not causal (the network is an
undirected co-expression graph).

If you only want to apply the method to your own data, you do not need this
directory; see the package README.

## Requirements

- The QUESTIONS package installed (`devtools::install_github("CarlyBobak/QUESTIONS")`).
- Bioconductor packages used by the reproduction (not all are package
  dependencies): `limma`, `sva`, `GEOquery`, `WGCNA`. Install via:
  ```r
  BiocManager::install(c("limma", "sva", "GEOquery", "WGCNA"))
  ```
  (`WGCNA` is needed only by the data-loading step, for probe-to-gene collapse.)
- CRAN helpers for loading: `data.table`, `dplyr`, `tidyr`, `stringr`, `tibble`, `ggplot2`.

## Data

The analysis runs from staged intermediate objects:

- `spillover_network.rds`  consensus co-expression network + module assignment
- `de_selection.rds`       responsive-gene expression matrix
- `week1_data.rds`         sample metadata + windowed expression (`meta_win`, `expr_all_win`)

These are provided as a release asset (see the repository Releases page). They can
also be regenerated from the public raw data; see "Regenerating from GEO" below.

## Run order

From this directory, with the package installed:

```r
# 1. (optional) regenerate the staged data from GEO
source("data-raw/build_from_geo.R")        # writes week1_data.rds

# 2. (optional) rebuild the network from the responsive-gene expression
#    (reproduce_table1.R loads the frozen network by default; rebuilding is only
#     needed to reconstruct it from scratch, and Infomap is randomized, so module
#     numbering will differ from the frozen object)
#    See the optional block inside reproduce_table1.R.

# 3. reproduce the primary result (Table 1)
source("reproduce_table1.R")               # prints Table 1, writes table1_reproduced.csv

# 4. reproduce the descriptive heatmaps (Figure 1)
source("spillover_heatmaps.R")             # writes fig_spillover_baseline/increment (pdf+png)
```

## What you should see

`reproduce_table1.R` prints a five-row Table 1. Spillover into the myeloid modules
M8 (day 42, month 3) and M91 (day 30, day 42), and into the platelet module M6
(day 56), clears Benjamini-Hochberg adjustment across the 105-test family, all in
the peak-response window. Coefficients and direct-term p-values are deterministic
and reproduce exactly; permutation and BH values are from a randomized test and
will vary slightly run to run (the fixed seed and settings match the reported
run). `p_direct` is reported unadjusted, so spillover is held to the stricter
multiple-testing standard.

`spillover_heatmaps.R` produces the two Figure 1 panels: mean cross-module
spillover per module and timepoint, for the baseline-referenced and
increment-referenced departures.

## Regenerating from GEO

`data-raw/build_from_geo.R` rebuilds `week1_data.rds` from the raw GenomeStudio
export on GEO: normexp background correction, log2, quantile normalization, a
detection-p filter, probe-to-gene collapse, and two-stage log-scale ComBat (hyb
chamber, then cohort). 

The consensus network (`spillover_network.rds`) and responsive-gene selection
(`de_selection.rds`) are rebuilt with the packaged `build_knn_network()` and
`select_responsive_genes()`; because Infomap is randomized, a rebuilt network
carries the same partition up to label numbering. Use the provided frozen objects
to match the published module numbering exactly.

## Files

- `code/reproduce_table1.R`        end-to-end reproduction of Table 1
- `code/build_from_geo.R` regenerate week1_data.rds from GEO
- `code/figures/spillover_heatmaps.R`      Figure 1 descriptive heatmaps
- `code/figures/submodule_panels.R`      Subnetwork figures
- `code/figures/network_figure_repel.R`      Full network figure with labels
- `code/pathway_enrichment/characterize_thematic.R` enrichment labelling
- `data/module_labels_curated.csv` curated module annotations (heatmap axis labels)
- `data/de_selection.rds` limma selected temporally associated genes
- `data/spillover_network.rds` network object
- `data/table1_full_screen.rds` spillover results object
- `data/table1_reproduced` spillover results object in .csv format
- `data/week1_data.rds` cleaned data files from GEO
