# QUESTIONS

**Q**uantifying **E**xpression **S**pillover **T**rajectories **i**n **O**mic **N**etwork **S**eries

QUESTIONS adapts the idea of *network spillover*, developed to quantify indirect
(peer) effects under interference in social networks, to gene co-expression
networks. Treating modules as units and the network as the conduit of influence,
it separates a module's own activity (a *direct* term) from the activity reaching
it through the network (a *spillover* term), and tests whether cross-module
spillover predicts an outcome beyond the module's own activity.

The associations QUESTIONS estimates are correlational, not causal: the
underlying network is an undirected co-expression graph.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("CarlyBobak/QUESTIONS")
```

QUESTIONS depends on **limma**, which is a Bioconductor package. If you do not
already have it, install Bioconductor support first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

---

## Quick start (your own data)

You need four things:

- `expr`: a gene x sample expression matrix (rows are genes),
- `group`: a factor giving each sample's condition or timepoint,
- `unit`: a factor giving each sample's subject/unit (for repeated measures),
- `outcome`: a named numeric outcome, one value per unit.

```r
library(QUESTIONS)

# 1. (optional) scope to responsive genes and build the consensus network.
#    Skip this if you already have a weighted network + module assignment.
de  <- select_responsive_genes(expr, group = group, block = unit, fdr = 0.01)
net <- build_knn_network(de$expr_de, k = 50, sim = "signed")

# module ids come back as integers named by gene; label them once
mem <- setNames(paste0("M", net$membership), names(net$membership))   # "M8"-style

# 2. build the direct + spillover terms (long format, one row per sample x module)
long <- prepare_spillover_terms(
  network        = net$W_graph,
  membership     = mem,
  expr           = de$expr_de,
  sample_unit    = unit,
  sample_time    = group,
  outcome        = outcome,
  baseline_times = "Pre"                # level(s) of `group` that define baseline
)

# 3. define the candidate family, independent of the outcome
eff  <- baseline_effect_vectors(de$expr_de, group, block = unit, baseline_levels = "Pre")
Mmat <- descriptive_module_spillover(net$W_graph, mem, eff)   # accepts mem or integer ids
cand <- select_candidates(long, n = 15, module_spillover = Mmat, exclude_times = "Pre")

# 4. screen with small-sample validation (permutation + BH, bootstrap CI, CV, LOO)
res  <- fit_spillover_screen(long, cand, n_perm = 5000, n_boot = 5000)
res
```

`res` has one row per candidate: the standardized spillover coefficient, a
bootstrap 95% confidence interval and sign-stability, a Freedman-Lane permutation
p-value adjusted across the family (`p_BH`), the cross-validated gain in
out-of-sample R-squared from adding spillover, the direct-term p-value, and a
leave-one-out worst-case p.

### A note on module labels

`build_knn_network()` returns integer module ids named by gene. Labeling them
once as `M<k>` (as above) is optional but makes `long$module` and the rownames of
`Mmat` read consistently; both `prepare_spillover_terms()` and
`descriptive_module_spillover()` accept either the integer ids or the labels.

---

## One model at a time, and adding covariates

To fit and validate a single module x condition combination directly:

```r
fit_spillover_outcome(y = outcome_vec, direct = direct_vec, spill = spill_vec)
```

If your sample size allows, you can adjust for covariates. They enter the full
model, the reduced (permutation and CV baseline) model, the bootstrap, and CV, so
the spillover test stays a *partial* test conditional on direct + covariates:

```r
fit_spillover_outcome(y, direct, spill, covariates = Z)   # Z: matrix/data.frame, rows = units
```

---

## What the method assumes, and what to tune

QUESTIONS needs only a weighted network over units, a per-unit quantity to
propagate, and an outcome, so it applies beyond the example above. Two settings
are tuned to a small-sample design (n = 37 in the source study) and should be
revisited for your data:

- **Permutation count** (`n_perm`): the smallest achievable p is `1 / (n_perm + 1)`.
  For small outcome p-values, or many candidates, raise it (20000 or more) so the
  BH-adjusted values at the top of the ranking are stable between runs.
- **Candidate rule** (`select_candidates`): candidates should be chosen
  *independently of the outcome* so the multiple-testing family is well defined.
  The default ranks modules by peak absolute descriptive spillover; substitute a
  rule appropriate to your design if needed, but never one that looks at `outcome`.

---

## Reproducing the paper

The analysis behind the manuscript (tuberculosis disease-severity case study) is
reproduced by a standalone script and its staged data. See the `paper/` directory
(or the companion reproduction repository) and its README.

---

## Citing QUESTIONS

If you use this package, please cite the methods paper (see `inst/CITATION`), and
the social-network spillover work it adapts:

> O'Malley AJ, Bobak CA, Barnato AE. Causal Inference for Intervention Spillover
> in a Stepped-Wedge Cluster-randomized Trial: Lessons from a Physician Network.
> *Social Networks* 2026; 86:431-445. doi:10.1016/j.socnet.2026.04.017

---

## License

See `LICENSE`. (Choose and add a license before the public push.)
