# =============================================================================
# data-raw/build_from_geo.R
# Regenerate the staged reproduction input from GEO accession GSE84152.
#
# Produces week1_data.rds with exactly what the current pipeline consumes:
#     meta_win      sample metadata (gsm, animal, timepoint, fdg, fdg_group, ...)
#     expr_all_win  ComBat-corrected log2 gene x sample matrix over the analysis window
#     animals_no_fdg animals lacking an FDG group (carried for bookkeeping)
#     params        the settings used, for provenance
#
# Downstream (select_responsive_genes -> build_knn_network -> the spillover
# pipeline) derives everything else; the earlier WGCNA/variable-gene outputs are
# not part of the current method and are not produced here.
#
# The series matrix stored linear, mean-scaled values with a linear-scale ComBat
# already applied, which is the wrong scale for co-expression. We reprocess from
# the raw GenomeStudio export: normexp background correction, log2, quantile
# normalization, detection-p filter, probe->gene collapse, then log-scale ComBat
# for the technical batches (hyb chamber, then cohort).
#
# TWO STEPS TO VERIFY (they cannot be fully checked automatically):
#   [VERIFY 1] the raw-array-header -> GSM mapping (spot-check is printed)
#   [VERIFY 2] the ComBat batch assignment (chamber/cohort must be complete)
#
# USAGE:  Rscript data-raw/build_from_geo.R   (writes week1_data.rds)
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(WGCNA); library(sva)
  library(data.table); library(dplyr); library(tidyr); library(stringr); library(tibble)
})
options(stringsAsFactors = FALSE)

PARAMS <- list(
  gse         = "GSE84152",
  cache_dir   = "data_cache",
  window      = c("Pre1","Pre2","D7","D10","D20","D30","D42","D56","M3"),  # D3 dropped
  detection_p = 0.01,          # keep probes detected (p < this) in >= 1 sample
  apply_combat = TRUE,         # correct hyb chamber + cohort on the log scale
  out_file    = "week1_data.rds"
)
dir.create(PARAMS$cache_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Series matrix (metadata + annotation) and the raw supplementary file
# -----------------------------------------------------------------------------
rds <- file.path(PARAMS$cache_dir, "gse_raw.rds")
gse <- if (file.exists(rds)) readRDS(rds) else {
  g <- getGEO(PARAMS$gse, GSEMatrix = TRUE, destdir = PARAMS$cache_dir)[[1]]; saveRDS(g, rds); g
}
supp_dir <- file.path(PARAMS$cache_dir, PARAMS$gse)
raw_path <- list.files(supp_dir, pattern = "non-normalized", full.names = TRUE)
if (length(raw_path) == 0) {
  getGEOSuppFiles(PARAMS$gse, baseDir = PARAMS$cache_dir)
  raw_path <- list.files(supp_dir, pattern = "non-normalized", full.names = TRUE)
}
cat("Raw file:", raw_path, "\n")

# -----------------------------------------------------------------------------
# 2. Metadata, with technical-batch factors for ComBat
# -----------------------------------------------------------------------------
pd <- pData(gse)
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
meta <- pd %>% rownames_to_column("gsm") %>% select(gsm, all_of(char_cols)) %>%
  pivot_longer(-gsm, values_to = "kv") %>%
  filter(kv != "", str_detect(kv, ": ")) %>%
  separate(kv, into = c("key","val"), sep = ": ", extra = "merge") %>%
  distinct(gsm, key, .keep_all = TRUE) %>% select(gsm, key, val) %>%
  pivot_wider(names_from = key, values_from = val) %>%
  mutate(animal = monkeyid,
         timepoint = factor(`infection time`,
           levels = c("Pre1","Pre2","D3","D7","D10","D20","D30","D42","D56","M3","M4","M5","M6")),
         clinical = `clinical status`, fdg_group = `FDG group`,
         fdg = suppressWarnings(as.numeric(`FDG avidity`)),
         chamber = `hyb chamber`, cohort = dataset)

modal <- function(x){ x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else names(sort(table(x), decreasing = TRUE))[1] }
lab <- meta %>% group_by(animal) %>%
  summarise(clinical = modal(clinical), fdg_group = modal(fdg_group), cohort = modal(cohort), .groups = "drop")
stopifnot(nrow(lab) == 38)
meta <- meta %>% select(-clinical, -fdg_group, -cohort) %>% left_join(lab, by = "animal")
animals_no_fdg <- lab %>% filter(is.na(fdg_group)) %>% pull(animal)

# -----------------------------------------------------------------------------
# 3. Raw signal + detection matrices
# -----------------------------------------------------------------------------
dt <- data.table::fread(raw_path)                          # reads .gz directly
ids <- dt[["ID_REF"]]
sig_cols <- grep("\\.AVG_Signal$", names(dt), value = TRUE)
det_cols <- grep("\\.Detection Pval$", names(dt), value = TRUE)
sids <- sub("\\.AVG_Signal$", "", sig_cols)
dids <- sub("\\.Detection Pval$", "", det_cols)
det_cols <- det_cols[match(sids, dids)]                    # align detection to signal order

E   <- as.matrix(dt[, ..sig_cols]); rownames(E) <- ids; colnames(E) <- sids
Det <- as.matrix(dt[, ..det_cols]); rownames(Det) <- ids; colnames(Det) <- sids
cat(sprintf("Raw: %d probes x %d arrays | signal %.1f to %.1f\n",
            nrow(E), ncol(E), min(E, na.rm = TRUE), max(E, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. normexp background correction -> log2 -> quantile normalize; detection filter
# -----------------------------------------------------------------------------
Ebc  <- backgroundCorrect(E, method = "normexp", offset = 16, verbose = FALSE)
En   <- normalizeBetweenArrays(log2(Ebc), method = "quantile")
cat(sprintf("Normalized log2 range: %.2f to %.2f (expect ~4-14)\n",
            min(En, na.rm = TRUE), max(En, na.rm = TRUE)))
keep <- rowSums(Det < PARAMS$detection_p, na.rm = TRUE) >= 1
En <- En[keep, , drop = FALSE]
cat(sprintf("Probes passing detection filter: %d\n", nrow(En)))

# -----------------------------------------------------------------------------
# 5. Map raw-array columns -> GSM.  [VERIFY 1]
# GEO stored each sample's verbatim raw-file column header in pData
# description.1, so we match raw signal headers against it and rename to GSM.
# -----------------------------------------------------------------------------
map_col <- "description.1"; stopifnot(map_col %in% colnames(pd))
gsm_for_raw <- rownames(pd)[match(sig_cols, pd[[map_col]])]
n_mapped <- sum(!is.na(gsm_for_raw))
cat(sprintf("Mapped %d / %d arrays to a GSM via '%s'\n", n_mapped, length(sig_cols), map_col))
stopifnot(n_mapped >= length(sig_cols) * 0.95)
colnames(En) <- gsm_for_raw
En <- En[, !is.na(colnames(En)), drop = FALSE]

cat("\n[VERIFY 1] array header -> gsm -> label spot-check (confirm these are sensible):\n")
print(data.frame(array = sig_cols[1:6], gsm = gsm_for_raw[1:6],
                 label = pd$description[match(gsm_for_raw[1:6], rownames(pd))]))

# -----------------------------------------------------------------------------
# 6. Collapse probes to genes (MaxMean)
# -----------------------------------------------------------------------------
fdat <- fData(gse)
sym_col <- c("Symbol","SYMBOL","ILMN_Gene","Gene Symbol")[c("Symbol","SYMBOL","ILMN_Gene","Gene Symbol") %in% colnames(fdat)][1]
if (is.na(sym_col)) { print(colnames(fdat)); stop("Set sym_col from the fData columns above.") }
sym <- fdat[rownames(En), sym_col]; ok <- !is.na(sym) & sym != ""
En <- En[ok, , drop = FALSE]; sym <- sym[ok]
expr_gene_all <- WGCNA::collapseRows(En, rowGroup = sym, rowID = rownames(En), method = "MaxMean")$datETcollapsed
cat(sprintf("Collapsed to %d genes x %d samples\n", nrow(expr_gene_all), ncol(expr_gene_all)))

# -----------------------------------------------------------------------------
# 7. ComBat on the log scale: hyb chamber, then cohort.  [VERIFY 2]
# -----------------------------------------------------------------------------
if (PARAMS$apply_combat) {
  m <- meta[match(colnames(expr_gene_all), meta$gsm), ]
  cham <- factor(m$chamber); coh <- factor(m$cohort)
  if (anyNA(cham) || anyNA(coh)) stop("[VERIFY 2] missing chamber/cohort for some samples; cannot ComBat.")
  expr_gene_all <- sva::ComBat(expr_gene_all, batch = cham, mod = NULL, par.prior = TRUE)
  expr_gene_all <- sva::ComBat(expr_gene_all, batch = coh,  mod = NULL, par.prior = TRUE)
  cat("Applied ComBat for hyb chamber and cohort.\n")
}

# -----------------------------------------------------------------------------
# 8. Restrict to the analysis window and save the lean object
# -----------------------------------------------------------------------------
meta_win <- meta %>% filter(gsm %in% colnames(expr_gene_all), timepoint %in% PARAMS$window) %>% droplevels()
cat("\nSamples per timepoint (window):\n"); print(table(meta_win$timepoint))
expr_all_win <- expr_gene_all[, meta_win$gsm, drop = FALSE]
cat(sprintf("Windowed matrix: %d genes x %d samples\n", nrow(expr_all_win), ncol(expr_all_win)))

saveRDS(list(meta_win = meta_win, expr_all_win = expr_all_win,
             animals_no_fdg = animals_no_fdg, params = PARAMS), PARAMS$out_file)
cat(sprintf("\nSaved %s (meta_win, expr_all_win, animals_no_fdg, params)\n", PARAMS$out_file))
