# =============================================================================
# spillover_heatmaps.R
# The two descriptive spillover heatmaps (manuscript Figure 1): mean cross-module
# spillover per module x timepoint, for the two ways of defining the propagated
# departure. Both matrices are built with the packaged functions
# (baseline_effect_vectors / increment via limma + descriptive_module_spillover);
# no reliance on the old spillover_timecourse path.
#
#   (a) baseline : each timepoint's departure from the pooled pre-infection mean
#   (b) increment: each timepoint's departure from the preceding timepoint
#
# Produces fig_spillover_baseline.pdf/.png and fig_spillover_increment.pdf/.png.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(limma); library(ggplot2) })
source("questions.R")
source("questions_prepare_terms.R")     # baseline_effect_vectors, descriptive_module_spillover

BASELINE <- c("Pre1", "Pre2"); MIN_SIZE <- 10

# ---- inputs ----
net  <- readRDS("spillover_network.rds")
sel  <- readRDS("de_selection.rds")
d    <- readRDS("week1_data.rds"); meta <- d$meta_win
expr <- as.matrix(sel$expr_de[net$genes, meta$gsm, drop = FALSE])
grab <- function(m, pats) { for (p in pats) { h <- grep(p, colnames(m), ignore.case = TRUE, value = TRUE)
  if (length(h)) return(h[1]) }; stop("column not found: ", pats[1]) }
animal <- as.character(meta[[grab(meta, c("^animal$", "monkey"))]])
tp     <- as.character(meta[[grab(meta, c("^timepoint$", "time.?point"))]])
mem <- setNames(paste0("M", net$membership), net$genes)

# curated labels for the y axis (optional; falls back to bare IDs if absent)
lab_path <- "module_labels_curated.csv"
mod_label <- if (file.exists(lab_path)) {
  L <- read.csv(lab_path, stringsAsFactors = FALSE)
  setNames(ifelse(is.na(L$label) | L$label == "", L$module, paste0(L$module, " ", L$label)), L$module)
} else NULL

# ---- increment effect vectors: each timepoint minus its predecessor ----------
increment_effect_vectors <- function(expr, group, block) {
  tp <- droplevels(factor(group)); lv <- levels(tp)
  design <- stats::model.matrix(~ 0 + tp); colnames(design) <- lv
  cc  <- limma::duplicateCorrelation(expr, design, block = factor(block))$consensus.correlation
  fit <- limma::lmFit(expr, design, block = factor(block), correlation = cc)
  cn  <- paste0(lv[-1], " - ", lv[-length(lv)])
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, limma::makeContrasts(contrasts = cn, levels = design)))
  eff <- lapply(seq_along(cn), function(j) stats::setNames(fit2$coefficients[, j], rownames(fit2)))
  stats::setNames(eff, lv[-1])
}

# ---- build both descriptive matrices via the package -------------------------
base_eff <- baseline_effect_vectors(expr, tp, block = animal, baseline_levels = BASELINE)
Mbase    <- descriptive_module_spillover(net$W_graph, mem, base_eff, min_size = MIN_SIZE)

inc_eff  <- increment_effect_vectors(expr, tp, block = animal)
Minc     <- descriptive_module_spillover(net$W_graph, mem, inc_eff, min_size = MIN_SIZE)

# ---- heatmap renderer --------------------------------------------------------
heatmap_plot <- function(M, title, time_order) {
  ord <- order(-apply(abs(M), 1, max))                 # modules by peak |spillover|
  M   <- M[ord, , drop = FALSE]
  tp_keep <- intersect(time_order, colnames(M))
  df <- expand.grid(module = rownames(M), timepoint = tp_keep, stringsAsFactors = FALSE)
  df$value <- M[cbind(df$module, df$timepoint)]
  df$module    <- factor(df$module, levels = rev(rownames(M)))
  df$timepoint <- factor(df$timepoint, levels = tp_keep)
  ylabs <- if (!is.null(mod_label)) mod_label[levels(df$module)] else levels(df$module)
  ylabs[is.na(ylabs)] <- levels(df$module)[is.na(ylabs)]
  lim <- max(abs(df$value), na.rm = TRUE)
  ggplot(df, aes(timepoint, module, fill = value)) +
    geom_tile(color = "grey92", linewidth = 0.2) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, limits = c(-lim, lim), name = "Mean\nspillover") +
    scale_y_discrete(labels = rev(ylabs)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

time_order <- c("Pre1","Pre2","D7","D10","D20","D30","D42","D56","M3")
p_base <- heatmap_plot(Mbase, "Spillover of baseline departure", time_order)
p_inc  <- heatmap_plot(Minc,  "Spillover of the timepoint increment", time_order)

ggsave("fig_spillover_baseline.pdf",  p_base, width = 5.2, height = 7)
ggsave("fig_spillover_baseline.png",  p_base, width = 5.2, height = 7, dpi = 300)
ggsave("fig_spillover_increment.pdf", p_inc,  width = 5.2, height = 7)
ggsave("fig_spillover_increment.png", p_inc,  width = 5.2, height = 7, dpi = 300)
cat("Wrote fig_spillover_baseline and fig_spillover_increment (pdf + png).\n")
