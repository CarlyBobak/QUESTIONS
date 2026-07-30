# =============================================================================
# submodule_panels.R  (zoom panels: M8, M91, M6)
# One function -> three consistent single-module gene-level panels. All genes as
# nodes, top hubs enlarged + labelled, internal edges faint with hub-hub edges
# emphasised so the within-module wiring is visible without the hairball.
# Palette matches the meta-network (M8 green, M91 purple; M6 its own rose).
# Saves each as a panel object (panel_<mod>.rds) for the composition step.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph); library(ggplot2); library(dplyr); library(ggrepel) })

HUB_N <- 10; CHARGE <- 0.015; EDGE_SAMPLE <- 1500
net <- readRDS("spillover_network.rds"); W <- net$W_graph; gname <- net$genes; mem <- net$membership[gname]
if (is.null(rownames(W))) { rownames(W) <- gname; colnames(W) <- gname }
A <- abs(W)
lab <- read.csv("module_labels_curated.csv", stringsAsFactors = FALSE); label_of <- setNames(lab$label, lab$module)
PAL <- c(M8 = "#55A868", M91 = "#8172B2", M6 = "#CC6677")

make_module_panel <- function(mm, seed = 1) {
  genes <- gname[paste0("M", mem) == mm]
  g <- graph_from_adjacency_matrix(A[genes, genes, drop = FALSE], mode = "upper", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, which(degree(g) <= 2))
  set.seed(seed); L <- layout_with_fr(g)
  nd <- data.frame(gene = V(g)$name, x = L[, 1], y = L[, 2], stringsAsFactors = FALSE)
  d <- Matrix::rowSums(A[nd$gene, nd$gene, drop = FALSE])
  hubs <- names(sort(d, decreasing = TRUE))[seq_len(min(HUB_N, length(d)))]; nd$hub <- nd$gene %in% hubs

  el <- as_edgelist(g, names = TRUE); ed <- data.frame(a = el[, 1], b = el[, 2], stringsAsFactors = FALSE)
  ed$hubedge <- ed$a %in% hubs & ed$b %in% hubs
  ix <- setNames(seq_len(nrow(nd)), nd$gene)
  ed$x <- nd$x[ix[ed$a]]; ed$y <- nd$y[ix[ed$a]]; ed$xend <- nd$x[ix[ed$b]]; ed$yend <- nd$y[ix[ed$b]]
  faint <- ed[!ed$hubedge, ]; if (nrow(faint) > EDGE_SAMPLE) faint <- faint[sample(nrow(faint), EDGE_SAMPLE), ]
  hube <- ed[ed$hubedge, ]
  col <- PAL[mm]

  p <- ggplot() +
    geom_segment(data = faint, aes(x, y, xend = xend, yend = yend), colour = col, linewidth = 0.12, alpha = 0.25) +
    geom_segment(data = hube, aes(x, y, xend = xend, yend = yend), colour = col, linewidth = 0.4, alpha = 0.65) +
    geom_point(data = filter(nd, !hub), aes(x, y), colour = col, size = 0.9, alpha = 0.55) +
    geom_point(data = filter(nd, hub), aes(x, y), colour = col, size = 2.6) +
    geom_text_repel(data = filter(nd, hub), aes(x, y, label = gene), size = 3, fontface = "bold",
      colour = "grey15", max.overlaps = Inf, seed = seed, segment.size = 0.2, min.segment.length = 0, box.padding = 0.5) +
    scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.18)) +
    coord_equal(clip = "off") + theme_void(base_size = 11) +
    ggtitle(paste0(mm, ": ", label_of[mm])) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, colour = col), plot.margin = margin(14, 14, 14, 14))

  saveRDS(p, paste0("panel_", mm, ".rds"))
  ggsave(paste0("fig_sub_", mm, ".pdf"), p, width = 6, height = 5.5)
  ggsave(paste0("fig_sub_", mm, ".png"), p, width = 6, height = 5.5, dpi = 300)
  cat(sprintf("%s: %d genes, hubs: %s\n", mm, nrow(nd), paste(hubs, collapse = ", ")))
}

for (mm in c("M8", "M91", "M6")) make_module_panel(mm, seed = 1)
