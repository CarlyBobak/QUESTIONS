# =============================================================================
# network_figure_repel.R  (Version A: full labels, aggressive repel)
# Same as network_figure.R but with stronger label placement: radial OUTWARD
# nudging (push each label away from the crowded centre), higher force, a
# harder-working solver, and clip-off + expansion so labels use the margin.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(dplyr) })

LABELS_FILE <- "module_labels_curated.csv"; EMPHASIZE <- c("M6", "M8", "M91")
COLOR_TIERS <- c("high", "medium"); LABEL_TIERS <- c("high", "medium", "tentative")
RECOMPUTE_LAYOUT <- FALSE; CHARGE <- 0.001; NITER <- 800; SEED <- 1; SHOW_ISOLATED <- FALSE
EDGE_SAMPLE <- 25000; NUDGE_FRAC <- 0.15; FORCE <- 3

net <- readRDS("spillover_network.rds"); W <- net$W_graph; mem <- net$membership
if (is.null(rownames(W))) { rownames(W) <- net$genes; colnames(W) <- net$genes }
if (!RECOMPUTE_LAYOUT && file.exists("network_layout.rds")) nodes <- readRDS("network_layout.rds") else {
  suppressPackageStartupMessages(library(igraph)); A0 <- as(forceSymmetric(abs(W), "U"), "generalMatrix")
  g <- graph_from_adjacency_matrix(A0, mode = "upper", weighted = TRUE, diag = FALSE); V(g)$module <- mem[V(g)$name]
  if (!SHOW_ISOLATED) g <- delete_vertices(g, which(degree(g) == 0 | is.na(V(g)$module)))
  set.seed(SEED); lay <- layout_with_graphopt(g, charge = CHARGE, niter = NITER)
  nodes <- data.frame(gene = V(g)$name, x = lay[, 1], y = lay[, 2], module = V(g)$module); saveRDS(nodes, "network_layout.rds") }
nodes$mlab <- paste0("M", nodes$module)

lab <- read.csv(LABELS_FILE, stringsAsFactors = FALSE); lab$label[is.na(lab$label)] <- ""
tier_of <- setNames(lab$confidence, lab$module); label_of <- setNames(lab$label, lab$module)
nodes$tier <- ifelse(nodes$mlab %in% names(tier_of), tier_of[nodes$mlab], "unannotated")
nodes$label <- ifelse(nodes$mlab %in% names(label_of), label_of[nodes$mlab], "")
csub <- lab[lab$confidence %in% COLOR_TIERS & lab$module %in% nodes$mlab, ]; csub <- csub[order(-csub$size), ]
colored <- csub$module; nodes$grp <- ifelse(nodes$mlab %in% colored, nodes$mlab, "other")

A <- as(forceSymmetric(abs(W), "U"), "generalMatrix"); sm <- Matrix::summary(A); sm <- sm[sm$i < sm$j, ]; gn <- rownames(A)
edf <- data.frame(a = gn[sm$i], b = gn[sm$j], stringsAsFactors = FALSE); edf <- edf[edf$a %in% nodes$gene & edf$b %in% nodes$gene, ]
ix <- setNames(seq_len(nrow(nodes)), nodes$gene)
addxy <- function(d) { d$x <- nodes$x[ix[d$a]]; d$y <- nodes$y[ix[d$a]]; d$xend <- nodes$x[ix[d$b]]; d$yend <- nodes$y[ix[d$b]]; d }
bg <- edf; if (nrow(bg) > EDGE_SAMPLE) bg <- bg[sample(nrow(bg), EDGE_SAMPLE), ]; bg <- addxy(bg)

palette16 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
               "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#3B7EA1","#8DA0CB","#B15928")
pal <- c(setNames(palette16[seq_along(colored)], colored), other = "grey82")

cent <- nodes %>% filter(label != "" & tier %in% LABEL_TIERS) %>% group_by(mlab) %>%
  summarise(x = median(x), y = median(y), .groups = "drop")
cent$label <- label_of[cent$mlab]; cent$tier <- tier_of[cent$mlab]
cent$txt <- paste0(cent$mlab, " ", cent$label); cent$emph <- cent$mlab %in% EMPHASIZE

# radial outward nudge from plot centre
cx <- mean(nodes$x); cy <- mean(nodes$y); span <- max(diff(range(nodes$x)), diff(range(nodes$y)))
addnudge <- function(d) { dx <- d$x - cx; dy <- d$y - cy; r <- sqrt(dx^2 + dy^2); r[r == 0] <- 1
  d$nx <- dx / r * NUDGE_FRAC * span; d$ny <- dy / r * NUDGE_FRAC * span; d }
firm_ne <- addnudge(cent %>% filter(tier %in% c("high", "medium") & !emph))
firm_e  <- addnudge(cent %>% filter(tier %in% c("high", "medium") &  emph))
tent    <- addnudge(cent %>% filter(tier == "tentative"))
RP <- function(d, ...) ggrepel::geom_label_repel(data = d, aes(x, y, label = txt, ...),
  nudge_x = d$nx, nudge_y = d$ny, force = FORCE, force_pull = 0.5, max.time = 2, max.iter = 1e5,
  point.padding = unit(2, "pt"), min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf, seed = 1)

g_fig <- ggplot() +
  geom_segment(data = bg, aes(x, y, xend = xend, yend = yend), colour = "grey88", linewidth = 0.08, alpha = 0.35) +
  geom_point(data = filter(nodes, grp == "other"), aes(x, y), colour = "grey82", size = 0.4) +
  geom_point(data = filter(nodes, grp != "other" & !(grp %in% EMPHASIZE)), aes(x, y, colour = grp), size = 0.7) +
  geom_point(data = filter(nodes, grp %in% EMPHASIZE), aes(x, y, colour = grp), size = 1.3) +
  scale_colour_manual(values = pal, guide = "none") +
  RP(tent, size = I(2.3), colour = I("grey55"), fontface = I("italic"), label.size = I(0), fill = I("white"), alpha = I(0.72)) +
  RP(firm_ne, colour = mlab, size = I(2.7), label.size = I(0.12), fill = I("white"), alpha = I(0.92)) +
  RP(firm_e, colour = mlab, size = I(3.4), fontface = I("bold"), label.size = I(0.4), fill = I("white")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) + scale_y_continuous(expand = expansion(mult = 0.2)) +
  coord_equal(clip = "off") + theme_void(base_size = 11) + theme(plot.margin = margin(30, 30, 30, 30), legend.position = "none")

ggsave("fig_network_repel.pdf", g_fig, width = 12, height = 12)
ggsave("fig_network_repel.png", g_fig, width = 12, height = 12, dpi = 300)
cat("Version A (repel) written: fig_network_repel.\n")
