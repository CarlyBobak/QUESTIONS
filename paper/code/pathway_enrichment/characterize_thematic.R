# =============================================================================
# characterize_thematic.R
# Better module labels: pool the top-N terms per reference (by RANK, so the harsh
# BH penalty against big collections stops mattering), then score keywords by
# TF-IDF across modules. TF-IDF auto-suppresses ubiquitous words (regulation,
# signaling, response) and surfaces each module's DISTINCTIVE theme (platelet,
# interferon, neutrophil). BTM's discrete cell-type call kept alongside, since
# TF-IDF names functions well but BTM cleanly says T cell / monocyte.
# References: GO:BP, Reactome, WikiPathways, Hallmark (broad) + BTM (specific).
# Size-capped sets, DE-scoped background, human symbols, offline.
# =============================================================================

suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr); library(dplyr) })
stopifnot(requireNamespace("tmod", quietly = TRUE))

net <- readRDS("spillover_network.rds"); gname <- net$genes; mem <- net$membership[gname]
MIN_SIZE <- 20; MINGS <- 5; MAXGS <- 500; N_TERMS <- 25; SIG <- 0.05
mods <- as.integer(names(table(mem)[table(mem) >= MIN_SIZE]))

m <- msigdbr(species = "Homo sapiens")
catc <- intersect(c("gs_collection", "gs_cat"), colnames(m))[1]
subc <- intersect(c("gs_subcollection", "gs_subcat"), colnames(m))[1]
symc <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(m))[1]
t2g <- function(mask) data.frame(term = m[["gs_name"]][mask], gene = m[[symc]][mask])
REFS <- list(
  "GO:BP" = t2g(m[[catc]] == "C5" & m[[subc]] == "GO:BP"),
  "Reactome" = t2g(m[[catc]] == "C2" & m[[subc]] == "CP:REACTOME"),
  "WikiPathways" = t2g(m[[catc]] == "C2" & m[[subc]] == "CP:WIKIPATHWAYS"),
  "Hallmark" = t2g(m[[catc]] == "H"))
clean <- function(x) tolower(gsub("_", " ", sub("^(GOBP|GOCC|GOMF|REACTOME|WP|HALLMARK)_", "", x)))

top_terms <- function(genes, ref, n = N_TERMS) {
  e <- tryCatch(enricher(genes, TERM2GENE = ref, universe = gname, pvalueCutoff = 1, qvalueCutoff = 1,
                         minGSSize = MINGS, maxGSSize = MAXGS), error = function(x) NULL)
  if (is.null(e)) return(list(terms = character(0), best = NA))
  d <- as.data.frame(e); if (!nrow(d)) return(list(terms = character(0), best = NA))
  d <- d[order(d$p.adjust), ]; list(terms = clean(head(d$ID, n)), best = d$p.adjust[1])
}

btm_call <- function(genes, n = N_TERMS) {
  r <- tryCatch(tmod::tmodHGtest(fg = genes, bg = gname, qval = 1, mset = "LI"),
                error = function(x) tryCatch(tmod::tmodHGtest(fg = genes, bg = gname, qval = 1), error = function(y) NULL))
  if (is.null(r) || !nrow(r)) return(list(terms = character(0), label = NA, best = NA))
  r <- r[order(r$adj.P.Val), ]
  list(terms = tolower(head(r$Title, n)), label = tolower(r$Title[1]), best = r$adj.P.Val[1])
}

STOP <- c("of","the","and","to","in","via","by","response","regulation","positive","negative","process",
          "pathway","signaling","signalling","activation","system","generic","cluster","enriched","mediated",
          "involved","biological","cell","cells","for","with","from","other","related","dependent","non")
toks <- function(t) { w <- unlist(strsplit(gsub("[^a-z ]", " ", tolower(t)), "\\s+"))
  unique(w[nchar(w) > 2 & !(w %in% STOP)]) }

# ---- pool top terms per module, build token bags ----------------------------
info <- lapply(mods, function(mm) {
  genes <- names(mem)[mem == mm]
  br <- lapply(REFS, function(r) top_terms(genes, r)); bt <- btm_call(genes)
  pooled <- c(unlist(lapply(br, `[[`, "terms")), bt$terms)
  bag <- table(unlist(lapply(pooled, toks)))     # tf = # of top-terms containing token
  bests <- c(sapply(br, `[[`, "best"), BTM = bt$best)
  list(module = paste0("M", mm), size = length(genes), bag = bag,
       btm = bt$label, btm_p = bt$best, n_sig = sum(bests < SIG, na.rm = TRUE), best = min(bests, na.rm = TRUE))
})
names(info) <- sapply(info, `[[`, "module")

# ---- TF-IDF across modules: distinctive themes ------------------------------
all_tok <- unique(unlist(lapply(info, function(z) names(z$bag))))
dfreq <- sapply(all_tok, function(tk) sum(sapply(info, function(z) tk %in% names(z$bag))))
idf <- log(1 + length(info) / dfreq); names(idf) <- all_tok

theme_of <- function(z, k = 10) {
  if (!length(z$bag)) return(NA)
  s <- as.numeric(z$bag) * idf[names(z$bag)]; names(s) <- names(z$bag)
  paste(names(sort(s, decreasing = TRUE))[seq_len(min(k, length(s)))], collapse = ", ")
}

tab <- do.call(rbind, lapply(info, function(z) data.frame(
  module = z$module, size = z$size, theme = theme_of(z),
  btm_celltype = ifelse(is.na(z$btm), "-", z$btm), btm_padj = z$btm_p,
  n_sig = z$n_sig, best_padj = z$best, stringsAsFactors = FALSE))) %>% arrange(desc(size))

cat("Thematic characterization (TF-IDF distinctive keywords + BTM cell type):\n\n")
for (i in seq_len(nrow(tab))) cat(sprintf("  %-5s n=%3d sig=%d | theme: %-34s | BTM: %-22s (p=%.1g)\n",
    tab$module[i], tab$size[i], tab$n_sig[i], substr(tab$theme[i], 1, 34),
    substr(tab$btm_celltype[i], 1, 22), tab$btm_padj[i]))

write.csv(tab, "module_labels_thematic.csv", row.names = FALSE)
cat("\nWrote module_labels_thematic.csv. 'theme' = TF-IDF keywords (function); 'btm_celltype' = cell type.\n")
cat("Curate a final label per module by combining the two where they cohere.\n")
