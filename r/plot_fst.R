library(ape)
out <- Sys.getenv("OUT_DIR")
mat <- as.matrix(read.csv(file.path(out, "fst_matrix.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
pairs <- read.csv(file.path(out, "pair_definitions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
results <- read.csv(file.path(out, "pairwise_fst_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
boot <- read.csv(file.path(out, "bootstrap_fst.csv"), check.names = FALSE, stringsAsFactors = FALSE)
groups <- rownames(mat)
ij <- cbind(match(pairs$Population1, groups), match(pairs$Population2, groups))
make_tree <- function(values) {
  d <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  d[ij] <- pmax(values, 0)
  d[cbind(ij[,2], ij[,1])] <- pmax(values, 0)
  nj(as.dist(d))
}
tree <- make_tree(results$FST)
trees <- lapply(seq_len(nrow(boot)), function(i) {
  make_tree(as.numeric(unlist(boot[i, as.character(pairs$Pair_ID), drop = FALSE], use.names = FALSE)))
})
class(trees) <- "multiPhylo"
support <- 100 * prop.clades(tree, trees, rooted = FALSE) / length(trees)
root <- setdiff(tree$edge[,1], tree$edge[,2])
support[root - Ntip(tree)] <- NA_real_
tree$node.label <- ifelse(is.na(support), "", sprintf("%.1f", support))
write.tree(tree, file.path(out, "trees", "fst_NJ_bootstrap.newick"))
write.tree(trees, file.path(out, "trees", "bootstrap_trees.newick"))
write.csv(data.frame(Node = Ntip(tree) + seq_len(tree$Nnode), Bootstrap_percent = support),
          file.path(out, "trees", "node_support.csv"), row.names = FALSE)
write.csv(data.frame(Metric = c("Negative_input_FST_pairs", "Negative_NJ_branch_lengths"),
                     Value = c(sum(results$FST < 0), sum(tree$edge.length < 0))),
          file.path(out, "trees", "tree_qc.csv"), row.names = FALSE)

png(file.path(out, "plots", "fst_NJ_bootstrap.png"),
    width = 12, height = 10, units = "in", res = 600)
plot(tree, type = "unrooted", use.edge.length = FALSE, show.node.label = FALSE,
     cex = 0.9, lab4ut = "horizontal", main = "FST-based Neighbor-Joining tree")
nodes <- which(!is.na(support))
nodelabels(sprintf("%.1f", support[nodes]), node = Ntip(tree)+nodes, frame = "none", cex = 0.75)
mtext(paste0("Unrooted topology; shared genomic block bootstrap (n = ", nrow(boot), "). Support: %"),
      side = 1, line = 1, cex = 0.8)
dev.off()

png(file.path(out, "plots", "pairwise_fst_heatmap.png"),
    width = 11, height = 10, units = "in", res = 600)
par(mar = c(8, 8, 5, 2))
colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
n <- nrow(mat)
span <- max(abs(mat))
if (span == 0) span <- 1
image(seq_len(n), seq_len(n), t(mat[n:1, , drop=FALSE]), axes = FALSE,
      xlab = "", ylab = "", col = colors, zlim = c(-span, span),
      main = "Pairwise Weir-Cockerham FST")
axis(1, seq_len(n), groups, las = 2, tick = FALSE, cex.axis = 0.85)
axis(2, seq_len(n), rev(groups), las = 2, tick = FALSE, cex.axis = 0.85)
for (i in seq_len(n)) for (j in seq_len(n)) {
  value <- mat[n+1-j, i]
  text(i, j, sprintf("%.3f", value), cex = 0.65,
       col = if (abs(value) > 0.65*span) "white" else "black")
}
box()
mtext("Cell labels show FST; diagonal fixed at zero", side = 3, line = 0.4, cex = 0.8)
dev.off()

ordered <- results[order(results$FST), ]
height <- max(8, 0.22*nrow(ordered)+2)
png(file.path(out, "plots", "pairwise_fst_CI95.png"),
    width = 12, height = height, units = "in", res = 600)
par(mar = c(5, 13, 4, 2))
y <- seq_len(nrow(ordered))
plot(ordered$FST, y, xlim = range(c(ordered$CI95_Lower, ordered$CI95_Upper, ordered$FST)),
     pch = 19, col = "#0072B2", yaxt = "n", ylab = "",
     xlab = "Weir-Cockerham FST (95% block-bootstrap CI)", main = "Pairwise population differentiation")
segments(ordered$CI95_Lower, y, ordered$CI95_Upper, y, col = "#0072B2")
axis(2, y, paste(ordered$Population1, "vs", ordered$Population2), las = 2, cex.axis = 0.65)
abline(v = 0, lty = 2, col = "gray60")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
