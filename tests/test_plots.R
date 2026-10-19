# Run from repository root: Rscript --vanilla tests/test_plots.R
source("r/treemix_plot_functions.R")
folder <- tempfile("treemix_fixture_")
dir.create(folder)
stem <- file.path(folder,"fixture")
write_gz <- function(text,path) { con <- gzfile(path,"wt"); writeLines(text,con); close(con) }
write_gz(c("A B C", "A 0.02 0.003 0.005", "B 0.003 0.03 0.004", "C 0.005 0.004 0.04"),paste0(stem,".cov.gz"))
file.copy(paste0(stem,".cov.gz"),paste0(stem,".modelcov.gz"))
write_gz(c("A B C","A 0.001 0.001 0.001","B 0.001 0.001 0.001","C 0.001 0.001 0.001"),paste0(stem,".covse.gz"))
writeLines(c("A","B","C"),file.path(folder,"order.txt"))
pdf(file.path(folder,"residual.pdf"))
res <- plot_resid(stem,file.path(folder,"order.txt"))
dev.off()
stopifnot(all(as.matrix(res)==0), abs(get_f(stem)-1)<1e-12)
# Small tree exercises empty flip vector, zero migration edges, and original drift units.
write_gz(c("0 root ROOT NOT_MIG INTERNAL -1 1 1 2 1", "1 A NOT_ROOT NOT_MIG TIP 0 -1 0 -1 0", "2 B NOT_ROOT NOT_MIG TIP 0 -1 0 -1 0"),paste0(stem,".vertices.gz"))
write_gz(c("0 1 0.1 1 NOT_MIG 0", "0 2 0.2 1 NOT_MIG 0"),paste0(stem,".edges.gz"))
pdf(file.path(folder,"tree.pdf"))
result <- plot_tree(stem,mbar=FALSE)
dev.off()
stopifnot(abs(result$d$x[result$d[,1]==1]-.1)<1e-12,
          abs(result$d$x[result$d[,1]==2]-.2)<1e-12)
unlink(folder,recursive=TRUE)
cat("PASS: residual zero range, covariance fit, m=0 tree, empty flips, drift units\n")
