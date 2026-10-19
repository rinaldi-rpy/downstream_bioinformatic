# Plot all models; model count remains a scientific decision.
out <- Sys.getenv("OUT_DIR")
summary <- read.csv(file.path(out, "model_summary.csv"))
best <- read.csv(file.path(out, "best_runs.csv"))
pdf(file.path(out, "plots", "likelihood_by_migration.pdf"), width = 8, height = 6)
plot(summary$M, summary$Best_LogLikelihood, type = "b", pch = 19,
     xlab = "Number of migration edges (m)", ylab = "Best log-likelihood",
     main = "TreeMix model fit")
dev.off()
helper <- Sys.getenv("PLOT_FUNCS")
status <- data.frame(M=best$M, Tree=FALSE, Residual=FALSE)
if (nzchar(helper) && file.exists(helper) && requireNamespace("RColorBrewer", quietly=TRUE)) {
  source(helper)
  for (i in seq_len(nrow(best))) {
    prefix <- best$BestPrefix[i]
    tree_files <- paste0(prefix, c(".vertices.gz", ".edges.gz", ".covse.gz"))
    if (all(file.exists(tree_files))) {
      pdf(file.path(out,"plots",paste0("treemix_m",best$M[i],".pdf")),width=12,height=8)
      plot_tree(prefix, mbar=best$M[i]>0)
      title(main=paste("TreeMix: migration edges =",best$M[i]))
      dev.off()
      status$Tree[i] <- TRUE
    }
    resid_files <- c(paste0(prefix,c(".cov.gz", ".modelcov.gz", ".covse.gz")), file.path(out,"population_order.txt"))
    if (all(file.exists(resid_files))) {
      pdf(file.path(out,"plots",paste0("residuals_m",best$M[i],".pdf")),width=10,height=9)
      plot_resid(prefix,file.path(out,"population_order.txt"))
      title(main=paste("Residuals: migration edges =",best$M[i]))
      dev.off()
      status$Residual[i] <- TRUE
    }
  }
} else {
  message("Tree/residual plots skipped: local treemix_plot_functions.R and RColorBrewer required.")
}
write.csv(status,file.path(out,"plot_status.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
