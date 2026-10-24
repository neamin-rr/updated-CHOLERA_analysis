## Recover the 95% interval for the class ARI across the 500 posterior draws.
## Run inside the Stage 3 session, right after Stage2_LPA_Louvain_final.R has run,
## while CL and cls_main are still in memory.
suppressPackageStartupMessages(library(mclust))

if (!exists("CL") || !exists("cls_main"))
  stop("CL and cls_main not found - source Stage2_LPA_Louvain_final.R first (do not clear the workspace).")

a_cls <- apply(CL, 1, function(v) adjustedRandIndex(cls_main, v))
cat(sprintf("Class ARI over %d draws: mean %.2f (95%% interval %.2f-%.2f)\n",
            length(a_cls), mean(a_cls, na.rm = TRUE),
            quantile(a_cls, 0.025, na.rm = TRUE), quantile(a_cls, 0.975, na.rm = TRUE)))

## --- permanent fix: make the main script save the interval next time -------
## In Stage2_LPA_Louvain_final.R, in the fwrite() near line 150, change
##     "ARI_class_mean"
## to
##     "ARI_class_mean", "ARI_class_lo", "ARI_class_hi"
## and in the matching value vector change
##     mean(a_cls)
## to
##     mean(a_cls), quantile(a_cls, .025), quantile(a_cls, .975)
