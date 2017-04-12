
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file we only save the output of sessioninfo

# Load packages that were used
## ----------------------------------------------------------------------------------------------------------------------------

rm(list=ls())

packages <- c("ggplot2", "prettyR", "RColorBrewer", "gridExtra", "cowplot", 
              "grid", "plyr", "reshape", "corrplot", "sm", "heatmap3", "NMF", 
              "XLConnect", "drc", "reshape2", "qdap", "stringr", "mvtnorm", "markdown")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)

sess <- sessionInfo()
saveRDS(sess, "Session_info.rds")
