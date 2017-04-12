

# Description
## ----------------------------------------------------------------------------------------------------------------------------
# This file we want to create a panel with 3 plots that were saved as RDS.


# wd, clear space
## ----------------------------------------------------------------------------------------------------------------------------
rm(list = ls())


# Load dependencies
## ----------------------------------------------------------------------------------------------------------------------------
packages <- c("ggplot2", "prettyR", "RColorBrewer", "gridExtra", "cowplot", "grid")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)


# Load data
## ----------------------------------------------------------------------------------------------------------------------------
regr <- readRDS("output/figures/regression_plot.RDS")
dev <- readRDS("output/figures/Etest_deviation_plot.RDS")
dens <- readRDS("output/figures/Density_plot.RDS")



# make plot
## ----------------------------------------------------------------------------------------------------------------------------
# ggdraw() +
#   draw_plot(regr$regression, 0, 0.3, .5, .68) +
#   draw_plot(dens$Density, .5, 0.3, .5, .68) +
#   draw_plot(dev$Etestdeviation, 0, 0, 1, .28) +
#   draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.5, 0), c(1, 1, 0.3), size = 20)
# ggsave("output/figures/3_panels_figure.pdf",width=15,height=10)
# colorvector <- brewer.pal(9,"Set1")[c(3,2,4,1,5,9,8)]

colorvector <- brewer.pal(9,"Set1")[c(3,2,1,4,5,9,8)]
per <- paste0("(EA = ", c(53, 31, 53, 57, 61,57, 58), "%)")
names(per) <- c("AZM", "CFM", "CRO", "CIP", "PEN", "SPT", "TET")

ggdraw() +
  draw_plot(regr$regression, 0, 0.35, .5, .63) +
  draw_plot(dens$Density, .5, 0.35, .5, .63) +
  draw_plot(dev$Etestdeviation, 0, 0, 1, .33) +
  draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.5, 0), c(1, 1, 0.35), size = 20)
cowplot::ggsave("output/figures/Figure2.pdf",width=16,height=12.5)
cowplot::ggsave("output/figures/Figure2.png",width=16,height=12.5)

# ggsave("output/figures/Figure2.png",width=16,height=12.5)

# ggsave("output/figures/3_panels_figure.jpeg",width=16,height=12.5)



