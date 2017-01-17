

# Description
## ----------------------------------------------------------------------------------------------------------------------------
# This file we want to create a panel with 3 plots that were saved as RDS.


# wd, clear space
## ----------------------------------------------------------------------------------------------------------------------------
rm(list = ls())


# Load dependencies
## ----------------------------------------------------------------------------------------------------------------------------
packages <- c("ggplot2", "prettyR", "corrplot", "RColorBrewer", "gridExtra", "cowplot", "grid")
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

corr <- readRDS("output/figures/FigureS3A.RDS")
caldomap <-  readRDS("output/figures/FigureS3B.RDS")

# make plot
## ----------------------------------------------------------------------------------------------------------------------------
# ggdraw() +
#   draw_plot(regr$regression, 0, 0.3, .5, .68) +
#   draw_plot(dens$Density, .5, 0.3, .5, .68) +
#   draw_plot(dev$Etestdeviation, 0, 0, 1, .28) +
#   draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.5, 0), c(1, 1, 0.3), size = 20)
# ggsave("output/figures/3_panels_figure.pdf",width=15,height=10)

ggdraw() +
  draw_plot(regr$regression, 0, 0.35, .5, .63) +
  draw_plot(dens$Density, .5, 0.35, .5, .63) +
  draw_plot(dev$Etestdeviation, 0, 0, 1, .33) +
  draw_plot_label(c("(a)", "(b)", "(c)"), c(0, 0.5, 0), c(1, 1, 0.35), size = 20)
ggsave("output/figures/Figure2.pdf",width=16,height=12.5)
# ggsave("output/figures/3_panels_figure.jpeg",width=16,height=12.5)



