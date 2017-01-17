
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file we want plot all fitted curves, in antibiotic specific panels


# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
options(scipen = 3)

packages <- c("ggplot2", "drc", "RColorBrewer")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)

# used to format axis labels
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}

# Load data
## ----------------------------------------------------------------------------------------------------------------------------

df <- read.csv("output/tables/categories_predicted_training.csv", stringsAsFactors = F)
m2 <- readRDS("output/tables/training_data_normalized_model_list.rds")
m3 <- readRDS("output/tables/validation_data_normalized_model_list.rds")
m <- c(m2, m3)

# we only keep entries with quality ok
## ----------------------------------------------------------------------------------------------------------------------------

keeps <- df[df$quality=="quality ok"&df$run!="reference",]$ID
m <- m[keeps]

# check if names are unique
summary(duplicated(names(m2)))
summary(duplicated(names(m3)))


# get all the predictions
## ----------------------------------------------------------------------------------------------------------------------------
dp <- data.frame()

# define a concentration
conc.range <- 2^seq(-15,10,0.1)
 # conc.range <- 2^seq(-15,10,1)
l.conc.range <- log(conc.range)

for(i in 1:length(m)){
  #debug: i <- 1
  
  predi <- PR(m[[i]], l.conc.range)
  tmp <- data.frame(conc=conc.range, predi=predi, ID=names(m)[i])
  
  dp <- rbind(dp, tmp)
}

# get all the predictions
## ----------------------------------------------------------------------------------------------------------------------------
dp$ID <- as.character(dp$ID)
# df$ID <- as.character(df$ID)


dp <- merge(dp, #!is.na(df$preSIR.Etest_predicted)
            df[, c("ID","EUCAST", "antibiotic","run", "preSIR.Etest_predicted", "strain")],
            by="ID")
dp$preSIR.Etest_predicted[is.na(dp$preSIR.Etest_predicted)] <- "no EUCAST category"

dp <- transform(dp,
                preSIR.Etest_predicted = factor(preSIR.Etest_predicted , levels = c("no EUCAST category","S","I","R")))

colori <- brewer.pal(4, "Set1")[c(3,2,1,4)]
# names(dp) <- c("ID","conc","predi","EUCAST","antibiotic","run","EUCAST category","strain")

ggplot(dp, aes(x = conc, y = predi, colour = preSIR.Etest_predicted, group = ID)) + theme_bw(base_size = 10) +
  geom_line()  + scale_x_log10("Concentration (mg/L)", labels = fmt_dcimals(5)) + scale_y_continuous("Viability (in %)") +
  facet_wrap(~ antibiotic, ncol = 2, scales = "free" ) +
  scale_color_manual(values=colori) +
  theme( strip.background = element_rect(fill = "snow2"), # this is the box with I,S,R
         legend.position="top",
         plot.margin = unit(c(0.5,0.5,1.0,1.0), "cm")) 


ggsave(file="output/figures/Figure1.pdf", width =20, height = 25, units = "cm")


