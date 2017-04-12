
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file we fit a linear regression to the training data and plot the model

# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("ggplot2", "RColorBrewer")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)


# load and subset data (to calculate the number of evaluable strains)
## ----------------------------------------------------------------------------------------------------------------------------
alldata <- read.csv("output/tables/alldata.csv",  stringsAsFactors = F)

#above limit of detection (resazurin method): 
alldata <- subset(alldata, alldata$quality != "above limit of detection")

#above or below limit of detection (Etest method): 
alldata$above <- grepl(">",alldata$MIC)
length(alldata[alldata$above==T,]$above)

alldata$below <- grepl("<",alldata$MIC)
length(alldata[alldata$below==T,]$below)

cleandata <- subset(alldata, alldata$above!=TRUE & alldata$below!=TRUE)
# cleandata_t <- subset(cleandata, cleandata$run == "training")
# cleandata_v <- subset(cleandata, cleandata$run == "validation")
# cleandata_r <- subset(cleandata, cleandata$run == "reference")

training <- subset(cleandata,cleandata$run == "training")
validation <- subset(cleandata,cleandata$run == "validation")


# estimate linear regression model
##-----------------------------------------------------------------------------------------------------------------------
etest <- as.numeric(training$MIC)
esti <-  as.numeric(training$Estimate)
print(summary(linreg <- lm(log(etest)~(esti))))
esti_fit_log <- linreg$fitted.values
esti_fit <- exp(linreg$fitted.values)
prediction <- data.frame(esti_fit)

# save parameters and parameters variance~covariance matrix in a named list, we will use it to bootstrap
vcov.list <- list(list("Estimates"=linreg, "Matrix"=list(vcov(linreg), "Summary"=summary(linreg))))
saveRDS(vcov.list, "output/tables/lm_parameters_variance_covariance_matrix_list.rds")

# save pearson correlation
pearson <- cor(log(etest),esti,method="pearson")
saveRDS(pearson, "output/tables/pearson.rds")

# Plot regression and save it as rds object, so that we can combine it with plots from
# other scripts
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# get fitted values
df4 <- as.data.frame(cbind(esti_fit,training))

# set plot options
par(mar = c(5,5,5,5))
par(oma = c(1,1,1,1))
options(scipen = 100)

# function for formatting axis label decimals
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}
# colorvector <- c( "#D95F02", "#7570B3", "#E7298A" ,"#66A61E", "#1F78B4", "#A6761D", "#666666")
colorvector <- brewer.pal(9,"Set1")[c(3,2,1,4,5,9,8)]
#colorvector <- brewer.pal(8,"Set2")

dfplot <- df4
dfplot$antibiotic <- gsub("Penicillin","Penicillin", dfplot$antibiotic ) # cosmetic changes

p1<-ggplot(dfplot, aes (x = exp(Estimate), y = as.numeric(MIC))) +
  #geom_text(data = data.frame(), aes(0.05, 0.0005, label = ""),cex=16, fontface="bold")+
  geom_point(size=2,alpha=0.6,aes(colour = factor(antibiotic)))+ #labs(aes(colour=antibiotic))+
  geom_smooth(method = "lm",color="steelblue2",show.legend = F,se = F,level=0.95)+
  scale_y_continuous("Etest MIC (mg/L)",trans='log10' , limits = c(0.00001, 5000),
                     breaks=c(0.00001,0.0001,0.001, 0.01,0.1,1,10,100,1000), labels = fmt_dcimals(5))+
  scale_x_continuous(expression("EC"[50]*" (mg/L)"), trans='log10' , limits = c(0.00001, 5000),
                     breaks=c(0.00001,0.0001,0.001, 0.01,0.1,1,10,100,1000), labels = fmt_dcimals(5)) +
  scale_color_manual(values=colorvector) + geom_abline(intercept = 0.0000001, slope = 1, color="black", lty=2)+
  theme_bw(18)+theme(legend.justification=c(-0.01,1.01),
                   legend.position=c(0,1),
                   legend.text = element_text(colour="black"),
                   legend.title = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.key = element_rect(colour = "white")) +
                   guides(colour = guide_legend(override.aes = list(size=8))) 
p1

regr_plot <- list("regression" = p1)
saveRDS(regr_plot, "output/figures/regression_plot.RDS")


# file end
## -----------------------------------------------------------------------------------------------------------




