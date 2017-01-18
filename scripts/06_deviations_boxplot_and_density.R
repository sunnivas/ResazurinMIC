
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# this file produces a couple of graphs:
# - boxplot with deviations from Etest
# - density of deviations from Etest


# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("plyr", "reshape2", "ggplot2", "qdap")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)


# custom functions
## ----------------------------------------------------------------------------------------------------------------------------

# function for formatting axis label decimals
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}


# load data and subset
## ----------------------------------------------------------------------------------------------------------------------------
bestEC <- read.csv("output/tables/categories_predicted.csv",stringsAsFactors = F)

bestEC$Estimate <- exp(bestEC$Estimate)
bestEC <- bestEC[bestEC$antibiotic!="Gentamicin",]
bestEC<-bestEC[bestEC$run!="reference",]

#theoretical samples 868
bestEC <- bestEC[bestEC$quality=="quality ok",]

training <- bestEC[bestEC$run == "training",]
validation <- bestEC[bestEC$run == "validation",]

# count samples with CI going over the cutoff
#----------------------------------------------------------------------------------------------
bestEC$CI_overlap <- ifelse(bestEC$CI_low < bestEC$susceptible & bestEC$CI_up > bestEC$susceptible |
                            bestEC$CI_low < bestEC$resistant & bestEC$CI_up > bestEC$resistant, 1,0  )
bestEC$correct <- grepl("^S_to_S$|^R_to_R$|^I_to_I$", bestEC$compare ) * 1

table(bestEC$CI_overlap, bestEC$correct)
summary(bestEC$CI_overlap)
write.csv(bestEC[, c("CI_overlap", "run", "antibiotic")], "output/tables/spanning_CI.csv", row.names = F )

# create cutoffs
#----------------------------------------------------------------------------------------------
bestEC$Etest.low <- as.numeric(bestEC$MIC) / 2
bestEC$Etest.high <- as.numeric(bestEC$MIC) * 2

#----------------------------------------------------------------------------------------------
bestEC1 <- cbind(training, dataset ="training")
bestEC2 <- cbind(validation, dataset ="validation")
bestEC3 <- rbind(bestEC1, bestEC2)
bestEC3$dataset <- c("validation")
#----------------------------------------------------------------------------------------------
bestEC1$deviation <- log2(as.numeric(bestEC1$Estimate))-log2(as.numeric(bestEC1$MIC))
bestEC2$deviation <- log2(as.numeric(bestEC2$Etest_predicted))-log2(as.numeric(bestEC2$MIC))
bestEC3$deviation <- log2(as.numeric(bestEC3$Etest_predicted))-log2(as.numeric(bestEC3$MIC))
bestEC$deviation <- log2(as.numeric(bestEC$Etest_predicted))-log2(as.numeric(bestEC$MIC))
#----------------------------------------------------------------------------------------------
df <- rbind(bestEC1, bestEC3)

# check if predicted is in Etest +- 1 dilutions
bestEC$prediction_in_ETest_CI <- ifelse(bestEC$Etest_predicted >= bestEC$Etest.low &
                                          bestEC$Etest_predicted <= bestEC$Etest.high, T, F  )
bestEC$overlap.lower <- ifelse(bestEC$CI_low >= bestEC$Etest.low , T, F)
bestEC$overlap.upper <- ifelse(bestEC$CI_up <= bestEC$Etest.high , T, F)
bestEC$Etest_in_prediction_CI <- ifelse(bestEC$MIC >= bestEC$CI_low &
                                        bestEC$MIC <= bestEC$CI_up, T, F)

# just take relevant information, so that it s easier to verify
# df <- bestEC[c("ETest_predicted","antibiotic", "Etest", "deviation","EC50")]
# df = df[df$run=="validation",]
## ----------------------------------------------------------------------------------------------------------------------------
BestEC.p <- bestEC
anti.sub <- c("Azithromycin","Cefixime","Ceftriaxone","Ciprofloxacin","Penicillin","Spectinomycin","Tetracycline")
anti.sub.new <- c("AZM","CFM","CRO","CIP","PEN","SPT","TET")
BestEC.p$antibiotic <- multigsub(anti.sub, anti.sub.new, BestEC.p$antibiotic)
#plot
p <- ggplot(BestEC.p, aes(x=antibiotic, y=deviation)) + geom_boxplot(aes(fill=antibiotic)) +  
  scale_y_continuous(expression("log"[2]*" deviation from Etest"),breaks=c(-10:10 ), lim=c(-10,10)) +
  scale_x_discrete("",limits = rev(anti.sub.new)) +
  geom_hline(yintercept = 0, colour="red", linetype = 3,size=1) +
  theme_bw(base_size = 18)+
  theme(legend.position="none", plot.margin = unit(c(0.5,0.3,0,1.04), "cm")) + coord_flip() 
p

p_out <-list("Etestdeviation" = p)
saveRDS(p_out, "output/figures/Etest_deviation_plot.RDS")

## ----------------------------------------------------------------------------------------------------------------------------
dev <- ggplot(df, aes(deviation, fill = dataset)) +
  geom_density(alpha = 0.2 ) +
  scale_x_continuous(expression("log"[2]*" deviation from Etest"),limits=c(-10,10),breaks=seq(-10,10,2))+
  scale_y_continuous("Density function",limits=c(0,0.3),breaks=c(0,0.05,0.10,0.15,0.20,0.25,0.3)) +
  scale_fill_manual( values = c("purple","blue"), labels=c(expression("EC"[50]*" (training)"), "Predicted MIC")) 

dev_Exp <- dev + theme_bw(base_size = 18)+
  theme(legend.justification=c(-0.01,1.01),
        legend.position=c(0,1),
        legend.text = element_text(colour="black"),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_rect(colour = "white"))+
  geom_vline(xintercept = -0.01459 ,colour="red",linetype = 3, size=1)+
  geom_vline(xintercept = -1.678,colour="red", linetype = 3,size=1 ) +
  guides(colour = guide_legend(override.aes = list(size=8)))

saveRDS(list("Density" = dev_Exp), "output/figures/Density_plot.RDS")

# calculate essential agreement
## ----------------------------------------------------------------------------------------------------------------------------

# first we need to count how many observations we have
bestEC$run<-"t"
bestEC$total_obs <- 1
cat_agr <- ddply(bestEC, ~ antibiotic + run, transform, total_obs=sum(total_obs))
# cat_agr <- ddply(df, ~ antibiotic, transform, total_obs=sum(total_obs))

# we just need to make two columns with dummies for absolute value lower than a given threshold
cat_agr$onedeviation <- ifelse( abs(cat_agr$deviation) <= 1, 1, 0)
cat_agr$twodeviation <- ifelse( abs(cat_agr$deviation) <= 2, 1, 0)
cat_agr$fourdeviation <- ifelse( abs(cat_agr$deviation) <= 4, 1, 0)
# count different agreements
cat_agr$total_cat_obs <- 1
cat_agr <- ddply(cat_agr, ~ antibiotic + run , transform, onedeviation=sum(onedeviation)/max(total_obs),
                 twodeviation=sum(twodeviation)/max(total_obs),fourdeviation=sum(fourdeviation)/max(total_obs) )

cat_agr_s <- cat_agr[,c( "run", "antibiotic", "onedeviation", "twodeviation","fourdeviation")]
cat_agr_s <- cat_agr_s[!duplicated(cat_agr_s),]
cat_agr_s$onedeviation <- round( cat_agr_s$onedeviation * 100 )
cat_agr_s$twodeviation <- round( cat_agr_s$twodeviation * 100)
cat_agr_s$fourdeviation <- round( cat_agr_s$fourdeviation * 100)
# reshape long
cat_agr_s <- melt(cat_agr_s, id.vars = c("run", "antibiotic"))

# reshape wide
cat_agr_s$total_obs <- NULL
cat_agr_s_w <- reshape(cat_agr_s, timevar = "antibiotic", idvar = c("variable","run"), direction = "wide")
cat_agr_s_w <- cat_agr_s_w[-c(1)]
names(cat_agr_s_w) <- c("Etest_deviation", "AZM","CFM","CRO","CIP","PEN","SPT","TET")
# save data
write.csv(cat_agr_s_w, "output/tables/deviations.csv", row.names = F )

# save outlier (more than four doubling dilutions)
outlier = subset(bestEC,bestEC$deviation <= -4 | bestEC$deviation>= 4)
outlier = outlier[c("ID","strain","antibiotic","MIC","Etest_predicted","deviation","compare")]
write.csv(outlier,"output/tables/outlier.csv")



