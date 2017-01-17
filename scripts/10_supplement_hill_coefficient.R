
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this we make a couple of graphs with the Hill coefficient  

# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("plyr", "ggplot2",  "reshape", "corrplot", "sm", "heatmap3", "NMF")

to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)

# load data
## ----------------------------------------------------------------------------------------------------------------------------
bestEC <- read.csv("output/tables/alldata.csv",  stringsAsFactors = F)
alldata <- read.csv("output/tables/alldata.csv",  stringsAsFactors = F)

#below limit of detection: 0
alldata <- subset(alldata, alldata$quality != "below limit of detection")

#above limit of detection: training: 33 validation:
alldata <- subset(alldata, alldata$quality != "above limit of detection")


# make Hill statistic
## ----------------------------------------------------------------------------------------------------------------------------

Hill <- bestEC[bestEC$run == "validation"|bestEC$run=="training",]
Hilldata <- Hill

#test if Hill is different per antibiotic
#test if Hill is normal distributed
Hill=Hilldata$Hill
shapiro.test(Hill)

#test with pairwise t-test if difference is significant
Hillsig <- pairwise.t.test(alldata$Hill, as.factor(alldata$antibiotic))
summary(is.na(alldata$Hill))
summary(is.na(alldata$antibiotic))

# sm.density.compare(alldata$Hill, as.factor(alldata$antibiotic), xlab="Hill")
qplot(Hill, data = Hilldata, geom = "density",
      color = antibiotic)

# compute summary statistics for Hill, by antibiotic (use ddpyl so that we can get what we want)
# Hillstat <- aggregate(Hill, data=alldata,by=list(factor(alldata$antibiotic)), summary)
Hillstat <- ddply(alldata, ~ antibiotic, summarise, mean=mean(Hill),
                   sd=sd(Hill),
                   quart.25=quantile(Hill, 0.25),
                   median=median(Hill),
                   quart.75=quantile(Hill, 0.75))                                          

d.m <- melt(Hillsig$p.value)
# newnames1 <- gsub("Penicillin","Penicillin G",d.m$X1)
# newnames2 <- gsub("Penicillin","Penicillin G",d.m$X2)
# d.m$X1 <- newnames1
# d.m$X2 <- newnames2
d.m$rownames <- rownames(d.m)
d.m <- merge(d.m, Hillstat[,c("antibiotic","mean")], by.x = "X2", by.y = "antibiotic", all.x = T, sort = F )
names(d.m) <- gsub("mean$", "mean2", names(d.m) )
d.m <- merge(d.m, Hillstat[,c("antibiotic","mean")], by.x = "X1", by.y = "antibiotic", all.x = T, sort = F)
names(d.m) <- gsub("mean$", "mean1", names(d.m) )
d.m$value <- round(d.m$value, 4)

d.m$diff <- ifelse(!is.na(d.m$value), as.character(round( d.m$mean1 - d.m$mean2, 3)), NA )
d.m$diffstar <- ifelse(d.m$value <= 0.00999999, paste0(d.m$diff, "***"), d.m$diff )

d.m.wide <- reshape(d.m[,c("X1","X2","diff")],
                    timevar = "X2",
                    idvar = "X1",
                    direction = "wide")
d.m.pvalue <- reshape(d.m[,c("X1","X2","value")],
                      timevar = "X2",
                      idvar = "X1",
                      direction = "wide")
row.names(d.m.wide) <- d.m.wide$X1 
row.names(d.m.pvalue) <- d.m.wide$X1
d.m.wide$X1 <- NULL
d.m.pvalue$X1 <- NULL
names(d.m.wide) <- gsub("diff\\.", "", names(d.m.wide))
names(d.m.pvalue) <- gsub("value\\.", "", names(d.m.pvalue))
d.m.wide <- d.m.wide[, c("Azithromycin", "Cefixime", "Ceftriaxone",
                         "Ciprofloxacin", "Penicillin", "Spectinomycin")]
d.m.pvalue <- d.m.pvalue[, c("Azithromycin", "Cefixime", "Ceftriaxone",
                             "Ciprofloxacin", "Penicillin", "Spectinomycin")]

d.m.wide[is.na(d.m.wide)] <- 1
d.m.diff <- data.matrix(d.m.wide)
d.m.pmatrix <- data.matrix(d.m.pvalue)

colnames(d.m.diff) <- c("CFM","CRO","CIP","PEN","SPT","TET")
rownames(d.m.diff) <- c("AZM","CFM","CRO","CIP","PEN","SPT")
pdf("output/figures/FigureS3A.pdf",width=10,height=8)
corrplot(d.m.diff,is.corr = FALSE, method = "number",type="lower",p.mat = d.m.pmatrix, sig.level = 0.05, rect.col = "black",tl.col = "black",cl.lim=c(-2,2),cl.cex=1.5,number.cex = 1.6)
dev.off()


# analyse seperately for R and S
## ----------------------------------------------------------------------------------------------------------------------------

#R
onlyR <- subset(Hilldata,Hilldata$EUCAST=="R")
Rdistribution <- sm.density.compare(onlyR$Hill, as.factor(onlyR$antibiotic), xlab="Hill")
HillR <- pairwise.t.test(onlyR$Hill, as.factor(onlyR$antibiotic))

#S
onlyS <- subset(Hilldata,Hilldata$EUCAST=="S")
Sdistribution <- sm.density.compare(onlyS$Hill, as.factor(onlyS$antibiotic), xlab="Hill")
HillS <- pairwise.t.test(onlyS$Hill, as.factor(onlyS$antibiotic))

#Compare R and S for each antibiotic
azmR <- subset(onlyR, onlyR$antibiotic =="Azithromycin")
azmS <- subset(onlyS, onlyS$antibiotic =="Azithromycin")
azmT <- t.test(azmR$Hill,azmS$Hill)

cipR <- subset(onlyR, onlyR$antibiotic =="Ciprofloxacin")
cipS <- subset(onlyS, onlyS$antibiotic =="Ciprofloxacin")
cipT <- t.test(cipR$Hill,cipS$Hill)

cxR <- subset(onlyR, onlyR$antibiotic =="Cefixime")
cxS <- subset(onlyS, onlyS$antibiotic =="Cefixime")
cxT <- t.test(cxR$Hill,cxS$Hill)

penR <- subset(onlyR, onlyR$antibiotic =="Penicillin")
penS <- subset(onlyS, onlyS$antibiotic =="Penicillin")
penT <- t.test(penR$Hill,penS$Hill)

EC50 <- Hilldata
AZM <- EC50[EC50$antibiotic=="Azithromycin",]$Hill
CIP <- EC50[EC50$antibiotic=="Ciprofloxacin",]$Hill
CT <- EC50[EC50$antibiotic=="Ceftriaxone",]$Hill
CX <- EC50[EC50$antibiotic=="Cefixime",]$Hill
TET <- EC50[EC50$antibiotic=="Tetracycline",]$Hill
PEN <- EC50[EC50$antibiotic=="Penicillin",]$Hill
SPEC <- EC50[EC50$antibiotic=="Spectinomycin",]$Hill

Hill.m=matrix(data=cbind(AZM,CIP,CT,CX,TET,PEN,SPEC),ncol=7)
colnames(Hill.m) <- c("AZM","CIP","CRO","CFM","TET","PEN","SPT")
rownames(Hill.m) <- c(1:124)


pdf("output/figures/FigureS3B.pdf")
aheatmap(Hill.m)
dev.off()

# try to plot both
## -----------------------------------------------------------------------------------------------------------
#par(mfrow=c(1,2))
pdf("output/figures/FigureS3AB.pdf", width = 20, height = 10)
#par(mfrow=c(1,2))
layout(matrix(c(1,2), 1, 2, byrow = TRUE),
       widths=c(4.5,3))
corrplot(d.m.diff,is.corr = FALSE, method = "number",type="lower",p.mat = d.m.pmatrix, sig.level = 0.05,
         rect.col = "black",tl.col = "black",cl.lim=c(-2,2),cl.cex=1.5,number.cex = 1.6)
aheatmap(Hill.m,labRow=NA)
dev.off()

# file end
## -----------------------------------------------------------------------------------------------------------
