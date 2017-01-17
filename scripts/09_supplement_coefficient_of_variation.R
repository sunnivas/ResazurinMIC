rm(list=ls())
library("plyr")
library("prettyR")
library("ggplot2")

bestEC <- read.csv("output/tables/alldata.csv")
bestEC <- subset(bestEC,bestEC$quality=="quality ok")
bestEC <- subset(bestEC,bestEC$run=="reference")
bestEC$Estimate <- exp(bestEC$Estimate)

#cosmetic changes
bestEC <- bestEC[order(bestEC$antibiotic, bestEC$strain),]
# colnames(df) <- gsub("Std..Error", "sd", colnames(df))


## some summaries
df1 <- ddply(bestEC,~strain+antibiotic+MIC,summarise,mean=mean(Estimate),sd=sd(Estimate))

## coefficient of variation
CV <- df1$sd/df1$mean
df1$CV <- CV
# df1$dilution <-  log2(df1$mean)-log2(df1$sd)
write.table(df1,file="output/figures/coefficient_variation.txt")

##antiMIC2robial after CV
p=ggplot(df1, aes(antibiotic, CV)) +
  geom_bar(stat="identity", aes(fill=factor(antibiotic))) + facet_wrap(~strain,ncol=4) + labs(x = "", y = "CV (%)") +
  theme_bw(base_size = 14)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
p+ guides(fill=guide_legend(title="Antimicrobial"))
ggsave("output/figures/FigureS2.pdf",width=15,height=10)

