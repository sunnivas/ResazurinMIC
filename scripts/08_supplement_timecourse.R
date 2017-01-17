
rm(list=ls())
require("XLConnect")
require("ggplot2")
require("reshape")

#--------------------------------------------------------------------- Load data and merge it together
fl <- list.files("data/timecourse_data/")
fl1 <- fl[grep("^WHOF-L_Resazurin_", fl)]
fl2 <- fl[grep("^WHOM-P_Resazurin_", fl)]
fl<-c(fl1,fl2)
df <- data.frame()
for (i in unique(fl)){
  i=paste("data/timecourse_data/",i,sep="")
  for(j in getSheets(loadWorkbook(i))){
    tmp <- readWorksheetFromFile(i, sheet = j, header = T, startCol = 1, startRow = 16)
    tmp$time <- i
    tmp$strain <- j
    df <- rbind(df, tmp)
  }
}
df$time <- regmatches(df$time, regexpr("\\d+",df$time) )  ## get time point in digits
df$strain[grepl("Fluorometric1",df$strain)] <- "WHO F"   ## names strains
df$strain[grepl("Fluorometric2",df$strain)] <- "WHO G"
df$strain[grepl("Fluorometric3",df$strain)] <- "WHO K"
df$strain[grepl("Fluorometric4",df$strain)] <- "WHO L"
df$strain[grepl("FluorometricM",df$strain)] <- "WHO M"   ## names strains
df$strain[grepl("FluorometricN",df$strain)] <- "WHO N"
df$strain[grepl("FluorometricO",df$strain)] <- "WHO O"
df$strain[grepl("FluorometricP",df$strain)] <- "WHO P"

anti <- c("Spectinomycin","^H$",
          "Penicillin G","^G$",
          "Tetracycline","^F$",
          "Cefixime","^E$",
          "Ceftriaxone","^D$",
          "Ciprofloxacin","^C$",
          "Gentamicin","^B$",
          "Azithromycin","^A$")
df$Value <- gsub("^\\s+|\\s+$", "", df$Value)
for (z in seq(1,15,2)){
  df$Value[grepl(anti[z+1],df$Value)] <- anti[z]  # get antibiotics names
}
#drop Gentamicin(not included in paper because no Eucast categories available)
df <- subset(df,df$Value!="Gentamicin")
#--------------------------------------------------------------------- reshape data for ggplot
names(df)<- c("Value","Inf","100","20","4","0.8","0.16","0.032","0.0064","0.00128","0.00026","0.00005","0","time","strain")
# names(df)[!grepl("X", names(df))]
df1 <- melt.data.frame(df, id=names(df)[!grepl("X", names(df))]) 
df2 <- melt(df1)
df2$time <- as.numeric(df2$time)
names(df2) <- c("antimicrobial","time","strain","conc.","fluorescence")
ggplot(data=df2, aes(x=time, y=log(fluorescence), group=conc., colour=conc.)) +
  facet_grid(antimicrobial~strain) +
  scale_x_continuous(breaks=c(0,3,6,9,12,15), lim=c(0,16)) +
  geom_line() +
  geom_point() + theme_bw()
ggsave(file="output/figures/FigureS1.pdf", width =30, height = 20, units = "cm")
#---------------------------------------------------------------------











