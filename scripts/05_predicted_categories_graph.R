
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# this file makes a graph with the categories predicted by the model.
# If you want to go back to having S,I,R instead of S,R search "changeItoS" in this file, and comment out
# those 4 lines of code. The name used to save the plots depends on whehter we have S,I,R or S,R.


# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("plyr", "stringr", "reshape2", "ggplot2", "gridExtra", "grid", "RColorBrewer")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)
options(scipen = 3)

# function for formatting axis label decimals
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}


# Load data and subset
## ----------------------------------------------------------------------------------------------------------------------------
bestEC <- read.csv("output/tables/Predicted_Etest_all_data.csv", stringsAsFactors = F)
out1 <- bestEC  # we need this to save the data frame at the end
bestEC <- bestEC[!bestEC$run=="reference" ,]
bestEC <- bestEC[!bestEC$antibiotic=="Gentamicin" ,]
bestEC <- subset(bestEC, grepl("^quality ok$|^Error calculation fails$|^above limit of detection$|^Etest above limit of detection$|^Etest below limit of detection$",
                               bestEC$quality))

# Confidence interval omission when too big
## ----------------------------------------------------------------------------------------------------------------------------

dci_drops <- bestEC$CI_up/bestEC$Etest_predicted
bestEC$CIrm <- dci_drops > quantile(dci_drops, 0.99, na.rm = T)
bestEC$CIrm[grepl("limit of detection",bestEC$quality)] <- T

bestEC$CI_low[bestEC$CIrm==T] <- bestEC$Etest_predicted[bestEC$CIrm==T] 
bestEC$CI_up[bestEC$CIrm==T] <- bestEC$Etest_predicted[bestEC$CIrm==T] 


# add information about cutoffs for susceptible and resistant (note that sometimes they are the same, in those cases there is no I)
## ----------------------------------------------------------------------------------------------------------------------------
bestEC$susceptible[bestEC$antibiotic=="Azithromycin"] <- 0.25
bestEC$susceptible[bestEC$antibiotic=="Cefixime"] <- 0.125
bestEC$susceptible[bestEC$antibiotic=="Ceftriaxone"] <- 0.125
bestEC$susceptible[bestEC$antibiotic=="Ciprofloxacin"] <- 0.03
bestEC$susceptible[bestEC$antibiotic=="Gentamicin"] <- NA
bestEC$susceptible[bestEC$antibiotic=="Penicillin"] <- 0.06
bestEC$susceptible[bestEC$antibiotic=="Spectinomycin"] <- 64
bestEC$susceptible[bestEC$antibiotic=="Tetracycline"] <- 0.5

bestEC$resistant[bestEC$antibiotic=="Azithromycin"] <- 0.5
bestEC$resistant[bestEC$antibiotic=="Cefixime"] <- 0.125
bestEC$resistant[bestEC$antibiotic=="Ceftriaxone"] <- 0.125
bestEC$resistant[bestEC$antibiotic=="Ciprofloxacin"] <- 0.06
bestEC$resistant[bestEC$antibiotic=="Gentamicin"] <- NA
bestEC$resistant[bestEC$antibiotic=="Penicillin"] <- 1
bestEC$resistant[bestEC$antibiotic=="Spectinomycin"] <- 64
bestEC$resistant[bestEC$antibiotic=="Tetracycline"] <- 1


# now we need to check whether the predicted value falls in the right confidence interval. We also do this using the uppuer and 
# lower estimates for the predicted value.
## ----------------------------------------------------------------------------------------------------------------------------

# first keep only columns of interest, we will merge the rest back later
d <- bestEC[!is.na(bestEC$EUCAST) ,c("ID","antibiotic","EUCAST","Etest_predicted", "CI_low", "CI_up", "susceptible", "resistant","MIC", "CIrm" )]

# to make computations easier we melt long so that estimate and CI are in the same column
to_be_melted <- c("Etest_predicted", "CI_low", "CI_up")
id.vars <- setdiff(names(d), to_be_melted)
d.m <- melt(d, id.vars = id.vars  )
d.m$preSIR <- ifelse(d.m$value <= d.m$susceptible, "S",
                              ifelse(d.m$value >= d.m$susceptible &
                                     d.m$value <= d.m$resistant &
                                     d.m$susceptible < d.m$resistant, "I", "R" ))

# now we need to reshape wide
d.w <- reshape(d.m, idvar = id.vars, timevar = "variable", direction = "wide" )
names(d.w) <- gsub("value\\.", "", names(d.w))


# in case we want to reduce categories to two only
## ----------------------------------------------------------------------------------------------------------------------------
d.w$resistant_g <- d.w$resistant

# changeItoS
# comment out the next 4 rows of code to get the graphs in S,I,R 3x3 matrix
# d.w$EUCAST[d.w$EUCAST=="I"] <- "S"
# d.w$preSIR.Etest_predicted[d.w$preSIR.Etest_predicted=="I"] <- "S"
# dropI <- d.w$susceptible!=d.w$resistant
# d.w$resistant[dropI] <- d.w$susceptible[dropI]


# count data
## ----------------------------------------------------------------------------------------------------------------------------

# count total obs per antibiotic + run + EUCAST and antibiotic + run
d.w$total_obs <- 1
d.w <- ddply(d.w, ~ antibiotic + EUCAST, transform, total_obs=sum(total_obs) )
d.w$anti_total_obs <- 1
d.w <- ddply(d.w, ~ antibiotic  , transform, anti_total_obs=sum(anti_total_obs) )

# count the cross classifications
d.w$compare <- paste0(d.w$EUCAST,"_to_",d.w$preSIR.Etest_predicted)

# add matrices with confidence interval classification
tmp <- data.frame(model.matrix( ~ preSIR.CI_low - 1 , data = d.w))  
d.w <- cbind(d.w , tmp)
tmp <- data.frame(model.matrix( ~ preSIR.CI_up - 1 , data = d.w))  
d.w <- cbind(d.w , tmp)

# compute the percentages for the paper
perce_overall <- ddply(d.w, ~ compare, summarise ,
                       n = length(compare), n_per = round(length(compare)/nrow(d.w),2) )
perce_overall_anti <- ddply(d.w, ~ compare + antibiotic, summarise ,
                            n = length(compare),n_per = round(length(compare)/nrow(d.w),3) )

#perce_overall_anti <- ddply(d.w, ~ antibiotic, transform , strain_count = length(antibiotic) )
# S_S <-length(d.w[d.w$compare=="S_to_S",]$compare)
# S_R <-length(d.w[d.w$compare=="S_to_R",]$compare)
# R_S <-length(d.w[d.w$compare=="R_to_S",]$compare)
# R_R <-length(d.w[d.w$compare=="R_to_R",]$compare)
# I_S <-length(d.w[d.w$compare=="I_to_S",]$compare)
# S_I <-length(d.w[d.w$compare=="S_to_I",]$compare)
# I_I <-length(d.w[d.w$compare=="I_to_I",]$compare)
# I_R <-length(d.w[d.w$compare=="I_to_R",]$compare)
# R_I <-length(d.w[d.w$compare=="R_to_I",]$compare)

# try to plot
## ----------------------------------------------------------------------------------------------------------------------------

for(anti in unique(d.w$antibiotic)){
  
  # debug: anti <- "Azithromycin"
  # debug: anti <- "Cefixime"
  # debug: anti <- "Ceftriaxone"
  # debug: anti <- "Ciprofloxacin"
  # debug: anti <- "Penicillin"
  # debug: anti <- "Spectinomycin"
  # debug: anti <- "Tetracycline"
  
  print(anti)
  tmp1 <- d.w[d.w$antibiotic==anti,c("antibiotic","EUCAST","preSIR.Etest_predicted","Etest_predicted",
                                     "CI_up", "CI_low", "resistant", "susceptible", "resistant_g", "CIrm")]

  #for(class in unique(tmp1$EUCAST)){
    # debug: class <- "S"
  
  #tmp <- tmp1[tmp1$EUCAST==class, ]
  tmp <- tmp1
  tmp <- tmp[order(tmp$Etest_predicted), ] # first sort
  
  # prepare string for title and output
  title <- unique(tmp$antibiotic)
  if(exists("dropI")){output.ext <- "_onlySR"}else{output.ext<-""} # this is added only if we decided to keep  S and R only
  title_out <- paste0("output/figures/", title,output.ext, ".pdf")

  tmp$order <- 1
  tmp <- ddply(tmp, ~ antibiotic + EUCAST + preSIR.Etest_predicted, transform, order=cumsum(order))

  # get categorical percentage for EUCAST
  tmp <- ddply(tmp, ~ antibiotic + EUCAST , transform, count=length(EUCAST))
  tmp$EUCASTn <- paste0(tmp$EUCAST, " (n = ", tmp$count, ")" )
  tmp <- ddply(tmp, ~ antibiotic + preSIR.Etest_predicted , transform, count=length(preSIR.Etest_predicted))
  tmp$preSIR.Etest_predictedn <- paste0(tmp$preSIR.Etest_predicted, " (n = ", tmp$count, ")" )
  
  # get the percentage to by displayed in each plot, like in a contingency table
  tmp$percento <- 1
  tmp <- ddply(tmp, ~ antibiotic + EUCAST + preSIR.Etest_predicted , transform, percento=sum(percento)/nrow(tmp))
  tmp$percento <- round(tmp$percento, digits = 2) * 100
  tmp <- ddply(tmp, ~ antibiotic + EUCAST + preSIR.Etest_predicted , transform, count=length(EUCAST))
  #tmp$label <- paste0("n=",tmp$count, "\n(", tmp$percento,"%)")
  tmp$label <- paste0("n=",tmp$count)
  
  
  # create a variable to spread out the observations on the y axis
  tmp$order <- round(tmp$order/tmp$count, digits = 4) * 100
  
  # create a variable for coloring of points
  tmp$class_d <- paste0(tmp$EUCAST, "_to_", tmp$preSIR.Etest_predicted)
  tmp$color <- "black"
  tmp$color[grepl("I_to_I|R_to_R|S_to_S",tmp$class_d)] <- "green"
  tmp$color[grepl("R_to_S|S_to_R",tmp$class_d)] <- "red"
  
  # we need to add a fake observation (and remove afterwards from plot) for Ciprofloxacin because there are no EUCAST I 
  if(tmp$antibiotic[1]=="Ciprofloxacin"){

    tmp.cip <- tmp[1,]
    tmp.cip$EUCAST <- "I"
    tmp.cip$preSIR.Etest_predicted <- "I"
    tmp.cip$Etest_predicted <- 0.0463
    tmp.cip$CI_up <- 1.1
    tmp.cip$CI_low <- 0.001
    tmp.cip$EUCASTn <- "I (n = 0)"
    tmp.cip$preSIR.Etest_predictedn <- "I (n = 1)"
    tmp.cip$label <- "n=0"
    tmp <- rbind(tmp, tmp.cip)
  }
  
  # need to order factors
  EU_levels <- NULL
  EU_levels <- unique(tmp$EUCASTn)
  EU_levels <- c(EU_levels[grep("^S", EU_levels)],
                 EU_levels[grep("^I", EU_levels)],
                 EU_levels[grep("^R", EU_levels)])
  Mo_levels <- unique(tmp$preSIR.Etest_predictedn)
  Mo_levels <- c(Mo_levels[grep("^S", Mo_levels)],
                 Mo_levels[grep("^I", Mo_levels)],
                 Mo_levels[grep("^R", Mo_levels)])
  
  tmp <- transform(tmp,
                   EUCASTn = factor(EUCASTn, levels = EU_levels),
                   preSIR.Etest_predictedn = factor(preSIR.Etest_predictedn , levels = Mo_levels),
                   class_percentage_f = runif(nrow(tmp)))
  
  # need to set colors manually, as we want to keep colors fixed
  #display.brewer.pal(3, "Set1")
  colori <- NULL
  if (length(EU_levels)==3 ) {
    colori <- brewer.pal(3, "Set1")[c(2,3,1)]
  } else {
    colori <- brewer.pal(3, "Set1")[c(3,1)]
  }
  print(colori)
  
  # cutoffs
  cuts <- unique(c(unique(tmp$resistant_g), unique(tmp$susceptible) ))

  # need to compute numbers for annotation of number of observation
  addobs <- tmp[,c("EUCASTn", "preSIR.Etest_predictedn","label")]
  addobs <- addobs[!duplicated(addobs),]

  addobs <- merge(addobs,
                  expand.grid(EUCASTn=unique(addobs$EUCASTn), preSIR.Etest_predictedn=unique(addobs$preSIR.Etest_predictedn)),
                  by=c("EUCASTn", "preSIR.Etest_predictedn"), all = T)
  addobs$label[is.na(addobs$label)] <- "n=0"
  n_annotate_y <- max(c(tmp$CI_u, tmp$Etest_predicted))
  n_annotate_x <- 0
  print(length(tmp$percento))
  # make plot
grafico <-   ggplot(tmp, aes(order, Etest_predicted)) + 
    geom_errorbar(aes(ymax = CI_up, ymin = CI_low, color=color,  width=0.0), size=0.5 ) + #
    geom_point(aes( shape=CIrm), size=0.5) + ggtitle(title) +
    coord_flip() + scale_y_log10("",  labels = fmt_dcimals(3)) + scale_x_continuous("", limits = c(-1,100)) + 
    facet_grid(EUCASTn ~ preSIR.Etest_predictedn, switch = "both", margins = F ) +
    geom_hline(yintercept=cuts, linetype="dotted", size=0.28) +
    scale_color_manual(values=colori) +
    theme_bw(base_size = 10 ) +
    geom_text(data=addobs, aes(x=n_annotate_x, y=n_annotate_y, label=label),
            colour="black", inherit.aes=FALSE, parse=F, size=3.5, hjust=1.0, vjust=0) +
    theme(axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         # axis.title.y = element_text("bla",size = rel(1.5), angle = 90),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        strip.background = element_rect(fill = "snow2"), # this is the box with I,S,R
         legend.position="bottomright",
        plot.margin = unit(c(0.1,0.105,-0.0,-0.2), "cm"))
grafico
  
  # save
  assign(title, grafico)
  ggsave(paste0("output/figures/single_", title, ".pdf"), grafico,  width = 10, height = 10, units = "cm" )

  }

# we want to save all pictures in a panel (2 times, one for the paper and one for summary stats html)
panel <- grid.arrange(arrangeGrob(
             Azithromycin,
             Ciprofloxacin,
             Ceftriaxone,
             Cefixime,
             Penicillin,
             Spectinomycin,
             Tetracycline,
             ncol=4,
             bottom = textGrob("Model classification", vjust = 0.3),
             left = textGrob("EUCAST classification", rot = 90, vjust = 1)))

title_out <- paste0("output/figures/Figure3",output.ext, ".pdf")
ggsave(file=title_out, panel, width =35, height = 20, units = "cm")
title_outp <- paste0("output/figures/Figure3",output.ext, ".png")
ggsave(file=title_outp, panel, width =35, height = 20, units = "cm")


panel <- grid.arrange(arrangeGrob(
  Azithromycin,
  Ciprofloxacin,
  Ceftriaxone,
  Cefixime,
  Penicillin,
  Spectinomycin,
  Tetracycline,
  ncol=2,
  bottom = textGrob("Model classification", vjust = 0.3),
  left = textGrob("EUCAST classification", rot = 90, vjust = 1)))

# title_out <- paste0("output/figures/Figure3_2col",output.ext, ".png")
# ggsave(file=title_out, panel, width =30, height = 60, units = "cm")
title_out <- paste0("output/figures/Figure3_2col",output.ext, ".pdf")
ggsave(file=title_out, panel, width =30, height = 60, units = "cm")
# we also want to save the data, will be used in script 9
## ----------------------------------------------------------------------------------------------------------------------------

#debug: names(out1)
#debug: names(d.w)

out1.vars <- c("ID", "strain", "antibiotic", "run", "quality", "EUCAST", "MIC",
               "Estimate","Etest_predicted", "CI_low" ,"CI_up", "Hill", "upper", "lower")

out2.vars <- c("ID","preSIR.Etest_predicted", "susceptible", "resistant", "compare")

out <- merge(out1[, out1.vars],
             d.w[, out2.vars],
             by=c("ID"),
             all.x = T)
# save data
write.csv(out, "output/tables/categories_predicted.csv", row.names = F )

