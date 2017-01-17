

# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file we want to validate the method by predicting the E-Test starting from the MIC and its sd from the dose-response 
# model. We do this by bootstrap, because we need to take in account the sd of MIC from the dose-response.
# For each antibiotic we use the linear regression parameters variance~covariance matrix to sample new estimates of the lm model.
# Similarly, we sample new estimates for the MIC. We combine the two and get confidence intervals for the predicted E-Test
# Linear regression parameters are the same for each antibiotic.


# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("plyr", "mvtnorm")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)
options(scipen = 4)

# Load data
## ----------------------------------------------------------------------------------------------------------------------------
df <- read.csv("output/tables/alldata.csv",stringsAsFactors = F)
parameters <- readRDS("output/tables/lm_parameters_variance_covariance_matrix_list.rds")

# count data just to make sure everything is alright
su1 <- ddply(df, ~antibiotic + strain, summarise, count=length(antibiotic))
su2 <- ddply(su1, ~antibiotic , summarise, count=length(antibiotic))


# we make compute the CI by resampling both the
# estimates of the MIC and the estimates of the linear regression.
# Note that we have the same regression parameters for all observations.
## ----------------------------------------------------------------------------------------------------------------------------
size <- 1e5 # size of the random samples for bootstrapping

# First we get the antibiotic specific parameters (parameters are no antibiotic specific anymore, as there is no max correlation anymore)
coefficients <- coef(parameters[[1]]$Estimates)
intercept <- coefficients[[1]]
slope <- coefficients[[2]]
variance.covariance <- parameters[[1]]$Matrix[[1]]

# check if everything is ok: Estimate and its sd cannot be NA!!
if(any(is.na(df$Estimate[df$quality=="quality ok"]))==T | 
   any(is.na(df$EstimateStderror[df$quality=="quality ok"])==T)
   )stop("There are NAs in either the Estimate or sd. These obs cannot be included in the loop")
if(any(duplicated(df$ID))==T)stop("ID is not unique!!")

# Create the vars where we will save the confidence intervals
df$Etest_predicted_log <- NA
df$Etest_predicted <- NA
df$CI_low <- NA
df$CI_up <- NA


for(ID in unique(df$ID[df$quality=="quality ok"])){ ## we only do the boostrap for the quality ok
  ## debug:
  # ID <- "27_Azithromycin_80strains12.txt"
  # ID <- "65_Penicillin_80strains8.txt"
  # ID <- "22_Gentamicin_80strains3.txt"  sd=NA
  # ID <- "WHO_L_Gentamicin_4.12.15.exp3.txt"  sd=NA
  # ID <- "60_Gentamicin_80strains8.txt"  sd=NA

  MIC.estimated <- df$Estimate[df$ID==ID]
  MIC.estimated.sd <- df$EstimateStderror[df$ID==ID] 
  # get predicted value

  ETest_predicted_log <- intercept + slope * (MIC.estimated)  
  ETest_predicted <- exp(ETest_predicted_log)
  
  # now we want the confidence interval:
  
  # sample the EC50
  MIC.sample_log <- rnorm(size, MIC.estimated, MIC.estimated.sd)

  # sample the coefficients
  linear_coeff.sample <-  rmvnorm(size, mean = coefficients, sigma = variance.covariance )

  ETest_sample_log <- linear_coeff.sample[,1] + linear_coeff.sample[,2] * MIC.sample_log
  ETest_sample <- exp(ETest_sample_log)
  #hist(ETest_sample)
  
  # now we can compute the confidence intervals via cut-off
  ETest_sample_CI <- quantile(ETest_sample, probs = c(0.025, 0.975),na.rm=T)
  
  # save the results
  df$Etest_predicted_log[df$ID==ID] <- ETest_predicted_log  
  df$Etest_predicted[df$ID==ID] <- ETest_predicted
  df$CI_low[df$ID==ID] <- ETest_sample_CI[1]
  df$CI_up[df$ID==ID] <-  ETest_sample_CI[2]
}
  

# fix double observations problem: should do this before putting data in the regression
## ----------------------------------------------------------------------------------------------------------------------------
df$quality[grepl("<", df$MIC)] <- "Etest below limit of detection"
df$quality[grepl(">", df$MIC)] <- "Etest above limit of detection"


# fill in values for above limit of detection
## ----------------------------------------------------------------------------------------------------------------------------
df$Etest_predicted[grepl("limit of detection", 
                         df$quality)] <- exp(df$Estimate[grepl("limit of detection", df$quality)])
df$Etest_predicted_log[grepl("limit of detection",
                             df$quality)] <- log(df$Etest_predicted[grepl("limit of detection", df$quality)])

# save
## ----------------------------------------------------------------------------------------------------------------------------
write.csv(df,file="output/tables/Predicted_Etest_all_data.csv", row.names = F)


# file end
## -----------------------------------------------------------------------------------------------------------














