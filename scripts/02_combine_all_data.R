
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file validation, reference, training data is merged

# set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

# Load data
## ----------------------------------------------------------------------------------------------------------------------------

#672
training <- read.csv("output/tables/training_data_normalized_parametertable.csv",stringsAsFactors = F)
#320
validation <- read.csv("output/tables/validation_data_normalized_parametertable.csv",stringsAsFactors = F)
#192
reference <- read.csv("output/tables/reference_data_normalized_parametertable.csv",stringsAsFactors = F)
#all data: 1184-1036=148
all <- rbind(training, validation, reference)
core <- rbind(training, validation)
# save data
## ----------------------------------------------------------------------------------------------------------------------------
all$quality[is.na(all$EstimateStderror) & all$quality=="quality ok"] <- "Error calculation fails"
write.csv(all,file="output/tables/alldata.csv", row.names = F)
write.csv(core,file="output/tables/training+validationdata.csv", row.names = F)