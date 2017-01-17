
# Description
## ----------------------------------------------------------------------------------------------------------------------------
# In this file we fit the dose-response model. 
# Validation, training and reference data are stored in different folders, and processed susequently in the loop. 

# Load dependencies and set options
## ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

packages <- c("drc")
to.install <- setdiff(packages, rownames(installed.packages()))
if (length(to.install) > 0) {
  install.packages(to.install)
}
lapply(packages, library, character.only = TRUE)


# outer loop: process different types of data
## ----------------------------------------------------------------------------------------------------------------------------

for(data_type in c("reference_data", "training_data", "validation_data")){
#for(data_type in c("validation_data")){
    
# debug:   data_type <- "reference_data"
# data_type <- "validation_data"
# data_type <- "validation_data"
# set data type variables  (mainly used to import and save data)
input_names_path <- paste0("data/", data_type )
nome <- gsub("_data$", "", data_type)
csv.out <- paste0("output/tables/", data_type, "_normalized_parametertable.csv") 
rds.out <- paste0("output/tables/", data_type, "_normalized_model_list.rds") 

# Load file names
## -----------------------------------------------------------------------------------------------------------
input=list.files(path=input_names_path, pattern=".txt")

# fit model 
## -----------------------------------------------------------------------------------------------------------
# The model is firstly fitted on non normalized data.
# Subsequently, the model parameters are used to normalized data.
# Finally, the EC50 is estimated by re-fitting the model on the normalized data.

# empty list where model fitter to normalized data will be saved and empty parmist
model_list_export <- list()
parmlist=NULL

for (i in unique(input)){
  # debug: i="WHO_F_Azithromycin_1.12.15.exp1.txt"
  # i="1_Ceftriaxone_40strains1.txt"
  # i="19_Ceftriaxone_40strains2.txt"
  print("-------------------------------------------------------------------------------------------------------------------")
  print(data_type)
  print(i)
  print("loop1")
  print("-------------------------------------------------------------------------------------------------------------------")

  newdf <- read.table(paste0(input_names_path, "/", i), header=T)

  nummer <- length(newdf$response)
  newdf$MIC <- as.character(newdf$MIC)
  newdf$concentration[1] <- newdf$concentration[2]*10
  newdf$concentration[nummer] <- newdf$concentration[nummer-1]/10
  #newdf$concentration[12] <- newdf$concentration[11]/10
  
  ab <- as.character(newdf$antibiotic[1])
  str <- as.character(newdf$strain[1])
  
  # subtract background
  posctr <- subset(newdf,newdf$concentration==max(newdf$concentration))
  newdf$response <- newdf$response-posctr$response
  
  # check limits of detection
  LLO_u <- newdf$response[nummer]*0.5 < newdf$response[2]
  LLO_l <- newdf$response[nummer]*0.5 > newdf$response[nummer-1]
  
  # fit model without normalization
  tryCatch({ 
    model=try(drm(formula = newdf$response ~ log(newdf$concentration), fct = L.4()))
    model_original=try(drm(formula = newdf$response ~ newdf$concentration, fct = LL.4()))
    if(isTRUE(class(model)=="try-error")){
      print("loop2")
      allconc = subset(newdf, newdf$concentration != Inf)
      Estimate1 = as.numeric(NA)
      quality="poor quality"
      if(LLO_u==TRUE){quality="above limit of detection"
                      Estimate1 = log(max(allconc$concentration))}
      if(LLO_l==TRUE){quality="poor quality"
                      Estimate1 = as.numeric(NA)}
      
      Liste1 = data.frame(ID=i,Estimate=Estimate1,
                          EstimateStderror=as.numeric(NA),
                          EC50=as.numeric(NA),
                          Hill=as.numeric(NA),
                          upper=as.numeric(NA),
                          lower=as.numeric(NA),
                          strain=str,
                          antibiotic=ab,
                          MIC=newdf$MIC[1],
                          EUCAST=newdf$eucast[1],
                          quality=quality,
                          stringsAsFactors = FALSE,
                          run=nome)
      
      parmlist = rbind.data.frame(Liste1,parmlist)
      Liste1 <- NULL
      next
    } else {
      print("loop3")
      testdf=data.frame(PR(model,newdf$concentration))
      names(testdf) <- "fitted"
      testdf$dose <- as.numeric(row.names(testdf))
      normalized=(testdf$fitted/coefficients(model)[3])*100
      normalized_points=(newdf$response/coefficients(model)[3])*100
      normalized_points_original=(newdf$response/coefficients(model_original)[3])*100
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # fit model on normalized data
  tryCatch({
    print("loop4")
    model_norm=try(drm(formula = normalized_points ~ log(newdf$concentration), fct = L.4()))
    model_norm_original=try(drm(formula = normalized_points_original ~ newdf$concentration, fct = LL.4()))
    summary(model_norm)
    if(isTRUE(class(model_norm)=="try-error")){
      print("loop5")
      Liste2 = data.frame(ID=i,
                          Estimate=as.numeric(NA),
                          EstimateStderror=as.numeric(NA),
                          Hill=as.numeric(NA),
                          upper=as.numeric(NA),
                          lower=as.numeric(NA),
                          strain=str,
                          antibiotic=ab,
                          MIC=newdf$MIC[1],
                          EUCAST=newdf$eucast[1],
                          quality="normalization failed",
                          run=nome)
      parmlist=rbind.data.frame(Liste2,parmlist)
      Liste2 <- NULL
      next
    } else {
      print("loop6")
      z=50
      ECX <- NULL
      original <- NULL
      allconc <- NULL
      allconc = subset(newdf, newdf$concentration != Inf)
      ECX=ED(model_norm,z,interval=c("delta"))
      original=ED(model_norm_original,z,interval=c("delta"))
      ECY=ED(model,z,interval=c("delta"),logBase = exp(1))
      Estimate3=ECX[1]
      originalEC=original[1]
      Estimate_raw=ECY[1]
      Error_raw=ECY[2]
      Error=ECX[2]
      originalECerror=original[2]
      CI = confint(model_norm, level = 0.95)
      CI_l=CI[4]
      CI_u=CI[8]
      Hill=coefficients(model_norm)[1]
      Hill_raw=coefficients(model)[1]
      upper=coefficients(model_norm)[3]
      lower=coefficients(model_norm)[2]
      EC50=coefficients(model_norm)[4]
      
      if(LLO_u==TRUE){quality="above limit of detection"
                      Estimate3 = log(max(allconc$concentration))}
      if(LLO_l==TRUE){quality="poor quality"
                      Estimate3 = as.numeric(NA)}
      if(LLO_l==FALSE & LLO_u==FALSE){quality="quality ok"}
      # if(upper<=80){quality="poor quality"}
      
      Liste3 = data.frame(ID=i,
                          Estimate=Estimate3,
                          EstimateStderror=Error,
                          Hill=Hill,
                          upper=upper,
                          lower=lower,
                          strain=str,
                          antibiotic=ab,
                          MIC=newdf$MIC[1],
                          EUCAST=newdf$eucast[1],
                          quality=quality,
                          run=nome)
      parmlist=rbind.data.frame(Liste3,parmlist)
      Liste3 <- NULL
      
      # save full model so that we can use it at a later stage
      model_list_export[i] <- list("model"=model_norm)
    }
    LLO_l <- NULL
    LLO_u <- NULL
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# save data
## -----------------------------------------------------------------------------------------------------------
write.csv(parmlist, file= csv.out,row.names = F)
saveRDS(model_list_export, rds.out)

}

# file end
## -----------------------------------------------------------------------------------------------------------



