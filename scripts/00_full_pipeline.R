##############################################################################################################

#-------------------------------------------------------------------------------------------------------------
#-----------------------------set working directory-----------------------------------------------------------

# setwd("~/your_folder")

# ############################################################################################################

#-------------------------------------------------------------------------------------------------------------
source('scripts/01_fit_dose-response_model.R')
# In this file we fit the dose-response model. 
# Validation, training and reference data are stored in different folders, and processed susequently in the loop. 
# input files are in folders:
  # data/reference_data
  # data/training_data
  # data/validation_data
# output files are in
  # output/tables/  and end with _normalized_parametertable.csv and _normalized_model_list.rds

#-------------------------------------------------------------------------------------------------------------
source('scripts/02_combine_all_data.R')
# this scripts just combines the data from the previous script
# output file are
 # output/tables/alldata.csv"
 # output/tables/training+validationdata.csv"

#-------------------------------------------------------------------------------------------------------------
source('scripts/03_linear_regression_training_data.R')
# linear regression using output/tables/alldata.csv (using training data only)
# output files
 # output/tables/lm_parameters_variance_covariance_matrix_list.rds
 # output/tables/pearson.rds
 # output/figures/regression_plot.RDS

#-------------------------------------------------------------------------------------------------------------
source('scripts/04_predict_MICs_all_data.R')
# predict MIC for all data and compute bootstrap confidence intervals
# input files
 # output/tables/alldata.csv
 # output/tables/lm_parameters_variance_covariance_matrix_list.rds
# output
 # output/tables/Predicted_Etest_all_data.csv

#-------------------------------------------------------------------------------------------------------------
source('scripts/05_predicted_categories_graph.R')
# graph with the categories predicted by the model: Figure 3 (graph 4 single anti are saved as well)
# also saves the file
# output/tables/categories_predicted.csv

#-------------------------------------------------------------------------------------------------------------
source('scripts/06_deviations_boxplot_and_density.R')
# this file produces a couple of graphs:
# - boxplot with deviations from Etest
# - density of deviations from Etest

#-------------------------------------------------------------------------------------------------------------
source('scripts/07_all_curves_plotted_graph.R')
# In this file we want plot all fitted curves, in antibiotic specific panels
# input
 # output/tables/categories_predicted_training.csv
 # output/tables/training_data_normalized_model_list.rds
 # output/tables/validation_data_normalized_model_list.rds
# output
 # output/figures/Figure1.pdf

#-------------------------------------------------------------------------------------------------------------
source('scripts/08_supplement_timecourse.R')
# time course figure in the supplementary material:
# output/figures/FigureS1.pdf

#-------------------------------------------------------------------------------------------------------------
source('scripts/09_supplement_coefficient_of_variation.R')
# figure 2 in supplementary material
# output/figures/FigureS2.pdf

#-------------------------------------------------------------------------------------------------------------
source('scripts/10_supplement_hill_coefficient.R')
# figure 3 in supplementary material 
# output/figures/FigureS3AB.pdf

#-------------------------------------------------------------------------------------------------------------
source('scripts/11_figures_panel.R')
# produces figure 2 in the paper
# output/figures/Figure2.pdf

#-------------------------------------------------------------------------------------------------------------
source('scripts/99_session_info.R')
# get session info for html markdown file
