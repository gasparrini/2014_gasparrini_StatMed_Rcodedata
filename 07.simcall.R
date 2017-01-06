################################################################################
# Updated version of the R code for the analysis in:
#
#   "Modeling exposure-lag-response associations with distributed lag 
#     non-linear models"
#   Antonio Gasparrini
#   Statistics in Medicine - 2014
#   http://www.ag-myresearch.com/statmed2014.html
#
# Update: 14 December 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

# LOAD PACKAGES
#library(dlnm,lib.loc='/users/emsuagas/R/library')
#library(PermAlgo,lib.loc='/users/emsuagas/R/library')
library(dlnm) ; library(PermAlgo)
library(survival) ; library(foreign)

################################################################################
# DEFINE THE SETTING FOR THE SIMULATION

# NUMBER OF SUBJECTS, NUMBER OF ITERATIONS
nsub <- 400
nsim <- 3

# NUMBER OF SAMPLES USED AS EXAMPLES OF INDIVIDUAL ESTIMATES
nsample <- min(nsim,25)

# NOMINAL VALUE
qn <- qnorm(0.975)

################################################################################
# RUN THE SIMULATION

source("functions.R")
source("simprep.R")
source("simrun.R")

################################################################################
# SAVE

save.image("simul_400.RData")

#
