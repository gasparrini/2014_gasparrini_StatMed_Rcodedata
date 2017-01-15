################################################################################
# Updated version of the R code for the analysis in:
#
#   "Modeling exposure-lag-response associations with distributed lag 
#     non-linear models"
#   Antonio Gasparrini
#   Statistics in Medicine - 2014
#   http://www.ag-myresearch.com/2014_gasparrini_statmed.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2014_gasparrini_StatMed_Rcodedata
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
