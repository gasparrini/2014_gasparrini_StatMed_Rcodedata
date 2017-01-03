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

####################################################################
# PREPARE THE DATA
####################################################################

# LOAD THE PACKAGES
library(survival) ; library(dlnm) ; library(splines)
library(foreign) ; library(xtable)

# CHECK VERSION OF THE PACKAGE
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")

# LOAD THE FUNCTIONS
source("functions.R")

# LOAD THE DATA
uminers <- read.csv("uminers.csv",row.names=1)
nrow(uminers)

# ORDER THE DATASET BY RECORD ID
uminers <- uminers[order(uminers$record),]
rownames(uminers) <- uminers$record

# CREATE MAIN DATASET
main <- uminers[,1:9]

################################################################################
# SPLITTING THE DATASET

# DETERMINE THE FAILURE TIMES
ftime <- sort(unique(main$ageexit[main$ind==1]))
# SPLIT THE DATASET
mainspl <- survSplit(Surv(agest,ageexit,ind)~., main, cut=ftime, start="agest")
nrow(mainspl)
# ORDER THE DATASET BY RECORD ID AND EXIT TIME
mainspl <- mainspl[order(mainspl$record,mainspl$ageexit),]

# CREATE THE RISK SETS
mainspl$riskset <- as.numeric(factor(mainspl$ageexit,levels=ftime))
# ELIMINATE THE CENSORED AT NOT-FAILURE TIMES
mainspl <- mainspl[mainspl$ageexit%in%ftime,]
nrow(mainspl)

################################################################################
# REDUCE THE DATA, SAMPLING FROM EACH RISK SET

# SAMPLE THE CASE + (MAX) 100 CONTROLS FROM EACH RISK SET
set.seed(13041975)
subset <- unlist(tapply(seq(nrow(mainspl)),mainspl$riskset, function(x) {
  case <- x[mainspl$ind[x]==1]
  controls <- sample(x[mainspl$ind[x]==0],min(length(x)-1,100))
  return(c(case,controls))
}))
mainsub <- mainspl[subset,]
# ORDER THE DATASET BY RECORD ID AND EXIT TIME
mainsub <- mainsub[order(mainsub$record,mainsub$ageexit),]
nrow(mainsub)

################################################################################
# CREATE THE FULL EXPOSURE HISTORIES FOR RADON AND SMOKING FROM 5-YEARS PERIODS

rexp <- t(apply(uminers,1,function(x) fexpfull(x[11],x[12],x[17:34])))
sexp <- t(apply(uminers,1,function(x) fexpfull(x[15],x[16],x[35:52])))

################################################################################
# COMPUTE THE MATRICES OF EXPOSURE HISTORIES FOR RADON AND SMOKING - LAG 2-40
# (NB: THE DATA ARE EXPECTED TO BE ALREADY ORDERED BY RECORD ID)

Qx <- do.call(rbind, lapply(seq(nrow(mainsub)),
  function(i) exphist(rexp[mainsub$record[i],],mainsub$ageexit[i],c(2,40))))

Qz <- do.call(rbind, lapply(seq(nrow(mainsub)),
  function(i) exphist(sexp[mainsub$record[i],],mainsub$ageexit[i],c(2,40))))

################################################################################
# GENERATE COVARIATES

# GENERATE DATE OF BIRTH (ASSUMED 15TH OF THE BIRTH MONTH)
bdate <- as.Date(paste(15,mainsub$bmon,paste(1,mainsub$byr,sep=""),sep="/"),
  format="%d/%m/%Y")
# GENERATE CALENDAR TIME (CENTERED ON 1970)
caltime <- as.numeric(bdate)/365.25+mainsub$ageexit
# GENERATE CROSS-BASIS FOR SMOKING: NS WITH 1 KNOT AT THE MEDIAN
# AND 2 STRATA FOR THE LAG (< and >= OF 20)
cbz <- crossbasis(Qz,lag=c(2,40),argvar=list(fun="ns",knots=2.5),
  arglag=list(fun="strata",breaks=20))

################################################################################
# CROSS-BASIS FOR RADON: DEFINE BASIS FOR LAG AND VAR FOR THE 14 MODELS

# DEFINE LIST OF ARGUMENTS FOR w(lag) 
laglist <- list(list(fun="strata",df=1,int=TRUE),
  list(fun="strata",breaks=1:3*10,int=TRUE),
  list(fun="bs",degree=2,knots=13.3,int=TRUE),
  list(fun="bs",degree=2,knots=13.3,int=FALSE))
laglist <- rep(laglist,2)

# DEFINE LIST OF ARGUMENTS FOR f(x)
varlist <- list(list(fun="lin"),list(fun="strata",breaks=c(26.7,60.2,122.2)),
  list(fun="bs",degree=2,knots=60.2))
varlist <- varlist[c(rep(1,4),3,2,3,3)]

#
