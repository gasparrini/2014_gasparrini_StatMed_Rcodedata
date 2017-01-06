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

################################################################################
# SIMULATE THE EXPOSURE

# SET THE SEED
set.seed(13041975)

# SIMULATE THE EXPOSURE PROFILES FROM 1 TO 100
# A SEQUENCE OF 5-14 EXPOSURE EVENTS WITH LENGTH 1-10 AT INTENSITY 0-10
# THE VALUES OF EXPOSURE OCCURRENCES ARE RESAMPLED IF >10
expsim <- t(sapply(seq(nsub), function(nsub) {
  hist <- 0
  for(i in seq(15)) {
    hist0 <- rep(0,100)
    start <- sample(0:95,1)
    length <- sample(1:5,1)
    val <- runif(1,0,10)
    hist0[start:(start+length)] <- val
    hist <- hist+hist0
  }
  return(hist)
}))

rownames(expsim) <- seq(nsub)
expsim[expsim>10] <- runif(sum(expsim>10),0,10)

################################################################################
# DEFINE THE BI-DIMENSIONAL ASSOCIATIONS USED FOR SIMULATING DATA

# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- cbind(x=rep(c("flin","fplat","fexp"),each=3),
  lag=rep(c("fconst","fdecay","fpeak"),3))
rownames(combsim) <- paste(rep(c("Linear","Plateau","Exponential"),each=3),
  rep(c("Constant","Decay","Peak"),3),sep="-")


# LIST WITH TRUE EFFECT SURFACES OF FOR EACH COMBINATIONS
effsimlist <- lapply(seq(nrow(combsim)), function(j) {
  effsim <- t(sapply(seq(0,10,0.25), function(x) {
    do.call(combsim[j,1],list(x)) * do.call(combsim[j,2],list(0:20*2))  
  }))
  dimnames(effsim) <- list(seq(0,10,0.25),paste("lag",0:20*2,sep=""))
  return(effsim)
})
names(effsimlist) <- rownames(combsim)

################################################################################
# DEFINTIONS OF ARGUMENTS FOR CROSS-BASIS FUNCTIONS USED FOR ESTIMATION

# DEFINE A LIST OF LISTS OF ARGUMENTS FOR f(x)
svarlist <- list(list(fun="lin"),
  list(fun="bs",degree=2,Bound=c(0,10)),
  list(fun="bs",degree=2,knots=3.3,Bound=c(0,10)),
  list(fun="bs",degree=2,knots=5,Bound=c(0,10)),
  list(fun="bs",degree=2,knots=6.7,Bound=c(0,10)),
  list(fun="bs",degree=2,knots=c(3.3,6.7),Bound=c(0,10)))
svarlist <- rep(svarlist,length(svarlist))

# DEFINE A LIST OF LISTS OF ARGUMENTS FOR w(lag)
slaglist <- list(list(fun="strata",df=1,int=T),
  list(fun="bs",degree=2,int=T),
  list(fun="bs",degree=2,knots=13.3,int=T),
  list(fun="bs",degree=2,knots=20,int=T),
  list(fun="bs",degree=2,knots=26.7,int=T),
  list(fun="bs",degree=2,knots=c(13.3,26.7),int=T))
slaglist <- rep(slaglist,each=length(slaglist))

################################################################################
# CREATE THE OBJECT TO STORE RESULTS

# LIST OF MODELS
modellist <- vector("list",length(svarlist))

# LISTS FOR STORING THE AVERAGE PREDICTION, COVERAGE AND SAMPLES
# FOR THE WHOLE SURFACE
aictotpred <- vector("list",nrow(combsim))
names(aictotpred) <- rownames(combsim)
bictotpred <- aictotcov <- bictotcov <- 
  aictotsample <- bictotsample <- aictotpred

# MATRICES FOR STORING PREDICTED VALUE, COVERAGE AND TRUE VALUE
aicpred <- bicpred <- aiccov <- biccov <- matrix(NA,nsim,nrow(combsim))
trueeff <- matrix(NA,nsim,nrow(combsim))

# MATRICES FOR STORING THE DF USED IN EACH DIMENSION FOR SELECTED MODELS
aicdfvar <- aicdflag <- bicdfvar <- bicdflag <- matrix(NA,nsim,nrow(combsim))

# MATRIX FOR STORING THE NUMBER OF EVENTS
navg <- matrix(NA,nsim,nrow(combsim))

#
