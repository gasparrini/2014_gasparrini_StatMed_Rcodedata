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


# FUNCTIONS TO COMPUTE THE LOGIT AND INVERSE LOGIT TRANSFORMATIONS
logit <- function(prob) log(prob/(1-prob))
invlogit <- function(linpred) exp(linpred)/(1+exp(linpred))

# METHOD FUNCTIONS AIC AND BIC FOR coxph OBJECTS
AIC.coxph <- function(object) -2*object$loglik[2] + 2*length(coef(object))
BIC.coxph <- function(object) -2*object$loglik[2] +
  length(coef(object))*log(object$nevent)

# FUNCTION TO PERFORM THE LIKELIHOOD RATIO TEST BETWEEN TWO COX MODELS
#   RETURNS A VECTOR WITH STATISTIC, DF AND P-VALUE
#   NO CHECK IS PERFORMED ABOUT THE NESTING
flrtest <- function(m1,m2) {
  ll <- abs(m1$loglik-m2$loglik)[2]
  df <- abs(length(coef(m1))-length(coef(m2)))
  return(round(c(stat=2*ll,df=df,pvalue=1-pchisq(2*ll,df)),3))
}

################################################################################
# FUNCTION TO RE-CREATE THE FULL EXPOSURE EXPERIENCE FROM 5-YEARS PERIODS, GIVEN:
#   - AGE AT START OF EXPOSURE (>5)
#   - AGE AT END OF EXPOSURE (<85)
#   - EXPOSURE LEVELS AT AGES 0-4, 5-9, ..., 85-89 (18 VALUES)
# STEPS:
#   - CREATE THE CUT-OFF POINTS
#   - CONVERT AGE AT START AND END IN YEARS WITH APPROPRIATE ROUNDING
#   - SUBSTITUTE THE APPROPRIATE CUT-OFF POINTS WITH AGE AT START AND END
#   - COMPUTE THE AVERAGE EXPOSURE FOR EACH YEAR (ADDING 10 MORE YEARS TO 99)
# NB: THE FUNCTION IS NOT VECTORIZED: ACCEPTS TWO SCALARS AND A VECTOR AS ARGS

fexpfull <- function(start,end,exp) {
  age <- 0:18*5
  start <- max(round(start-0.5),6)
  end <- min(round(end+0.5),84)
  age[(age-start)<0&(start-age)<5] <- start
  age[(age-end)>0&(age-end)<5] <- end
  expfull <- c(rep(as.numeric(exp)/diff(age),diff(age)),rep(0,10))
  names(expfull) <- 0:99+0.5
  return(expfull)
}

################################################################################

# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT OF EXPOSURE GIVEN:
#   - THE FULL EXPOSURE EXPERIENCE exp FROM TIME 1
#   - THE lag CONSIDERED
#   - THE TWO FUNCTIONS fvar AND flag DEFINING THE BI-DIMENSIONAL EFFECT
fcumeff <- function(exp,lag,fvar,flag) {
  exphistlist <- lapply(seq(length(exp)),fexphist,exp,lag)
  cumeff <- sapply(exphistlist, function(exphist){
    sum(do.call(fvar,list(exphist)) * do.call(flag,list(seq(lag[1],lag[2]))))
  })
  return(cumeff)
}

# FUNCTION TO SAMPLE THE EXIT TIME AND EVENT TYPE FOR A SUBJECT GIVEN:
#   - THE CUMULATIVE EFFECT AT EACH TIME
#   - BASELINE AND CENSORING RISK PROBABILITIES base AND censor (CONSTANT OR NOT)
# NB: THE FUNCTION FIRST COMPUTES THE CONDITIONAL PROBABILITIES AT EACH TIME,
#   THEN SAMPLES THE FIRST CENSORING OR EVENT TIME
# NB: THIS FUNCTION HAS BEEN SUPERSEDED BY permalgorithm IN PACKAGE PermAlgo
frsurv <- function(cumeff,pbase,pcensor) {
  pevent <- exp(log(pbase)+cumeff)
  if(any(pevent>1)) stop("Simulated effects generate event probability >1")
  tevent <- which.max(rbinom(101,1,c(pevent,1)))
  tcensor <- min(which.max(rbinom(101,1,c(pcensor,1))),100)
  return(if(tevent<tcensor) c(event=1,exit=tevent) else c(event=0,exit=tcensor))
}

################################################################################

# FUNCTIONS TO SIMULATE THE EXPOSURE-RESPONSE: LINEAR, PLATEAU, EXPONENTIAL
flin <- function(x) x/18
fplat <- function(x) (1-(1+x/1.5)/((1+x/1.5)^2))/1.9
fexp <- function(x) (exp(x/3.5)-1)/exp(10/3.5)/1.1

# FUNCTIONS TO SIMULATE THE LAG STRUCTURE: CONSTANT, DECAY, WITH PEAK
fconst <- function(lag) lag-lag+0.20
fdecay <- function(lag) exp(-lag/9)*0.95
fpeak <- function(lag) 8.5*dnorm(lag,10,7)

################################################################################

# ARGUMENTS FOR 3D PLOTS
arg3D <- list(x=seq(0,10,0.25),y=0:20*2,ticktype="detailed",theta=230,
  ltheta=200,phi=30,lphi=30,xlab="Exposure",ylab="Lag",zlab="HR",
  shade = 0.75,r=sqrt(3),d=5,cex.axis=0.7,cex.lab=0.8,border=grey(0.3),
  col=grey(0.99))

#
