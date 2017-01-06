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
# SIMPLE EXAMPLE
####################################################################

# PREPARE THE DATA
source("01.prep.R")

# THE MAIN OBJECTS CREATED FOR THE ANALYSIS ARE:
# - THE DATASET: main 
# - A MATRIX OF EXPOSURE HISTORIES FOR RADON: Qx
# - A MATRIX OF EXPOSURE HISTORIES FOR SMOKING: Qz
# - THE CROSS-BASIS FOR SMOKING (IDENTICAL TO CUMULATIVE EXPOSURE): cbz
# - THE VARIABLE CALENDAR TIME: caltime
# - METHOD FUNCTIONS TO COMPUTE AIC AND BIC: AIC.coxph, BIC.coxph

####################################################################
# SEE:

help(dlnm)
help(crossbasis)
help(crosspred)
help (plot.crosspred)

####################################################################
# DLMs (WITH A LINEAR FUNCTION FOR f(x))

# MODEL FOR UN-WEIGHTED CUMULATIVE EXPOSURE (MODEL1 IN THE MANUSCRIPT)
# DEFINE THE CROSS-BASIS
cbx1 <- crossbasis(Qx,lag=c(2,40),argvar=list(fun="lin"),
  arglag=list(fun="strata",df=1))
summary(cbx1)

# RUN THE MODEL INCLUDING ALSO THE CROSS-BASIS FOR SMOKING AND CALENDAR YEAR
model1 <- coxph(Surv(agest,ageexit,ind)~cbx1+cbz+caltime,mainsub)
AIC(model1)

# OBTAIN THE PREDICTED RISK FOR A SEQUENCE OF RADON LEVELS
pred1 <- crosspred(cbx1,model1,at=0:25*10,cen=0)

# 3D GRAPH FOR MODEL1
plot(pred1,"3d",xlab="WLM/year",ylab="Lag (years)",zlab="HR")

# NOW A MODEL WITH A QUADRATIC B-SPLINE FOR w(lag), NO INTERCEPT (MODEL4)
cbx4 <- crossbasis(Qx,lag=c(2,40),argvar=list(fun="lin"),
  arglag=list(fun="bs",degree=2,knots=13.3,int=FALSE))
summary(cbx4)
model4 <- coxph(Surv(agest,ageexit,ind)~cbx4+cbz+caltime,mainsub)
AIC(model4)
pred4 <- crosspred(cbx4,model4,at=0:25*10,cen=0)

# 3D GRAPH FOR MODEL4, WITH NON-DEFAULT PERSPECTIVE
plot(pred4,"3d",xlab="WLM/year",ylab="Lag (years)",zlab="HR",theta=240)

# COMPARISON OF LAG CURVES FOR 100 WLM/year (WITH NON-DEFAULT LINES CI)
#   (SIMILAR TO FIGURE 2 IN THE MANUSCRIPT)
plot(pred4,var=100,xlab="Lag (years)",ylab="HR for 100 WLM/year",ci="lines")
lines(pred1,var=100,col=4,ci="lines")
legend("top",c("Model 4","Model 1"),lty=1,col=c(2,4),inset=0.1)

# COMPARISON OF EXPOSURE RESPONSE AT LAG 15 (WITH NON-DEFAULT LINES CI)
plot(pred4,lag=15,xlab="WLM/years",ylab="HR at lag 15",ci="lines",ylim=c(.95,1.15))
lines(pred1,lag=15,col=4,ci="lines")
legend("topleft",c("Model 4","Model 1"),lty=1,col=c(2,4),inset=0.1)

# HR FOR 100 WLM/year AT LAG 11 (WITH CI) FROM MODEL4
pred4$matRRfit["100","lag11"]
pred4$matRRlow["100","lag11"]
pred4$matRRhigh["100","lag11"]

####################################################################
# DLNMs (WITH A QUADRATIC B-SPLINE FUNCTION FOR f(x))

# MODEL8 IN THE MANUSCRIPT
cbx8 <- crossbasis(Qx,lag=c(2,40),argvar=list(fun="bs",degree=2,knots=59.4),
  arglag=list(fun="bs",degree=2,knots=13.3,int=F))
summary(cbx8)
model8 <- coxph(Surv(agest,ageexit,ind)~cbx8+cbz+caltime,mainsub)
AIC(model8)
pred8 <- crosspred(cbx8,model8,at=0:25*10,cen=0)

# PLOTS AS ABOVE (SEE FIGURE 3 IN THE MANUSCRIPT)
plot(pred8,"3d",xlab="WLM/year",ylab="Lag (years)",zlab="HR")
plot(pred8,var=100,xlab="Lag (years)",ylab="HR for 100 WLM/year")
plot(pred8,lag=15,xlab="WLM/years",ylab="HR at lag 15",ylim=c(.9,1.35))

# HR FOR 100 WLM/year AT LAG 11 (WITH CI) FROM MODEL8
pred8$matRRfit["100","lag11"]

# LR TEST FOR MODEL4 VS. MODEL8
flrtest(model4,model8)

####################################################################
# PREDICTIONS

# CREATE AN EXPOSURE PROFILE: 20 WLM/year FOR 10 YEARS, UNEXPOSED FOR 10 MORE
exp <- rep(c(20,0),c(10,10))

# CREATE THE EXPOSURE HISTORY USED FOR PREDICTION
hist <- exphist(exp,time=20,lag=c(2,40))

# CREATE PREDICTION
predhist8 <- crosspred(cbx8,model8,at=hist,cen=0)

# OVERALL CUMULATIVE EFFECT (SEE TABLE 3 IN THE MANUSCRIPT)
with(predhist8,cbind(allRRfit,allRRlow,allRRhigh))

# CREATE A MATRIX OF EXPOSURE HISTORIES FOR A SEQUENCE OF YEARS
# EXPOSURE PROFILE: 20 WLM/year FOR THE FIRST 15 YEARS
# EXPOSURE HISTORIES COMPUTED FOR YEAR 1 TO 60
profile <- rep(c(20,0),c(15,45))
dynhist <- exphist(profile,lag=c(2,40))

# CREATE PREDICTION
predhist8 <- crosspred(cbx8,model8,at=dynhist,cen=0)

# PLOT FOR THE TREND IN RISK ALONG TIME (FIGURE 5 IN THE MANUSCRIPT)
plot(predhist8,"overall",xlab="Years",ylab="HR",ylim=c(0,4.5))
abline(v=15,lty=2)

#
