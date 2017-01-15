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

####################################################################
# CONSTRAINTS ON THE LAG CURVE
####################################################################

# DEFINE A NEW FUNCTION FOR RIGHT-CONSTRAINED SPLINES
bsrc <- function(x,red=2,df=NULL,knots=NULL,degree=3,intercept=FALSE,
   Boundary.knots=range(x)) {
  # CREATE THE B-SPLINE BASIS
  basis <- bs(x,df=df,knots=knots,degree=degree,intercept=intercept,
   Boundary.knots=Boundary.knots)
  # SELECT THE ARGUMENTS DEFINING bs
  attr <- attributes(basis)[match(names(formals(bs)),names(attributes(basis)),
    nomatch=0)]
  if(ncol(basis)<=red) stop("dimension of the basis not compatible with 'red'")
  # REMOVE THE SELECTED COLUMNS
  basis <- basis[,-seq(ncol(basis)-red+1,ncol(basis)),drop=FALSE]
  # RE-GENERATE ATTRIBUTES
  attributes(basis) <- c(attributes(basis),list(red=red),attr)
  return(basis)
}

# GENERATES THE CROSS-BASIS USING THE USER-DEFINED FUNCTION
cbxrc <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[8]],
  arglag=list(fun="bsrc",degree=2,knots=c(13.3,26.7),int=F,red=2))
summary(cbxrc)

# VISUALIZE THE B-SPLINE
pdf("bspline.pdf",height=5,width=8)
matplot(0:40,bs(0:40,degree=2,knots=c(13.3,26.6)),type="l",lty=c(1,2,4,5),lwd=1.5,
  xlab="Lag",ylab="Value")
dev.off()
# THE LAST TWO VARIABLES NEED TO BE EXCLUDED

# RUN THE MODEL EXCLUDING THE LAST 6 CROSS-BASIS VARIABLES
modelrc <- coxph(Surv(agest,ageexit,ind)~cbxrc+cbz+caltime,mainsub)
predrc <- crosspred(cbxrc,modelrc,at=0:25*10,cen=0)

AIC(modelrc) ; AIC(model8)

################################################################################
# PLOT

pdf("constr.pdf",height=4.1,width=12)

layout(matrix(1:2,ncol=2,byrow=TRUE))
par(mar=c(2,2,1.5,0.5))

plot(predrc,theta=230,cex.axis=0.7,cex.lab=0.8,nticks=9,ltheta=230,lphi=45,
  border=grey(0.3),xlab="WLM/year",ylab="Lag (years)",zlab="HR",col=grey(0.99),
  zlim=c(0.95,1.3))
mtext("Right and left-constrained",cex=0.8)

par(mar=c(4,4,0.5,0.5))

plot(predrc,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5)
lines(pred8,var=100,lwd=1.5,lty=2,col=1)
legend("topright",c("Right and left-constrained","Left-constrained (Model 8)"),
  lty=1:2,lwd=1.5,col=1,cex=0.9,bty="n",inset=0.1)

dev.off()

################################################################################
# SENSITIVITY ANALYSIS: FITTING MODELS ON THE SUB-GROUP WITH MAX WLM/YEAR < 300
################################################################################

# DEFINE THE SUB-GROUP
sub <- apply(Qx,1,max)<300
# SUBJECT IN THE SUB-GROUP
length(unique(mainsub[sub,]$record))
length(unique(mainsub[sub,]$record))/nrow(main)

# RUN THE MODELS WITH THE SUBSET OF DATA
cbx4sub <- crossbasis(Qx[sub,],lag=c(2,40),argvar=varlist[[4]],
  arglag=laglist[[4]])
model4sub <- coxph(Surv(agest,ageexit,ind)~cbx4sub+cbz[sub,]+caltime[sub],
  mainsub[sub,])
pred4sub <- crosspred(cbx4sub,model4sub,at=0:25*10,cen=0)

cbx8sub <- crossbasis(Qx[sub,],lag=c(2,40),argvar=varlist[[8]],
  arglag=laglist[[8]])
model8sub <- coxph(Surv(agest,ageexit,ind)~cbx8sub+cbz[sub,]+caltime[sub],
  mainsub[sub,])
pred8sub <- crosspred(cbx8sub,model8sub,at=0:25*10,cen=0)

AIC(model4sub) ; AIC(model8sub)

################################################################################
# PLOT

pdf("subset.pdf",height=4.1,width=12)

layout(matrix(1:2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))

plot(pred8,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5)
lines(pred8sub,var=100,lwd=1.5,lty=2,col=1)
lines(pred4sub,var=100,lwd=1.5,lty=4,col=1)
legend("topright",c("Model 8","Model 8 (subset)","Model 4 (subset)"),
  lty=c(1:2,4),lwd=1.5,col=1,cex=0.9,bty="n",inset=0.02)

plot(pred8,lag=15,xlab="WLM/year",ylab="HR at lag 15",ylim=c(0.9,1.3),
  lab=c(8,5,5),col=1,lwd=1.5)
lines(pred8sub,lag=15,lwd=1.5,lty=2,col=1)
lines(pred4sub,lag=15,lwd=1.5,lty=4,col=1)
legend("topleft",c("Model 8","Model 8 (subset)","Model 4 (subset)"),
  lty=c(1:2,4),lwd=1.5,col=1,cex=0.9,bty="n",inset=0.02)

dev.off()

################################################################################
# SENSITIVITY ANALYSIS: TAKING A LOG TRANSFORMATION OF WLM
################################################################################

# DEFINE LOG FUNCTION
mylog <- function(x) log(x+1)

# GENERATES THE CROSS-BASIS USING THE USER-DEFINED FUNCTION
cbxlog <- crossbasis(Qx,lag=c(2,40),argvar=list(fun="mylog"),
  arglag=laglist[[4]])
summary(cbxlog)

# RUN THE MODEL AND OBTAIN PREDICTIONS
modellog <- coxph(Surv(agest,ageexit,ind)~cbxlog+cbz+caltime,mainsub)
predlog <- crosspred(cbxlog,modellog,at=0:25*10,cen=0)

AIC(model8) ; AIC(modellog)

################################################################################
# PLOT

pdf("logscale.pdf",height=4.1,width=12)

layout(matrix(1:2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))

plot(pred8,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5)
lines(predlog,var=100,lwd=1.5,lty=2,col=1)
legend("topright",c("Model 8","Log function"),lty=c(1:2),lwd=1.5,col=1,
  cex=0.9,bty="n",inset=0.02)

plot(pred8,lag=15,xlab="WLM/year",ylab="HR at lag 15",ylim=c(0.9,1.3),
  lab=c(8,5,5),col=1,lwd=1.5)
lines(predlog,lag=15,lwd=1.5,lty=2,col=1)
legend("topleft",c("Model 8","Log function"),lty=c(1:2,4),lwd=1.5,col=1,
  cex=0.9,bty="n",inset=0.02)

dev.off()

################################################################################
# SENSITIVITY ANALYSIS ON KNOT LOCATION
################################################################################

# RUN VARIATIONS OF MODEL 8 WITH KNOTS AT DIFFERENT POSITIONS
laglist[[8]]
cbx8b <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[8]],
  arglag=list(fun="bs",df=3,degree=2,knots=20,int=FALSE))
summary(cbx8b)
model8b <- coxph(Surv(agest,ageexit,ind)~cbx8b+cbz+caltime,mainsub)
pred8b <- crosspred(cbx8b,model8b,at=0:25*10,cen=0)

cbx8c <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[8]],
  arglag=list(fun="bs",df=3,degree=2,knots=26.6,int=FALSE))
summary(cbx8c)
model8c <- coxph(Surv(agest,ageexit,ind)~cbx8c+cbz+caltime,mainsub)
pred8c <- crosspred(cbx8c,model8c,at=0:25*10,cen=0)

cbx8d <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[8]],
  arglag=list(fun="bs",df=3,degree=2,knots=c(13.3,26.6),int=FALSE))
summary(cbx8d)
model8d <- coxph(Surv(agest,ageexit,ind)~cbx8d+cbz+caltime,mainsub)
pred8d <- crosspred(cbx8d,model8d,at=0:25*10,cen=0)

AIC(model8) ; AIC(model8b) ; AIC(model8c) ; AIC(model8d)

################################################################################
# PLOT

pdf("knotloc.pdf",height=4.1,width=12)

layout(matrix(1:2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))

plot(pred8,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5)
lines(pred8b,var=100,lwd=1.5,lty=2,col=1)
lines(pred8c,var=100,lwd=1.5,lty=4,col=1)
lines(pred8d,var=100,lwd=1.5,lty=5,col=1)
legend("topright",c("Knots: 13.3","Knots: 20","Knots: 26.6","Knots: 13.3,26.6"),
  lty=c(1:2,4:5),lwd=1.5,col=1,cex=0.9,bty="n",inset=0.02)

plot(pred8,lag=15,xlab="WLM/year",ylab="HR at lag 15",ylim=c(0.9,1.3),
  lab=c(8,5,5),col=1,lwd=1.5)
lines(pred8b,lag=15,lwd=1.5,lty=2,col=1)
lines(pred8c,lag=15,lwd=1.5,lty=4,col=1)
lines(pred8d,lag=15,lwd=1.5,lty=5,col=1)
legend("topleft",c("Knots: 13.3","Knots: 20","Knots: 26.6","Knots: 13.3,26.6"),
  lty=c(1:2,4:5),lwd=1.5,col=1,cex=0.9,bty="n",inset=0.02)

dev.off()

#
