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
# MODELS WITH LINEAR EXPOSURE-RESPONSE RELATIONSHIPS (DLM)
####################################################################

################################################################################
# BUILD THE FIRST PART OF TABLE 2 IN THE MANUSCRIPT

tab2a <- matrix(NA,4,5)
colnames(tab2a) <- c("f(x)","w(ell)","df","AIC","BIC")
rownames(tab2a) <- paste("Model",1:4)

# RUN THE MODELS
for(i in 1:4) {
  # DEFINE THE FUNCTIONS f(x) and w(lag)
  tab2a[i,1:2] <- c(varlist[[i]]$fun,laglist[[i]]$fun)
  # CREATE THE CROSS-BASIS
  cbx <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[i]],arglag=laglist[[i]])
  tab2a[i,3] <- ncol(cbx)
  # RUN THE MODEL
  model <- coxph(Surv(agest,ageexit,ind)~cbx+cbz+caltime,mainsub)
  # SAVE THE AIC AND BIC VALUES
  tab2a[i,4] <- formatC(AIC(model),digits=1,format="f")
  tab2a[i,5] <- formatC(BIC(model),digits=1,format="f")
}
tab2a

################################################################################
# RE-FIT THE MODELS

cbx1 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[1]],arglag=laglist[[1]])
model1 <- coxph(Surv(agest,ageexit,ind)~cbx1+cbz+caltime,mainsub)
pred1 <- crosspred(cbx1,model1,at=0:25*10)

cbx2 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[2]],arglag=laglist[[2]])
model2 <- coxph(Surv(agest,ageexit,ind)~cbx2+cbz+caltime,mainsub)
pred2 <- crosspred(cbx2,model2,at=0:25*10)

cbx3 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[3]],arglag=laglist[[3]])
model3 <- coxph(Surv(agest,ageexit,ind)~cbx3+cbz+caltime,mainsub)
pred3 <- crosspred(cbx3,model3,at=0:25*10)

cbx4 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[4]],arglag=laglist[[4]])
model4 <- coxph(Surv(agest,ageexit,ind)~cbx4+cbz+caltime,mainsub)
pred4 <- crosspred(cbx4,model4,at=0:25*10)

################################################################################
# PLOT FROM MODEL4, MODEL2 AND MODEL1 (FIGURE 1)

# PLOT
pdf("dlm.pdf",height=4.1,width=6)
par(mar=c(4,4,0.5,0.5))
plot(pred4,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  lab=c(8,5,5),col=1,lwd=1.5)
lines(pred2,var=100,lwd=1.5,lty=2,col=1)
lines(pred1,var=100,lwd=1.5,lty=4,col=1)
legend("top",c("Model 4","Model 2","Model 1"),lty=c(1:2,4),lwd=1.5,col=1,
  cex=0.9,bty="n",inset=0.1)
dev.off()

################################################################################
# RISK CONTRIBUTION FOR 100 WLM AT LAG 10:12 FROM MODEL4 AND MODEL1

with(pred4,matRRfit["100",paste("lag",10:12,sep="")])
with(pred4,matRRlow["100",paste("lag",10:12,sep="")])
with(pred4,matRRhigh["100",paste("lag",10:12,sep="")])

with(pred1,matRRfit["100",paste("lag",10:12,sep="")])
with(pred1,matRRlow["100",paste("lag",10:12,sep="")])
with(pred1,matRRhigh["100",paste("lag",10:12,sep="")])

################################################################################
# LR TESTS 

# BETWEEN MODEL 3 AND MODEL 1 (NON-LINEARITY)
flrtest(model3,model1)

# BETWEEN MODEL 3 AND MODEL 4 (NON-ZERO EFFECT AT LAG 2)
flrtest(model3,model4)

#
