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
# MODELS WITH NON-LINEAR EXPOSURE-RESPONSE RELATIONSHIP (DLNM)
####################################################################

################################################################################
# BUILD THE SECOND PART OF TABLE 2 IN THE MANUSCRIPT

tab2b <- matrix(NA,4,5)
colnames(tab2b) <- c("f(x)","w(ell)","df","AIC","BIC")
rownames(tab2b) <- paste("Model",5:8)

# RUN THE MODELS
for(i in 5:8) {
  # DEFINE THE FUNCTION IN EACH OF THE TWO DIMENSIONS
  tab2b[i-4,1:2] <- c(varlist[[i]]$fun,laglist[[i]]$fun)
  # CREATE THE CROSS-BASIS
  cbx <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[i]],arglag=laglist[[i]])
  tab2b[i-4,3] <- ncol(cbx)
  # RUN THE MODEL
  model <- coxph(Surv(agest,ageexit,ind)~cbx+cbz+caltime,mainsub)
  # SAVE THE AIC VALUE
  tab2b[i-4,4] <- formatC(AIC(model),digits=1,format="f")
  tab2b[i-4,5] <- formatC(BIC(model),digits=1,format="f")
}
tab2b

table2 <- rbind(tab2a,tab2b)
table2
xtable(table2,align="lllccc")

################################################################################
# RE-FIT THE MODELS

cbx5 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[5]],arglag=laglist[[5]])
model5 <- coxph(Surv(agest,ageexit,ind)~cbx5+cbz+caltime,mainsub)
pred5 <- crosspred(cbx5,model5,at=0:25*10,cen=0)

cbx6 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[6]],arglag=laglist[[6]])
model6 <- coxph(Surv(agest,ageexit,ind)~cbx6+cbz+caltime,mainsub)
pred6 <- crosspred(cbx6,model6,at=0:25*10,cen=0)

cbx7 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[7]],arglag=laglist[[7]])
model7 <- coxph(Surv(agest,ageexit,ind)~cbx7+cbz+caltime,mainsub)
pred7 <- crosspred(cbx7,model7,at=0:25*10,cen=0)

cbx8 <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[8]],arglag=laglist[[8]])
model8 <- coxph(Surv(agest,ageexit,ind)~cbx8+cbz+caltime,mainsub)
pred8 <- crosspred(cbx8,model8,at=0:25*10,cen=0)

################################################################################
# LR TESTS

# BETWEEN MODEL 8 AND MODEL 4 (NON-LINEARITY)
flrtest(model8,model4)

# BETWEEN MODEL 8 AND MODEL 5 (CONSTANT EFFECT)
flrtest(model8,model5)

# BETWEEN MODEL 8 AND MODEL 7 (NON-ZERO EFFECT AT LAG 2)
flrtest(model8,model7)

################################################################################
# 3D AND SLICES PLOTS FROM MODEL8 AND MODEL6 (FIGURE 2)

# PLOT
pdf("dlnm.pdf",height=7,width=9.5)

layout(matrix(1:4,ncol=2,byrow=TRUE),height=c(1,0.7))
par(mar=c(2,2,1.5,0.5))

d3 <- plot(pred8,theta=230,cex.axis=0.7,cex.lab=0.8,nticks=9,ltheta=230,lphi=45,
  border=grey(0.3),xlab="WLM/year",ylab="Lag (years)",zlab="HR",col=grey(0.99),
  zlim=c(0.95,1.3))
mtext("Model 8",cex=0.8)
lines (trans3d(x=100,y=10:40,z=pred8$matRRfit["100",-c(1:8)],pmat=d3),
  col=1,lwd=3)
lines (trans3d(x=0:25*10,y=15,z=pred8$matRRfit[,"lag15"],pmat=d3),
  col=1,lwd=3)

d3 <- plot(pred6,theta=230,cex.axis=0.7,cex.lab=0.8,nticks=9,ltheta=230,lphi=45,
  border=grey(0.3),xlab="WLM/year",ylab="Lag (years)",zlab="HR",col=grey(0.99),
  zlim=c(0.95,1.3))
mtext("Model 6",cex=0.8)
lines (trans3d(x=100,y=10:40,z=pred6$matRRfit["100",-c(1:8)],pmat=d3),
  col=1,lwd=3)
lines (trans3d(x=0:25*10,y=15,z=pred6$matRRfit[,"lag15"],pmat=d3),
  col=1,lwd=3)

par(mar=c(4,4,0.5,0.5))

plot(pred8,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5)
lines(pred6,var=100,lwd=1.5,lty=2,col=1)
legend("topright",c("Model 8","Model 6"),lty=1:2,lwd=1.5,col=1,cex=0.9,
  bty="n",inset=0.1)

plot(pred8,lag=15,xlab="WLM/year",ylab="HR at lag 15",ylim=c(0.9,1.3),
  lab=c(8,5,5),col=1,lwd=1.5)
lines(pred6,lag=15,lwd=1.5,lty=2,col=1)
lines(pred4,lag=15,lwd=1.5,lty=4,col=1)
legend("topleft",paste("Model",c(8,6,4)),lty=c(1:2,4),lwd=1.5,col=1,cex=0.9,
  bty="n",inset=0.02)

dev.off()

################################################################################
# EXPOSURE AND LAG-SPECIFIC PLOTS FROM MODEL8 (FIGURE 3)

pdf("dlnm2.pdf",height=4.1,width=12)

layout(matrix(1:2,ncol=2,byrow=TRUE))
par(mar=c(4,4,0.5,0.5))

plot(pred8,var=100,xlim=c(0,40),xlab="Lag (years)",ylab="HR for 100 WLM/year",
  ylim=c(0.9,1.3),lab=c(8,5,5),col=1,lwd=1.5,ci="n")
lines(pred8,var=20,lwd=1.5,lty=2,col=1)
lines(pred8,var=50,lwd=1.5,lty=3,col=1)
lines(pred8,var=200,lwd=1.5,lty=4,col=1)
legend("topright",paste(c(20,50,100,200),"WLM/year"),lty=c(2,3,1,4),lwd=1.5,
  col=1,cex=0.9,bty="n",inset=0.1)

plot(pred8,lag=15,xlab="WLM/year",ylab="HR at lag 15",ylim=c(0.9,1.3),
  lab=c(8,5,5),col=1,lwd=1.5,ci="n")
lines(pred8,lag=5,lwd=1.5,lty=2,col=1)
lines(pred8,lag=10,lwd=1.5,lty=3,col=1)
lines(pred8,lag=25,lwd=1.5,lty=4,col=1)
legend("topleft",paste("Lag",c(5,10,15,25)),lty=c(2,3,1,4),lwd=1.5,col=1,cex=0.9,
  bty="n",inset=0.03)

dev.off()

################################################################################
# RISK CONTRIBUTION FOR 100 WLM AT LAG 10:12 FROM MODEL8

with(pred8,matRRfit["100",paste("lag",10:12,sep="")])
with(pred8,matRRlow["100",paste("lag",10:12,sep="")])
with(pred8,matRRhigh["100",paste("lag",10:12,sep="")])

# PEAK RISK AT VARIOUS EXPOSURE VALUES FROM MODEL 8
with(pred8,apply(matRRfit,1,function(x) colnames(matRRfit)[which.max(x)]))

#
