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
# PREDICTIONS
####################################################################


################################################################################
# PREDICTION FOR SPECIFIC YEARS GIVEN AN EXPOSURE HISTORY (TABLE 3)

# CREATE THE LAYOUT OF THE RESULTS
tab3 <- matrix(NA,5,12)

# CREATE THE MATRIX OF EXPOSURE HISTORIES USED FOR PREDICTION
hist <- rbind(exphist(rep(20,10),time=10,lag=c(2,40)),
  exphist(rep(100,10),time=10,lag=c(2,40)),
  exphist(rep(20,20),time=20,lag=c(2,40)),
  exphist(rep(c(20,0),c(10,10)),time=20,lag=c(2,40)),
  exphist(rep(c(20,0),c(10,30)),time=40,lag=c(2,40))
)

# DEFINE THE MODELS USED FOR PREDICTION
modpred <- c(1,4,5,8)

# RUN THE MODELS, PREDICT AND FILL THE TABLE
for(i in seq(modpred)) {
  # CREATE THE CROSS-BASIS
  cbx <- crossbasis(Qx,lag=c(2,40),argvar=varlist[[modpred[i]]],
    arglag=laglist[[modpred[i]]])
  # RUN THE MODEL
  model <- coxph(Surv(agest,ageexit,ind)~cbx+cbz+caltime,mainsub)
  pred <- crosspred(cbx,model,at=hist,cen=0)
  # FILL THE TABLE
  tab3[,1:3+(i-1)*3] <- with(pred,cbind(allRRfit,allRRlow,allRRhigh))
}
  
# FORMAT THE TABLE
tab3 <- formatC(tab3,digits=2,format="f")
table3 <- matrix(NA,nrow(tab3),length(modpred))
for(i in seq(modpred)) table3[,i] <- paste(tab3[,1+(i-1)*3]," (",
  paste(tab3[,2+(i-1)*3],tab3[,3+(i-1)*3],sep="--"),")",sep="")
colnames(table3) <- paste("Model",modpred)
rownames(table3) <- paste("Scenario",seq(nrow(table3)))

table3
xtable(table3,align="lcccc")

################################################################################
# PREDICTION ALONG TIME (FIGURE 5 IN THE MANUSCRIPT)

# MATRIX OF EXPOSURE HISTORIES COMPUTED FOR YEAR 1-60 ALONG THE EXPOSURE PROFILE
profile <- rep(c(20,0),c(15,45))
dynhist <- exphist(profile,lag=c(2,40))

# PREDICTIONS (MODELS ALREADY FITTED)
predhist4 <- crosspred(cbx4,model4,at=dynhist,cen=0)
predhist5 <- crosspred(cbx5,model5,at=dynhist,cen=0)
predhist8 <- crosspred(cbx8,model8,at=dynhist,cen=0)

# PLOT 
pdf("pred.pdf",height=4.1,width=6)
layout(1)
par(mar=c(4,4,0.5,0.5))
plot(predhist8,"overall",xlab="Years",ylab="HR",ylim=c(0,4.5),
  lab=c(8,5,5),col=1,lwd=1.5)
lines(predhist5,"overall",lwd=1.5,lty=2,col=1)
lines(predhist4,"overall",lwd=1.5,lty=4,col=1)
legend("topright",paste("Model",c(8,5,4)),lty=c(1:2,4),lwd=1.5,col=1,cex=0.9,
  bty="n",inset=0.1)
segments(15,-1,15,4)
dev.off()

# PREDICTIONS AT YEAR 20
with(predhist4,cbind(allRRfit,allRRlow,allRRhigh))[20,]
with(predhist5,cbind(allRRfit,allRRlow,allRRhigh))[20,]
with(predhist8,cbind(allRRfit,allRRlow,allRRhigh))[20,]

#
