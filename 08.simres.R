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

# LOAD THE DATA WITH THE RESULTS OF THE SIMULATIONS
load("simul_400.RData")

################################################################################
# DISTRIBUTION OF THE EXPOSURE OCCURRENCES AND EXAMPLES OF EXPOSURE HISTORIES

summary(as.numeric(expsim))

pdf("expdistr.pdf",height=7,width=10)

layout(matrix(1:4,ncol=2,byrow=TRUE))
hist(as.numeric(expsim),col=2,xlab="Exposure level",
  main="Distribution of exposure occurrences")
par(mex=0.8)
for(i in 1:3*5) plot(expsim[i,],type="s",xlab="Time",ylab="Exposure",
  ylim=c(0,11),col=4,cex=0.6,main=paste("Subject",i))

dev.off()

################################################################################
# PLOT OF THE FUNCTIONS USED TO SIMULATE THE EFFECT SURFACE

pdf("simulfun.pdf",height=5,width=14)

layout(matrix(1:2,ncol=2))
par(mar=c(5,4,1,1))

plot(0:10,type="n",xlab="x",ylab=expression(f(x)),
  xlim=c(-0.5,10.5),ylim=c(-0.005,0.8))
curve(flin,0,10,col=2,lty=2,lwd=2,add=TRUE)
curve(fplat,0,10,col=3,lty=4,lwd=2,add=TRUE)
curve(fexp,0,10,col=4,lty=5,lwd=2,add=TRUE)
legend("topleft",c("Linear","Plateau","Exponential"),lty=c(2,4,5),lwd=2,
  col=2:4,cex=0.8,bty="n")

# PLOT OF THE FUNCTIONS USED TO SIMULATE THE LAG STRUCTURE
plot(0:40,0:40/40,type="n",xlab="Lag",ylab=expression(w(l)),
  xlim=c(-0.5,40.5),ylim=c(-0.05,0.8))
curve(fconst,0,40,col=2,lty=2,lwd=2,add=TRUE)
curve(fdecay,0,40,col=3,lty=4,lwd=2,add=TRUE)
curve(fpeak,0,40,col=4,lty=5,lwd=2,add=TRUE)
legend("topright",c("Constant","Decay","Peak"),lty=c(2,4,5),lwd=2,
  col=2:4,cex=0.8,bty="n")

dev.off()

################################################################################
# GRAPHICAL REPRESENTATION OF THE SIMULATED EFFECT SURFACES

pdf("simuleff.pdf",height=10,width=10)

layout(matrix(1:9,ncol=3,byrow=TRUE))
par(mar=c(1,1,3,1))

for(j in seq(effsimlist)) {
  do.call(persp,modifyList(arg3D,list(z=exp(effsimlist[[j]]))))
  title(rownames(combsim)[j])
}

dev.off()


################################################################################
# TABLE 4 IN MANUSCRIPT

tab4 <- formatC(cbind(abs(colSums(aicpred-trueeff))/colSums(trueeff),
  abs(colSums(bicpred-trueeff))/colSums(trueeff),
  colMeans(aiccov),colMeans(biccov),
  colSums((aicpred-trueeff)^2)/colSums(trueeff),
  colSums((bicpred-trueeff)^2)/colSums(trueeff)),
  format="f",digits=2)
table4 <- rbind(rep(c("AIC","BIC"),3),tab4)
rownames(table4) <- c("Model",rownames(combsim))

table4
library(xtable)
xtable(table4,align="lcccccc")


################################################################################
# TABLE 5 IN MANUSCRIPT

tab5 <- formatC(cbind(colMeans(aicdfvar),colMeans(bicdfvar),colMeans(aicdflag),
  colMeans(bicdflag),colMeans(aicdfvar!=1),colMeans(bicdfvar!=1),
  colMeans(aicdflag!=1),colMeans(bicdflag!=1)),
  format="f",digits=2)
colnames(tab5) <- rep(c("df","test"),each=4)

table5 <- rbind(rep(c("f(x)","w(lag)"),2,each=2),rep(c("AIC","BIC"),4),tab5)
rownames(table5) <- c("","Model",rownames(combsim))

table5
xtable(table5,align="lcccccccc")


################################################################################
# GRAPH OF SIMULATION RESULTS (FIGURE 5)

var3d <- rep(c("5","5","8"),3)
var3d2 <- rep(c(5,5,8),3)
lag3d <- rep(paste("lag",c(20,4,12),sep=""),3)
lag3d2 <- rep(c(20,4,12),3)

ylimhigh <- vector("list",9)
ylimhigh[[1]] <- rep(1.3,4)
ylimhigh[[5]] <- rep(1.6,4)
ylimhigh[[9]] <- rep(c(1.5,1.6),c(2,2))

pdf("simulmain.pdf",height=12,width=9.5)

layout(matrix(1:15,5),heights=c(1,rep(0.8,4)))
for(j in c(1,5,9)) {
  
  par(mar=c(1,1,3,1))
  d3 <- do.call(persp,modifyList(arg3D,list(z=exp(effsimlist[[j]]))))
  title(rownames(combsim)[j])
  lines (trans3d(x=var3d2[j],y=0:20*2,z=exp(effsimlist[[j]])[var3d[j],],pmat=d3),
    col=1,lwd=2)
  lines (trans3d(x=0:40/4,y=lag3d2[j],z=exp(effsimlist[[j]])[,lag3d[j]],pmat=d3),
    col=1,lwd=2)
  
  par(mar=c(5,4,0.5,1))
  plot(0:20*2,exp(effsimlist[[j]][var3d[j],]),type="n",xlab="Lag",ylab="HR",
    ylim=c(0.9,ylimhigh[[j]][1]),frame.plot=FALSE)
  for(m in aictotsample[[j]]) lines(0:20*2,exp(m[var3d[j],]),
    lwd=1.5,col=grey(0.8))
  abline(h=1)
  lines(0:20*2,exp(aictotpred[[j]][var3d[j],]),lty=5,lwd=1.5)
  lines(0:20*2,exp(effsimlist[[j]][var3d[j],]),lwd=1.5)
  legend("top",c("True","AIC avg","AIC samples"),lty=c(1,5,1),cex=0.7,
    col=c(1,1,grey(0.7)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
  
  par(mar=c(5,4,0.5,1))
  plot(0:20*2,exp(effsimlist[[j]][var3d[j],]),type="n",xlab="Lag",ylab="HR",
    ylim=c(0.9,ylimhigh[[j]][2]),frame.plot=FALSE)
  for(m in bictotsample[[j]]) lines(0:20*2,exp(m[var3d[j],]),
    lwd=1.5,col=grey(0.8))
  abline(h=1)
  lines(0:20*2,exp(bictotpred[[j]][var3d[j],]),lty=5,lwd=1.5)  
  lines(0:20*2,exp(effsimlist[[j]][var3d[j],]),lwd=1.5)
  legend("top",c("True","BIC avg","BIC samples"),lty=c(1,5,1),cex=0.7,
    col=c(1,1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)
  
  plot(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),type="n",xlab="Exposure",
    ylab="HR",ylim=c(0.9,ylimhigh[[j]][3]),frame.plot=FALSE)
  for( m in aictotsample[[j]]) lines(seq(0,10,0.25),exp(m[,lag3d[j]]),
    lwd=1.5,col=grey(0.8))  
  abline(h=1)
  lines(seq(0,10,0.25),exp(aictotpred[[j]][,lag3d[j]]),lty=5,lwd=1.5)
  lines(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),lwd=1.5)
  legend("top",c("True","AIC avg","AIC samples"),lty=c(1,5,1),cex=0.7,
    col=c(1,1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)

  plot(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),type="n",xlab="Exposure",
    ylab="HR",ylim=c(0.9,ylimhigh[[j]][4]),frame.plot=FALSE)
  for( m in bictotsample[[j]]) lines(seq(0,10,0.25),exp(m[,lag3d[j]]),
    lwd=1.5,col=grey(0.8))
  abline(h=1)
  lines(seq(0,10,0.25),exp(bictotpred[[j]][,lag3d[j]]),lty=5,lwd=1.5)  
  lines(seq(0,10,0.25),exp(effsimlist[[j]][,lag3d[j]]),lwd=1.5)
  legend("top",c("True","BIC avg","BIC samples"),lty=c(1,5,1),cex=0.7,
    col=c(1,1,grey(0.7)),horiz=T,bty="n",x.intersp=0.2,inset=0.03)
}

dev.off()

################################################################################
# GRAPH OF COVERAGE IN THE BI-DIMENSIONAL SURFACE (FIGURE 6)

col <- colorRampPalette(c("black","white"))

for(j in c(1,5,9)) {
  
  pdf(paste("aiccov",combsim[j,1],combsim[j,2],".pdf",sep=""),width=7,height=5)
  
  filled.contour(seq(0,10,0.25),0:20*2,aictotcov[[j]],zlim=c(0,1),color.palette=col,
    xlab="Exposure",ylab="Lag",key.axes=axis(4,at=0:10/10),nlevels=10)
  title(paste(rownames(combsim)[j],"- AIC"))
  
  dev.off()
  
  pdf(paste("biccov",combsim[j,1],combsim[j,2],".pdf",sep=""),width=7,height=5)
  
  filled.contour(seq(0,10,0.25),0:20*2,bictotcov[[j]],zlim=c(0,1),color.palette=col,
    xlab="Exposure",ylab="Lag",key.axes=axis(4,at=0:10/10),nlevels=10)
  title(paste(rownames(combsim)[j],"- BIC"))
  
  dev.off()
}

#
