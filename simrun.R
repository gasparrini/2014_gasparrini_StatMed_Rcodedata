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
# START THE ITERATIONS
# 3 LEVELS:
#   - TYPE OF SURFACE (COMBINATIONS OF FUNCTIONS USED FOR SIMULATE EFFECTS)
#   - RANDOM GENERATED DATA
#   - TYPE OF MODEL FOR ESTIMATION (COMBINATIONS OF FUNCTIONS IN CROSS-BASIS)

# NB: EACH ITERATION OF THE SIMULATION ROUTINE TAKES FROM 2 TO 15 MINS FOR 
#   SIMULATED DATA SETS WITH 200 AND 800 SUBJECTS, RESPECTIVELY
# THE SIMULATION OF 500 ITERATIONS TAKES FROM 17 HOURS TO MORE THAN 5 DAYS
# TESTED IN A PC WITH 2.66GHz AND 4GB OF RAM

time <- proc.time()
# LOOP ACROSS SIMULATED EXPOSURE-LAG-RESPONSES
for(j in seq(nrow(combsim))) {
 
  # PRINT
  cat("\n\n Scenario",j,"\n")

  # TEMPORARY OBJECTS TO STORE THE RESULTS FROM EACH ITERATION BELOW
  aictotpredtemp <- bictotpredtemp <- aictotcovtemp <- bictotcovtemp <- 0
  
  # FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
  fcumeff <- function(hist,lag,f1,f2) sum(do.call(f1,list(hist)) * 
      do.call(f2,list(lag[1]:lag[2])))
  # FUNCTION TO COMPUTE CUMULATIVE EFFECTS AT EACH POINT OF AN EXPOSURE PROFILE
  fcumeffexp <- function(exp,lag,f1,f2) apply(exphist(exp,lag=lag),1,fcumeff,
    lag,f1,f2)
    
  # COMPUTE THE CUMULATIVE EFFECT AT EACH TIME FOR EACH SUBJECT
  cumeffmat <- t(apply(expsim,1,fcumeffexp,lag=c(0,40),
    f1=combsim[j,1],f2=combsim[j,2]))

################################################################################
# LOOP ACROSS RANDOMLY SIMULATED DATA

  for(i in seq(nsim)) {
    
    # PRINT
    cat(i,"")
 
    # SET THE SEED
    seed <- 100000 + i
    set.seed(seed)
    
    # SAMPLE THE SURVIVAL TIMES AND EVENTS
    # - USE CUMULATIVE EFFECT AS TIME-DEPENDENT VARIABLE, WITH COEF=1
    # - SET CENSORING PROBABILITY TO 20%
    data <- permalgorithm(nsub,100,Xmat=as.numeric(t(cumeffmat)),
      censorRandom=runif(nsub,1,100*2),betas=1)
    # RESTRICT DATA TO EVENT/CENSOR TIMES
    data <- data[cumsum(tapply(data$Id,data$Id,length)),c("Event","Fup")]
    names(data) <- c("event","exit")
    # ADD AN ID
    data <- cbind(id=seq(nrow(expsim)),as.data.frame(data))
    
    # SPLIT THE DATA
    ftime <- sort(unique(data$exit[data$event==1]))
    dataspl <- survSplit(Surv(exit,event)~., data, cut=ftime, start="enter")
    dataspl <- dataspl[order(dataspl$id),]
    
    # CREATE THE MATRIX Q OF LAGGED EXPOSURES USING THE FUNCTION flaghist
    Q <- do.call(rbind, lapply(seq(nrow(dataspl)),
      function(i) exphist(expsim[dataspl$id[i],],dataspl$exit[i],c(0,40))))
    
    # LOOP ACROSS MODELS WITH DIFFERENT FUNCTIONS
    for(k in seq(length(svarlist))) {
    
      # DEFINE THE CROSS-BASIS
      cb <- crossbasis(Q,lag=c(0,40),argvar=svarlist[[k]],arglag=slaglist[[k]])
      
      # RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
      modellist[[k]] <- coxph(Surv(enter,exit,event)~cb,dataspl,
        y=FALSE,ties="efron")
    }
    
    # DETERMINE THE BEST MODEL FOR AIC AND BIC
    bestaic <- which.min(sapply(modellist,AIC))
    bestbic <- which.min(sapply(modellist,BIC))

################################################################################
# STORE RESULTS FOR BEST FITTING MODELS

    # STORE SUM OF PREDICTIONS FOR THE WHOLE SURFACE
    cb <- crossbasis(Q,lag=c(0,40),argvar=svarlist[[bestaic]],
      arglag=slaglist[[bestaic]])
    cpaic1 <- crosspred(cb,modellist[[bestaic]],from=0,to=10,by=0.25,
      bylag=2,cen=0)
    cb <- crossbasis(Q,lag=c(0,40),argvar=svarlist[[bestbic]],
      arglag=slaglist[[bestbic]])
    cpbic1 <- crosspred(cb,modellist[[bestbic]],from=0,to=10,by=0.25,
      bylag=2,cen=0)
    aictotpredtemp <- aictotpredtemp + cpaic1$matfit
    bictotpredtemp <- bictotpredtemp + cpbic1$matfit
    
    # STORE COVERAGE FOR THE WHOLE SURFACE
    aictotcovtemp <- aictotcovtemp + (effsimlist[[j]] >= cpaic1$matfit-qn*
      cpaic1$matse & effsimlist[[j]] <= cpaic1$matfit+qn*cpaic1$matse)
    bictotcovtemp <- bictotcovtemp + (effsimlist[[j]] >= cpbic1$matfit-qn*
      cpbic1$matse & effsimlist[[j]] <= cpbic1$matfit+qn*cpbic1$matse)    
    
    # STORE SAMPLES
    if(i<=nsample) {
      aictotsample[[j]] <- c(aictotsample[[j]],list(cpaic1$matfit))
      bictotsample[[j]] <- c(bictotsample[[j]],list(cpbic1$matfit))
    }
    
    # SAMPLE RANDOM EXPOSURE HISTORY USED FOR ASSESSMENT
    rsub <- sample(seq(nsub),1)
    rtime <- sample(41:100,1)
    history <- exphist(expsim[rsub,],rtime,c(0,40))
    
    # TRUE VALUE
    trueeff[i,j] <- cumeffmat[rsub,rtime]
    
    # STORE PREDICTED VALUE
    cb <- crossbasis(Q,lag=c(0,40),argvar=svarlist[[bestaic]],
      arglag=slaglist[[bestaic]])
    cpaic2 <- crosspred(cb,modellist[[bestaic]],at=history,cen=0)
    cb <- crossbasis(Q,lag=c(0,40),argvar=svarlist[[bestbic]],
      arglag=slaglist[[bestbic]])
    cpbic2 <- crosspred(cb,modellist[[bestbic]],at=history,cen=0)
    aicpred[i,j] <- cpaic2$allfit
    bicpred[i,j] <- cpbic2$allfit
    
    # STORE COVERAGE
    aiccov[i,j] <- trueeff[i,j] >= cpaic2$allfit-qn*cpaic2$allse & 
      trueeff[i,j] <= cpaic2$allfit+qn*cpaic2$allse
    biccov[i,j] <- trueeff[i,j] >= cpbic2$allfit-qn*cpbic2$allse & 
      trueeff[i,j] <= cpbic2$allfit+qn*cpbic2$allse
    
    # STORE DF
    
    aicdfvar[i,j] <- rep(c(1,2,3,3,3,4),6)[bestaic]
    aicdflag[i,j] <- rep(c(1,2,3,3,3,4),each=6)[bestaic]
    bicdfvar[i,j] <- rep(c(1,2,3,3,3,4),6)[bestbic]
    bicdflag[i,j] <- rep(c(1,2,3,3,3,4),each=6)[bestbic]
    
    # STORE NUMBER OF EVENTS
    navg[i,j] <- sum(data$event==1)
  }

################################################################################
  
  # COMPUTE AVERAGE PREDICTION AND COVERAGE FOR THE WHOLE SURFACE
  aictotpred[[j]] <- aictotpredtemp/nsim
  bictotpred[[j]] <- bictotpredtemp/nsim
  aictotcov[[j]] <- aictotcovtemp/nsim
  bictotcov[[j]] <- bictotcovtemp/nsim
}
proc.time()-time

rm(modellist,Q,dataspl,cb)

#
