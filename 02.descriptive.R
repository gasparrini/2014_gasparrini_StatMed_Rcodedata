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
# DESCRIPTIVE STATISTICS (TABLE 1 IN THE MANUSCRIPT)
####################################################################

table1 <- matrix(NA,18,10)
colnames(table1) <- c("Full cohort",rep(" ",4),"Lung cancer cases",rep(" ",4))
rownames(table1) <- c(
  "Subjects",
  "Deaths (%)",
  "Ever smokers (%)",
  "Age at enter",
  "Follou-up time  (years)",
  "Exposure to radon",
  "Exposure period (years)",
  "Total cumulative exposure  (WLM/year)",
  "Yearly exposure  (WLM/year)",
  "- All",
  "- Lag 0--9",
  "- Lag 10--19",
  "- Lag 20--29",
  "- Lag 30--40",
  "Smoking",
  "Exposure period (years) ",
  "Total cumulative exposure (pack x 100)",
  "Yearly exposure (pack x 100)"
)

################################################################################

fvtab <- function(v,ind) c(median(v),min(v),quantile(v,0.25),quantile(v,0.75),
  max(v),median(v[ind]),min(v[ind]),quantile(v[ind],0.25),quantile(v[ind],0.75),
  max(v[ind]))
fmtab <- function(m,ind) c(median(m[m>0]),min(m[m>0]),quantile(m[m>0],0.25),
  quantile(m[m>0],0.75),max(m[m>0]),median(m[ind,][m[ind,]>0]),
  min(m[ind,][m[ind,]>0]),quantile(m[ind,][m[ind,]>0],0.25),
  quantile(m[ind,][m[ind,]>0],0.75),max(m[ind,][m[ind,]>0]))

vcase <- uminers$ind==1
mcase <- mainsub$ind==1
vsmk <- uminers$totsmk>0

################################################################################

table1[1,c(1,2,6,7)] <- c(nrow(uminers),100,sum(uminers$ind),100)
table1[2,c(1,2,6,7)] <- with(uminers,c(sum(status==1),sum(status==1)/
  table1[1,1]*100,sum((status==1)[vcase]),sum((status==1)[vcase])/
  table1[1,6]*100))
table1[3,c(1,2,6,7)] <- with(uminers,c(sum(vsmk),sum(vsmk)/table1[1,1]*100,
  sum(vsmk[vcase]),sum(vsmk[vcase])/table1[1,6]*100))

table1[4,] <- fvtab(uminers$agest,vcase)
table1[5,] <- fvtab(uminers$ageexit-uminers$agest,vcase)

table1[7,] <- fvtab(uminers$rendage-uminers$rdnstar,vcase)
table1[8,] <- fvtab(uminers$totrdn,vcase)

table1[10,] <- fmtab(Qx,mcase)
table1[11,] <- fmtab(Qx[,1:8],mcase)
table1[12,] <- fmtab(Qx[,9:18],mcase)
table1[13,] <- fmtab(Qx[,19:28],mcase)
table1[14,] <- fmtab(Qx[,29:39],mcase)

table1[16,] <- with(uminers[vsmk,],fvtab(sendage-smkstar,vcase[vsmk]))
table1[17,] <- fvtab(uminers$totsmk[vsmk],vcase[vsmk])
table1[18,] <- fmtab(sexp[vsmk,],vcase[vsmk])

################################################################################

library(xtable)
xtable(table1,digits=1,display=rep("f",11),align="l|rrrrr|rrrrr")

#
