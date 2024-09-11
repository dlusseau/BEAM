##### ICES WGBYC 2023
setwd("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/")

library(glmmTMB)
library(metafor)
library(emmeans)
library(ggeffects)


monitor.effort<-read.csv("D2_Bycatch_monitoring_effort_Last5Years_v2.csv",header=T)

event<-read.csv("D3_Bycatch_event_Last5Years_v2.csv",header=T)

region<-unique(event$AreaCode)
species<-unique(event$Species)
ML4<-unique(monitor.effort$MetierL4)



for (i in 1:length(region)) {
for (j in 1:length(species) {
for (k in 1:length(ML4)) {


data.temp<-subset(table,Species==species[j] & AreaCode==region[i] & MetierL4==ML4[k])

if (dim(data.temp)[1]<2) {
BPUE<-data.temp$TotalBycatch/data.temp$DaysAtSeaOb
BPUE.LCL

}

if (dim(data.temp)[1]<5) {

bpue.mod<-glmmTMB(TotalBycatch~1,offset=log(DaysAtSeaOb),family=nbinom2,data=data.temp)

bpue.r<-as.data.frame(emmeans(bpue.mod,~1,type="response",offset=log(1))) # BPUE

#### here we should simply predict for the total effort DaS ?

BPUE<-bpue.r$response
BPUE.LCL<-BPUE.r$asymp.LCL
BPUE.UCL<-BPUE.r$asymp.UCL #asymptotic 95% CI

#place in estimate table

} else {

glmmTMB(number1+(1|nation)+(1|year)+(1|L5)+(1|protocol)+(1|sampling))



} #end else small sample size

} # end k
} #end j
} # end i

data.temp<-data.frame(TotalBycatch=sample.int(10,5),DaysAtSeaOb=sample.int(1000,5))

data.temp.re<-data.frame(TotalBycatch=sample.int(10,30,replace=T),DaysAtSeaOb=sample.int(1000,30,replace=T),nation=rep(c("a","b","c"),10))

bpue<-glmmTMB(TotalBycatch~1,offset=log(DaysAtSeaOb),family=nbinom2,data=data.temp.re)

bpue<-glmmTMB(TotalBycatch~1+(1|nation),offset=DaysAtSeaOb,family=nbinom2,data=data.temp.re)
bpue2<-glmmTMB(TotalBycatch~1+(1|nation),offset=DaysAtSeaOb,family=nbinom2,data=data.temp.re)




library(ggeffects)

tau.est <- sqrt(as.numeric(summary(bpue)$varcor)) ### from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7398826/
I2<-tau.est^2/(+tau.est^^2)
#Q is chi2 distributed with K-1 df
#nbinom2 is log link

  wi    <- 1/vi
  vt <- 1/mean(wi)
  tau2   <- glmmTMB::VarCorr(bpue)[[1]][[1]][[1]]
  I2  <- 100 * tau2 / (vt + tau2)
  
  p=1
  beta   <- cbind(glmmTMB::fixef(bpue)$cond[seq_len(p)])
                  vb     <- as.matrix(vcov(bpue)$cond)[seq_len(p),seq_len(p),drop=FALSE]
                  tau2   <- glmmTMB::VarCorr(bpue)[[1]][[1]][[1]]

  
escalc(measure="PR", data=data.temp.re,xi=TotalBycatch, ti=DaysAtSeaOb,addvi=addvi)


X.fit + study + (group - 1 | study), offset=dat.off, family=dat.fam, verbose=verbose, data=NULL, control=do.call(glmmTMB::glmmTMBControl, glmerCtrl)), silent=!verbose)