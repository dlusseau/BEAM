
calc_bpue <- function(Ecoregion, ml4, Species, min_re_obs = 2,ices=THEDATA) {

### I am thinking we put a line here tolower all colnames to avoid inter-annual changes
### and perhaps use amatch to 'find' the right column

  needle <- data.table(ecoregion = Ecoregion, metierL4 = ml4, species = Species)
  #  monitor <- data$monitor[needle, on = .(ecoregion, ml4), nomatch = 0]
 #   bycatch <- data$bycatch[needle, on = .(ecoregion, ml4, species), nomatch = 0]
#dat <- THEDATA[Ecoregion==ecoregion&MetierL4==ml4&Species==species,]
dat <- ices[needle, on = .(ecoregion, metierL4,species), nomatch = 0]
    if (nrow(dat)==0) {
        return(list(bpue = NA, lrw = NA, upr = NA, model = "None",replicates=nrow(dat),
		 base_model_heterogeneity=NA,
		 bpue.cond.name=NA,
		 bpue.cond=NA,
		 bpue.cond.lwr=NA,
		 bpue.cond.upr=NA))
    }
    
    
    # cop out when we only have 1 row of data
    if (nrow(dat) == 1) {
        bpue <- dat$n_ind / dat$daysAtSea
        lwr <- bpue -1.96 * sqrt(dat$n_ind/dat$daysAtSea^2)
        upr <- bpue +1.96 * sqrt(dat$n_ind/dat$daysAtSea^2)
        return(list(bpue = bpue, lwr = lwr, upr = upr, model= "only one",replicates=nrow(dat),
		 base_model_heterogeneity=NA,
		 bpue.cond.name=NA,
		 bpue.cond=NA,
		 bpue.cond.lwr=NA,
		 bpue.cond.upr=NA))
    }
    
    base_model <- tryCatch(glmmTMB(n_ind ~ 1, offset = log(daysAtSea), family = nbinom2, data = dat),error=function(e) e$message)
    # If we have less than 5 rows of monitoring data, fit a simple model and we're done.

    if (nrow(dat) < 5) {
        
        best <- base_model
        # fit rate estimate model? go back to using a different assumption on family
        # use offset as denomitator, so that we can get variance of bpue estimates?
              heterogeneity_base<-tryCatch((rma.glmm(xi=n_ind,ti=(daysAtSea),measure="IRLN",data=dat)$QEp.Wld<0.05),error=function(e) e$message)# I^2 (or other mesaure) to select model?# I^2 (or other mesaure) to select model?
   
        
    } else {
    
        # If we have more than 5 rows of monitoring data, fit all possible
        # combinations of the base model and one or more of the terms in the
        # vector below, added to the model as random effects.
        re <- c("country", "areaCode","year", "metierLevel5", "vesselLength_group","samplingProtocol","monitoringMethod")
		#re <- c("Country", "Year", "MetierLevel5", "VesselLength_group")
        # but only consider terms where the number of unique values
        # is greater than min_re_obs
        re <- re[sapply(re, function(x) length(unique(dat[[x]]))) >= min_re_obs]
        re.i <- do.call(CJ, replicate(length(re), c(TRUE, FALSE), simplify = FALSE))

        # fit all candidate models 
        candidates <- apply(re.i, 1, function(i) {
            
            if (all(i == FALSE)) {
                return(base_model)
            }
            
            re <- sprintf("(1|%s)", re[unlist(i)])
            form <- sprintf("n_ind ~ 1 + %s", paste(re, collapse = " + "))
            form <- as.formula(form)
            tryCatch(glmmTMB(formula = form, offset=log(daysAtSea), family = nbinom2, data = dat),error=function(e) e$message,warning=function(w) w$message)
        })
        candidates.converged<-candidates[which(lapply(candidates,class)!="character")]
        scores <- sapply(candidates.converged, AIC) # or something other than AIC?
		if (length(scores)>1&sum(!is.na(scores))>0) {
        best <- candidates.converged[[which.min(scores)]]} else { # potential issue with ties? #should be ok, default to the first global min
                                              # which.min returns the FIRST min val
        best=base_model}
		
		
         heterogeneity_base<-tryCatch((rma.glmm(xi=n_ind,ti=(daysAtSea),measure="IRLN",data=dat)$QEp.Wld<0.05),error=function(e) e$message)# I^2 (or other mesaure) to select model?# I^2 (or other mesaure) to select model?
    
	}
    
    if(class(best)!="character") {
    bpue.r <- as.data.frame(emmeans(best,~1, type="response", offset=log(1))) # BPUE for one DaS we can push to total effort prediction instead next
	} else{
	bpue.r <-data.frame(response=NA,asymp.LCL=NA,asymp.UCL=NA)
	}
	
	if(formula(base_model)!=formula(best)) {
	bpue.re <- as.data.frame(ggpredict(best,terms=c(paste(names(best$frame)[2:(length(names(best$frame))-2)],sep=",")),condition=c(daysAtSea=1),type="re")) # BPUE for one DaS we can push to total effort prediction instead next
	#for glmmTMB need to do ggpredict instead
	
    return(list(bpue = bpue.r$response, 
         lwr = c(bpue.r$lower.CL,bpue.r$asymp.LCL), 
         upr = c(bpue.r$upper.CL,bpue.r$asymp.UCL),
         model = format(formula(best)),
		 replicates=nrow(dat),
		 base_model_heterogeneity=heterogeneity_base,
		 bpue.cond.name=NA,
		 bpue.cond=NA,
		 bpue.cond.lwr=NA,
		 bpue.cond.upr=NA))
		 
		 } else {
		 return(list(bpue = bpue.r$response, 
         lwr = c(bpue.r$lower.CL,bpue.r$asymp.LCL), 
         upr = c(bpue.r$upper.CL,bpue.r$asymp.UCL), #trick to grab either depending on model
         model = format(formula(best)),
		 replicates=nrow(dat),
		 base_model_heterogeneity=heterogeneity_base,
		 bpue.cond.name=NA,
		 bpue.cond=NA,
		 bpue.cond.lwr=NA,
		 bpue.cond.upr=NA))
		 
		 }
}





library(ggeffects)
library(data.table)
library(glmmTMB)
library(metafor)
library(emmeans)

THEDATA<-fread("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/monitor_effort_bycatch_ecoregion_area_species_v2THELIST.csv")


THEDATA<-subset(THEDATA,daysAtSea>0)

THEDATA<-subset(THEDATA,taxa_bycatch_monitor_ok==TRUE)



ecoregions<-unique(THEDATA$ecoregion)
ML4<-unique(THEDATA$metierL4)
species<-unique(THEDATA$species)

ALL.df<-aggregate(cbind(n_ind,daysAtSea)~ecoregion+metierL4+species,data=THEDATA,sum)
ALL.df[c("bpue","lwr","upr","model","replicates","base_model_heterogeneity","bpue.cond.name","bpue.cond","bpue.cond.lwr","bpue.cond.upr")]<-NA

#5:14

tic<-Sys.time()
for (i in 1:2000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
print(Sys.time())
flush.console()
}
print(Sys.time())
print(i)
flush.console()
}
Sys.time()-tic
al1<-ALL.df[1:2000,]
save(al1,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf1_2000.Rdata")


tic<-Sys.time()
for (i in 2001:4000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al2<-ALL.df[2001:4000,]
save(al2,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf2001_4000.Rdata")


tic<-Sys.time()
for (i in 4001:6000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al3<-ALL.df[4001:6000,]
save(al3,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf4001_6000.Rdata")


tic<-Sys.time()
for (i in 6001:8000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al4<-ALL.df[6001:8000,]
save(al4,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf6001_8000.Rdata")


tic<-Sys.time()
for (i in 8001:10000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al5<-ALL.df[8001:10000,]
save(al5,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf8001_10000.Rdata")


tic<-Sys.time()
for (i in 10001:12000) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al6<-ALL.df[10001:12000,]
save(al6,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf10001_12000.Rdata")


tic<-Sys.time()
for (i in 12001:13242) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,6:15]<-calc_bpue(ALL.df$ecoregion[i],ALL.df$metierL4[i],ALL.df$species[i],ices=THEDATA)
}
print(i)
flush.console()
}
Sys.time()-tic
al7<-ALL.df[12001:13242,]
save(al7,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean/alldf12001_13242.Rdata")


################################################################################################################################
################################################################################################################################
###################################################################################################################################
############### pre-BEAM BPUE
################################################################################################################################
################################################################################################################################
###################################################################################################################################



library(ggeffects)
library(data.table)
library(glmmTMB)
library(metafor)
library(emmeans)
setwd("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc2024/data/clean//")

THEDATA<-fread("monitor_effort_bycatch_ecoregion_area_species_v2THELIST.csv")

#####
##### some errors

THEDATA<-subset(THEDATA,DaysAtSea>0)
THEDATA<-subset(THEDATA,taxa_bycatch_monitor_ok==TRUE)

load("alldf1_2000.Rdata")
load("alldf2001_4000.Rdata")
load("alldf4001_6000.Rdata")
load("alldf6001_8000.Rdata")
load("alldf8001_10000.Rdata")
load("alldf10001_12000.Rdata")
load("alldf12001_13242.Rdata")


BPUE.est<-rbind(al1,al2,al3,al4,al5,al6,al7)

write.csv(BPUE.est,file="BPUE_est_2024.csv")



