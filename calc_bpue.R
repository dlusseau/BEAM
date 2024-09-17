
library(ggeffects)
library(data.table)
library(glmmTMB)
library(metafor)
library(emmeans)
setwd("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/")

THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")

#####
##### some errors

THEDATA<-subset(THEDATA,DaysAtSea>0)


#setwd("C:/Users/a21997/Documents/Møter/2023-09 ICES WGBYC")

#event <- fread("D3_Bycatch_event_Last5Years_v2.csv", dec = ".")
#monitor <- fread("D2_Bycatch_monitoring_effort_Last5Years_v2.csv")

# convert to numeric
#event[IndividualsWithoutPingers %in% c("-9", "NULL"), IndividualsWithoutPingers := NA]
#event[IndividualsWithPingers %in% c("-9", "NULL"), IndividualsWithPingers := NA]
#event[, IndividualsWithoutPingers := as.numeric(IndividualsWithoutPingers)]
#event[, IndividualsWithPingers := as.numeric(IndividualsWithPingers)]

# grab ecoregions from monitoring data by matching up
# corresponding area codes. This results in incorrect ecoregions
# for an unknown number of rows, but there is no easy way to deal
# with this for the time being.
#event[monitor, on = "AreaCode", Ecoregion := i.Ecoregion]

# reclassify vessel lengths into a binary factor variable with 
# two levels: small and large.
#event[, VS := factor(VesselLengthRange, 
                     # levels = c("NK", "VL0006", "VL0608", "VL0612", "VL0810", 
                                # "VL0010", "VL1012", "VL0015", "VL1215", "VL1218",
                                # "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
                     # labels = rep(c("small", "large"), times = c(7,9)))]

# monitor[, VS := factor(VesselLengthRange, 
                     # levels = c("NK", "VL0006", "VL0608", "VL0612", "VL0810", 
                                # "VL0010", "VL1012", "VL0015", "VL1215", "VL1218",
                                # "VL1518", "VL15XX", "VL1824", "VL2440", "VL40xx", "VL40XX"),
                     # labels = rep(c("small", "large"), times = c(7,9)))]

# # aggregate monitoring and bycatch data separately, over some columns,
# # currently ecoregion, ml4, ml5, country, and year
# data <- list(monitor = monitor[, .(das = sum(DaysAtSeaOb)), 
                              # .(ecoregion = Ecoregion, ml4 = MetierL4, vs = VS, 
                                # country = Country, ml5 = MetierL5, year = Year)],
            # bycatch = event[, .(totalbycatch = sum(IndividualsWithPingers) + sum(IndividualsWithoutPingers)), 
                            # .(ecoregion = Ecoregion, ml4 = MetierL4, species = Species, 
                              # vs = VS, country = Country, ml5 = MetierL5, year = Year)])

ecoregions<-unique(THEDATA$Ecoregion)
ML4<-unique(THEDATA$MetierL4)
species<-unique(THEDATA$Species)




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
         model = paste(format(formula(best)),collapse=""),
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
         model = paste(format(formula(best)),collapse=""),
		 replicates=nrow(dat),
		 base_model_heterogeneity=heterogeneity_base,
		 bpue.cond.name=NA,
		 bpue.cond=NA,
		 bpue.cond.lwr=NA,
		 bpue.cond.upr=NA))
		 
		 }
}





ecoregions<-unique(THEDATA$Ecoregion)
ML4<-unique(THEDATA$MetierL4)
species<-unique(THEDATA$Species)

ALL.df<-aggregate(n_ind~Ecoregion+MetierL4+Species,data=THEDATA,sum)
ALL.df[c("bpue","lwr","upr","model","replicates","base_model_heterogeneity","bpue.cond.name","bpue.cond","bpue.cond.lwr","bpue.cond.upr")]<-NA

#5:14

tic<-Sys.time()
for (i in 1:2468) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,5:14]<-calc_bpue(ALL.df$Ecoregion[i],ALL.df$MetierL4[i],ALL.df$Species[i])
}
print(i)
flush.console()
}
Sys.time()-tic
al1<-ALL.df[1:2468,]
save(al1,file="alldf1_2468.Rdata")


tic<-Sys.time()
for (i in 2469:4937) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,5:14]<-calc_bpue(ALL.df$Ecoregion[i],ALL.df$MetierL4[i],ALL.df$Species[i])
}
print(i)
flush.console()
}
Sys.time()-tic
al2<-ALL.df[2469:4937,]
save(al2,file="alldf2469_4937.Rdata")


tic<-Sys.time()
for (i in 4938:7405) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,5:14]<-calc_bpue(ALL.df$Ecoregion[i],ALL.df$MetierL4[i],ALL.df$Species[i])
}
print(i)
flush.console()
}
Sys.time()-tic
al3<-ALL.df[4938:7405,]
save(al3,file="alldf4938_7405.Rdata")


tic<-Sys.time()
for (i in 7406:9874) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,5:14]<-calc_bpue(ALL.df$Ecoregion[i],ALL.df$MetierL4[i],ALL.df$Species[i])
}
print(i)
flush.console()
}
Sys.time()-tic
al4<-ALL.df[7406:9874,]
save(al4,file="alldf7406_9874.Rdata")


tic<-Sys.time()
for (i in 9875:12338) {
if(ALL.df$n_ind[i]>0) {
ALL.df[i,5:14]<-calc_bpue(ALL.df$Ecoregion[i],ALL.df$MetierL4[i],ALL.df$Species[i])
}
print(i)
flush.console()
}
Sys.time()-tic
al5<-ALL.df[9875:12338,]
save(al5,file="alldf9875_12338.Rdata")





library(ggeffects)
library(data.table)
library(glmmTMB)
library(metafor)
library(emmeans)
setwd("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/")

THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")

#####
##### some errors

THEDATA<-subset(THEDATA,DaysAtSea>0)

load("alldf1_2468.Rdata")
load("alldf2469_4937.Rdata")
load("alldf4938_7405.Rdata")
load("alldf7406_9874.Rdata")
load("alldf9875_12338.Rdata")


BPUE.est<-rbind(al1,al2,al3,al4,al5)

write.csv(BPUE.est,file("BPUE_est_2023.csv")


#############################
#let's start to weed out those we can't say anything about
BPUE.est$col<-"red"

BPUE.retain<-BPUE.est[!is.na(BPUE.est$model),]
BPUE.retain<-BPUE.retain[-which(BPUE.retain$model=="n_ind ~ 1" & BPUE.retain$base_model_heterogeneity!=FALSE),]
#BPUE.retain<-BPUE.retain[-which(is.na(BPUE.retain$bpue.cond)&((BPUE.retain$model!="n_ind ~ 1")|(BPUE.retain$model!="only one"))),]
# we have values for all

BPUE.retain$colour<-"green"
BPUE.retain$colour[BPUE.retain$replicates<5]<-"yellow"
BPUE.retain$QC1<-"green"
BPUE.retain$QC2<-"green"
BPUE.retain$QC3<-"green"



### the writing to column of the list element did not work :( we have to refit the 132 models

# 552 estimates
# 49 yellows (less than 5 replicates)
#132 random effect models to check against effort
#
# check are there monitor 0 and effort non.zero values?
# aggregate effort at Ecoregion, ML4
#refit each model and predict for the total effort offset for the n_ind ~1 models
# refit the 132 models with random effects 




# calc bpue for all ecoregion x ml4 x species combinations

calc_bpue("Greater North Sea", "GNS", "Phocoena phocoena")
calc_bpue("Greater North Sea", "GNS", "Alca torda")
calc_bpue("Greater North Sea", "GNS", "Chelidonichthys lucerna")

# bulk processing
bpue <- data$bycatch[ecoregion=="Greater North Sea" & ml4 == "GNS", 
                     calc_bpue(ecoregion, ml4, species), .(ecoregion, ml4, species)]
					 
					 

					 
# bpue<-tryCatch(glmmTMB(TotalBycatch~1+(1|row.names(data.temp.re)),offset=log(DaysAtSeaOb),family=nbinom2,data=data.temp.re),
# error=function(e) e$message,
                # warning=function(w) w$message)
