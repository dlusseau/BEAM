

#######################################################################################
######################################################################################
setwd("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/")

effort<-read.csv("Ag_effort_by_EcoR_Country_ML4_VS.csv")
BPUE.final<-read.csv("BPUE_INTERMEDIATE.csv")

library(data.table)
library(ggeffects)
library(glmmTMB)


THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")
THEDATA<-subset(THEDATA,DaysAtSea>0)



#variables retained in models
#Country Year (for year we just want the 2022 random intercept),VesselLength_group
#monitoring methods and sampling method - we need to go back to the monitored effort - estimate the ratio fo method and apply this ratio to total effort
BPUE.final$mean_TB<-NA
BPUE.final$lwr_TB<-NA
BPUE.final$upr_TB<-NA
BPUE.final<-subset(BPUE.final,model!="only one")

for (i in 1:dim(BPUE.final)[1]) {
dat<-THEDATA[THEDATA$Ecoregion==BPUE.final$Ecoregion[i]& THEDATA$MetierL4==BPUE.final$MetierL4[i] & THEDATA$Species==BPUE.final$Species[i],]
form<-as.formula(BPUE.final$model[i])
best<-glmmTMB(formula = form, offset=log(DaysAtSea), family = nbinom2, data = dat)
effort.sub<-effort[effort$Ecoregion==BPUE.final$Ecoregion[i] & effort$MetierL4==BPUE.final$MetierL4[i] ,]
if (form=="n_ind ~ 1") {
TE<-sum(effort.sub$FishingEffort)
bpue.re <-  as.data.frame(emmeans(best,~1,offset=log(TE),type="response"))
BPUE.final$mean_TB[i]<-(bpue.re$response) 
BPUE.final$lwr_TB[i]<-(bpue.re$asymp.LCL)
BPUE.final$upr_TB[i]<-(bpue.re$asymp.UCL)

} else {

if (!(("SamplingProtocol" %in% term )|("MonitoringMethod" %in% term ))) {

if (form=="n_ind ~ 1 + (1 | Country)") {
countries<-unique(dat$Country)
countries.ef<-unique(effort.sub$Country)

if ((sum(countries.ef%in%countries)/length(countries.ef)==1)&(sum(countries%in%countries.ef)/length(countries)==1)) {
TE<-aggregate(FishingEffort~Country,data=effort.sub,sum)

bpue.re<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[1],"]",sep="")),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[1]]),type="re"))

for (j in 2:length(countries)) {
bpue.rex<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[j],"]",sep="")),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[j]]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

} #end for countries
BPUE.final$mean_TB[i]<-(bpue.re$predicted)
BPUE.final$lwr_TB[i]<-(bpue.re$conf.low)
BPUE.final$upr_TB[i]<-(bpue.re$conf.high)
} #end country match check
} #end country

if (form=="n_ind ~ 1 + (1 | Year)") {

TE<-sum(effort.sub$FishingEffort)
bpue.re<-as.data.frame(ggpredict(best,terms=c("Year [2022]"),condition=c(DaysAtSea=TE),type="re"))
BPUE.final$mean_TB[i]<-(bpue.re$predicted) 
BPUE.final$lwr_TB[i]<-(bpue.re$conf.low)
BPUE.final$upr_TB[i]<-(bpue.re$conf.high)

}

if (form=="n_ind ~ 1 + (1 | VesselLength_group)") {

TE<-aggregate(FishingEffort~VS,data=effort.sub,sum)

bpue.re<-as.data.frame(ggpredict(best,terms=c("VesselLength_group [below_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$VS=="small"]),type="re"))

bpue.rex<-as.data.frame(ggpredict(best,terms=c("VesselLength_group [above_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$VS=="large"]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

BPUE.final$mean_TB[i]<-(bpue.re$predicted)
BPUE.final$lwr_TB[i]<-(bpue.re$conf.low)
BPUE.final$upr_TB[i]<-(bpue.re$conf.high)

}

if (form=="n_ind ~ 1 + (1 | Country) + (1 | Year)") {

countries<-unique(dat$Country)
countries.ef<-unique(effort.sub$Country)

if ((sum(countries.ef%in%countries)/length(countries.ef)==1)&(sum(countries%in%countries.ef)/length(countries)==1)) {
TE<-aggregate(FishingEffort~Country,data=effort.sub,sum)

bpue.re<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[1],"]",sep=""),"Year [2022]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[1]]),type="re"))

for (j in 2:length(countries)) {
bpue.rex<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[j],"]",sep=""),"Year [2022]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[j]]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

} #end for countries
BPUE.final$mean_TB[i]<-(bpue.re$predicted)
BPUE.final$lwr_TB[i]<-(bpue.re$conf.low)
BPUE.final$upr_TB[i]<-(bpue.re$conf.high)
} #end country match check


}

if (form=="n_ind ~ 1 + (1 | Country) + (1 | Year) + (1 | VesselLength_group)") {

countries<-unique(dat$Country)
countries.ef<-unique(effort.sub$Country)

if ((sum(countries.ef%in%countries)/length(countries.ef)==1)&(sum(countries%in%countries.ef)/length(countries)==1)) {
TE<-aggregate(FishingEffort~Country+VS,data=effort.sub,sum)

bpue.re<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[1],"]",sep=""),"Year [2022]","VesselLength_group[below_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[1]&TE$VS=="small"]),type="re"))

for (j in 2:length(countries)) {
bpue.rex<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[j],"]",sep=""),"Year [2022]","VesselLength_group[below_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[j]&TE$VS=="small"]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

} #end for countries

bpue.rex<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[1],"]",sep=""),"Year [2022]","VesselLength_group[above_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[1]&TE$VS=="large"]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

for (j in 2:length(countries)) {
bpue.rex<-as.data.frame(ggpredict(best,terms=c(paste("Country [",countries[j],"]",sep=""),"Year [2022]","VesselLength_group[above_12]"),condition=c(DaysAtSea=TE$FishingEffort[TE$Country==countries[j]&TE$VS=="large"]),type="re"))
bpue.re$predicted<-bpue.re$predicted+bpue.rex$predicted
bpue.re$conf.low<-bpue.re$conf.low+bpue.rex$conf.low
bpue.re$conf.high<-bpue.re$conf.high+bpue.rex$conf.high

} #end for countries


BPUE.final$mean_TB[i]<-(bpue.re$predicted)
BPUE.final$lwr_TB[i]<-(bpue.re$conf.low)
BPUE.final$upr_TB[i]<-(bpue.re$conf.high)
} #end country match check


}


} # end if term SampProt and MonMet exception

} #end else fe v re
print(i)
flush.console()
} #end i

write.csv(BPUE.final,file="BPUE_final_w_TB.csv")

BPUE.final$CV.rel<-(BPUE.final$upr_TB-BPUE.final$lwr_TB)/BPUE.final$mean_TB # a measure of relative uncertainty (CI/mean)
BPUE.final$monitor.effort<-NA
BPUE.final$total.effort<-NA

for (i in 1:dim(BPUE.final)[1]) {
dat<-THEDATA[THEDATA$Ecoregion==BPUE.final$Ecoregion[i]& THEDATA$MetierL4==BPUE.final$MetierL4[i] & THEDATA$Species==BPUE.final$Species[i],]
effort.sub<-effort[effort$Ecoregion==BPUE.final$Ecoregion[i] & effort$MetierL4==BPUE.final$MetierL4[i] ,]
BPUE.final$monitor.effort[i]<-sum(dat$DaysAtSea)
BPUE.final$total.effort[i]<-sum(effort.sub$FishingEffort)

}


write.csv(BPUE.final,file="BPUE_final_w_TB.csv")

#############################################################################
#############################################################################
####### visualise analytical summary

tried<-aggregate(DaysAtSea~Ecoregion+MetierL4+Species,data=subset(THEDATA,n_ind>0),"sum")

achieved<-subset(BPUE.final,is.na(upr_TB)==F)

tried$achieved<-0
for (i in 1:nrow(achieved)) {
tried$achieved[which(tried$Ecoregion==achieved$Ecoregion[i] & tried$MetierL4==achieved$MetierL4[i] & tried$Species==achieved$Species[i])]<-1

}

glm0<-glm(achieved~DaysAtSea,data=tried,family=binomial)
library(ggeffects)

pred0<-ggpredict(glm0)
plot(pred0)

tried$Ecoregion<-factor(tried$Ecoregion)
tried$MetierL4<-factor(tried$MetierL4)
tried$Species<-factor(tried$Species)

tab.E.M<-table(tried[tried$achieved==1,]$Ecoregion,tried[tried$achieved==1,]$MetierL4)/table(tried$Ecoregion,tried$MetierL4)
tab.E.M[which(table(tried$Ecoregion,tried$MetierL4)==0)]<-NA

library(reshape2)
po=melt(tab.E.M)
head(po)
names(po)<-c("Ecoregion","MetierL4","value")

library(ggplot2)
library(scales) # for muted function
gg.E.M<-ggplot(po, aes(MetierL4,Ecoregion)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes(fill = po$value, label = round(po$value, 2))) + # write the values
  scale_fill_gradient2(name="Species proportion",
					   low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0.5) + 
	xlab("Metier level 4")+
	theme_minimal()
	
tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/summary_table_prop_Eco_metier.tiff",width=30,height=10,units="cm",res=200)
gg.E.M	
dev.off()
##################################################################


tab.S.M<-table(tried[tried$achieved==1,]$Species,tried[tried$achieved==1,]$MetierL4)/table(tried$Species,tried$MetierL4)
tab.S.M[which(table(tried$Species,tried$MetierL4)==0)]<-NA

library(reshape2)
po=melt(tab.S.M)
head(po)
names(po)<-c("Species","MetierL4","value")

library(ggplot2)
library(scales) # for muted function
gg.S.M<-ggplot(po, aes(MetierL4,Species)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes(fill = po$value, label = round(po$value, 2))) + # write the values
  scale_fill_gradient2(name="Ecoregion proportion",
					   low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0.5) + 
	xlab("Metier level 4")+
	theme_minimal()
	
tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/summary_table_prop_Species_metier.tiff",width=25,height=32,units="cm",res=200)
gg.S.M	
dev.off()


########################
###### BPUE beam table
BPUE.final

library(taxize)
library(rlang)
species.class<-tax_name(BPUE.final$Species,get="class",db="itis")
species.df<-data.frame(species=unique(BPUE.final$Species),common=NA,class=NA)
species.df$common<-sci2comm(species.df$species,db="itis")
missing<- names(unlist(lapply(species.df$common,function(x) which(is_empty(x)))))
species.df$species%in%missing

species.df$common[species.df$species%in%missing]<-sci2comm(species.df$species[species.df$species%in%missing],db="ncbi")
species.df$common[species.df$species=="Dipturus intermedius"]<-"flapper skate"
species.df$common.main<-unlist(lapply(species.df$common,function(x) unlist(split(x,","))[1] ))
## quite a mess
species.df$common.main[1]<-unlist(species.df$common[1])[4]
species.df$common.main[5]<-unlist(species.df$common[5])[3]
species.df$common.main[21]<-unlist(species.df$common[21])[3]
species.df$common.main[30]<-unlist(species.df$common[30])[2]
species.df$common.main[40]<-unlist(species.df$common[40])[2]
species.df$common.main[41]<-unlist(species.df$common[41])[2]
species.df$common.main[66]<-unlist(species.df$common[66])[2]
species.df$common.main[69]<-unlist(species.df$common[69])[2]
species.df$common.main[77]<-unlist(species.df$common[77])[2]
species.df$common.main[80]<-unlist(species.df$common[80])[5]

BPUE.final$taxon<-"Fish"
BPUE.final$taxon[which(species.class$class=="Aves")]<-"Bird"
BPUE.final$taxon[which(species.class$class=="Mammalia")]<-"Mammals"
BPUE.final$taxon[which(species.class$class=="Reptilia")]<-"Reptiles"

BPUE.final$common<-species.df$common.main[match(BPUE.final$Species,species.df$species)]


achieved<-subset(BPUE.final,is.na(upr_TB)==F)
achieved$label<-paste(achieved$Ecoregion,achieved$MetierL4,achieved$common,sep=" and ")


mammals<-ggplot(subset(achieved,taxon=="Mammals"& label!="Norwegian Sea and GNS and Harbor Porpoise"& label!="Azores and LHM and Short-beaked Common Dolphin"), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()+
  ylim(0,5)

birds<-ggplot(subset(achieved,taxon=="Bird"), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()


reptiles<-ggplot(subset(achieved,taxon=="Reptiles"), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()+
  ylim(0,4)
 
 fish<-ggplot(subset(achieved,taxon=="Fish"&log10(upr_TB)<12&log10(lwr_TB)>(-1)), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()
 

	
tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_fish.tiff",width=20,height=30,units="cm",res=200)
fish
dev.off()


tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_mammals.tiff",width=20,height=10,units="cm",res=200)
mammals
dev.off()


tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_birds.tiff",width=20,height=10,units="cm",res=200)
birds
dev.off()


tiff("C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_reptiles.tiff",width=20,height=10,units="cm",res=200)
reptiles
dev.off()

######

table.print<-achieved[,c("taxon","Ecoregion","MetierL4","common","monitor.effort","total.effort","model","bpue","lwr","upr","lwr_TB","upr_TB")]

table.print$monitor.effort<-round(table.print$monitor.effort,0)
table.print$total.effort<-round(table.print$total.effort,0)
table.print$bpue<-round(table.print$bpue,6)
table.print$lwr<-round(table.print$lwr,7)
table.print$upr<-round(table.print$upr,7)
table.print$lwr_TB<-round(log10(table.print$lwr_TB),2)
table.print$upr_TB<-round(log10(table.print$upr_TB),2)


table.print$BPUE<-paste0(formatC(table.print$bpue,digits=5,format="f")," [",table.print$lwr," ; ",table.print$upr,"]")
#table.print$upr_TB<-formatC(as.numeric(table.print$upr_TB),format="d")

table.print$representability<-"a constant BPUE appears to be representative"
table.print$representability[table.print$model=="n_ind ~ 1 + (1 | Year)"]<-"there is between-year variability in BPUE"
table.print$representability[table.print$model=="n_ind ~ 1 + (1 | VesselLength_group)"]<-"there is between-vessel length category variability in BPUE"

table.final<-table.print[,c(1:6,13:14,11:12)]

write.csv(table.final,file="C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/retainedBPUE.csv")



#################################################################################################################
###### let's table up the QC process

library(taxize)
library(rlang)

THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")
THEDATA<-subset(THEDATA,DaysAtSea>0)



table.all<-aggregate(n_ind ~ Ecoregion + MetierL4 + Species,data=THEDATA,sum)
table.all<-aggregate(cbind(n_ind,DaysAtSea) ~ Ecoregion + MetierL4 + Species,data=THEDATA,sum)

no.bycatch<-subset(table.all,n_ind==0) # here we have 584 non-zero count for which we can engage with bpue estimation (at least for now the zeros are not dealt with, however in the future we will concentrate on finding a ways to extract information for some of theme_minimal
nob.species.df<-data.frame(species=unique(no.bycatch$Species),common=NA,class=NA)
nob.species.df$class<-tax_name(nob.species.df$species,get="class",db="itis")
#nob.species.df$common<-sci2comm(nob.species.df$species,db="ncbi")
no.bycatch$taxon<-"Fish"

no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Mammalia"&is.na(nob.species.df$class$class)==FALSE]]<-"Mammals"
no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Aves"&is.na(nob.species.df$class$class)==FALSE]]<-"Birds"
no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Reptilia"&is.na(nob.species.df$class$class)==FALSE]]<-"Reptilia"

no.bycatch$DaysAtSea<-ceiling(no.bycatch$DaysAtSea)

write.csv(no.bycatch,file="C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/no_bycatch_detected.csv")

matchingmonitor<-read.csv("bpue_retain_taxa_issue.csv")

#QC1 can we fit a model and have enough monitoring
bpue.first<-read.csv("BPUE_est_2023.csv")
bpue.first<-subset(bpue.first,is.na(bpue)==FALSE)

df.QC1<-bpue.first[bpue.first$model=="n_ind ~ 1" & bpue.first$ base_model_heterogeneity==T,]
write.csv(df.QC1,"QC1_interceptonlymodel_heterogeneous.csv")

df.QC1replicates<-bpue.first[bpue.first$model=="n_ind ~ 1" & bpue.first$base_model_heterogeneity==F&bpue.first$replicates<5,]
write.csv(df.QC1replicates,"QC1_interceptonlymodel_homogeneous_lowreplicates.csv")

df.QC1replicates.table<-df.QC1replicates[,c(2,3,4,6,7,8)]
df.QC1replicates.table$bpue<-round(df.QC1replicates.table$bpue,3)
df.QC1replicates.table$lwr<-round(df.QC1replicates.table$lwr,4)
df.QC1replicates.table$upr<-round(df.QC1replicates.table$upr,4)
write.table(df.QC1replicates.table,"QC1_interceptonlymodel_homogeneous_lowreplicates_table.txt",dec = ".",sep=";")


df.QC2.QC3<-read.csv("BPUE_retain_QC2_QC3.csv")

########
#finally the taxonomic mismatch between species and observatoin method

#########OCT 2023 check why grey seal celtic and Phocoena North sea were removed?
df.taxa<-read.csv("bpue_retain_taxa_issue.csv")
df.retained<-df.taxa[df.taxa$taxa_bycatch_monitor_ok==FALSE&is.na(df.taxa$taxa_bycatch_monitor_ok)==FALSE&df.taxa$QC2=="green"&df.taxa$QC3=="green",]
write.csv(df.retained[,c(2:4)]y,file="bycatch_rejected_mismatch_species_method.csv")

##################################################################################
# the Oct 2023 problem is indeed from the NAs
#####################################################################################

monitoring.original<-read.csv("D2_Bycatch_monitoring_effort_Last5Years_v2.csv")

missing<-df.taxa[is.na(df.taxa$taxa_bycatch_monitor_ok)==TRUE,]

missing.levels<-unique(missing[c("Ecoregion","MetierL4")])
missing.levels$Y1<-NA
missing.levels$Y2<-NA
missing.levels$Y3<-NA
missing.levels$Y4<-NA
missing.levels$Y5<-NA

for (i in 1:dim(missing.levels)[1]) {
tryCatch(missing.levels$Y1[i]<-names(which.max(table(unlist(c(monitoring.original$SamplingProtocol[monitoring.original$Ecoregion==missing.levels$Ecoregion[i] & monitoring.original$MetierL4==missing.levels$MetierL4[i] & monitoring.original$SubmissionYear==2019]))))),error=function(e) NA)
tryCatch(missing.levels$Y2[i]<-names(which.max(table(unlist(c(monitoring.original$SamplingProtocol[monitoring.original$Ecoregion==missing.levels$Ecoregion[i] & monitoring.original$MetierL4==missing.levels$MetierL4[i] & monitoring.original$SubmissionYear==2020]))))),error=function(e) NA)
tryCatch(missing.levels$Y3[i]<-names(which.max(table(unlist(c(monitoring.original$SamplingProtocol[monitoring.original$Ecoregion==missing.levels$Ecoregion[i] & monitoring.original$MetierL4==missing.levels$MetierL4[i] & monitoring.original$SubmissionYear==2021]))))),error=function(e) NA)
tryCatch(missing.levels$Y4[i]<-names(which.max(table(unlist(c(monitoring.original$SamplingProtocol[monitoring.original$Ecoregion==missing.levels$Ecoregion[i] & monitoring.original$MetierL4==missing.levels$MetierL4[i] & monitoring.original$SubmissionYear==2022]))))),error=function(e) NA)
tryCatch(missing.levels$Y5[i]<-names(which.max(table(unlist(c(monitoring.original$SamplingProtocol[monitoring.original$Ecoregion==missing.levels$Ecoregion[i] & monitoring.original$MetierL4==missing.levels$MetierL4[i] & monitoring.original$SubmissionYear==2023]))))),error=function(e) NA)

}

#so essentially for the NAs, we have different SO for different Metiers at levels higher than 4, in most cases it is ok because the primary SO is All, so detected bycatch were looked for.
missing.levels$LC<-NA
for (i in 1:dim(missing.levels)[1]) {
missing.levels$LC<-names(which.max(table(unlist(c(missing.levels[i,3:7])))))

}

###############################################################################################################################################################################################################################
###############################################################################################################################################################################################################################
###### replicate the whole process for the Med sublist
#######################################################################################
######################################################################################
setwd("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/")

# effort<-read.csv("Ag_effort_by_EcoR_Country_ML4_VS.csv")
# BPUE.final<-read.csv("BPUE_INTERMEDIATE.csv")

library(data.table)
library(ggeffects)
library(glmmTMB)
library(metafor)
library(emmeans)


THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")
THEDATA<-subset(THEDATA,DaysAtSea>0)

####

THEDATAMED<-fread("Med_monitor_effort_bycatch_ecoregion_species.csv", dec=".")
THEDATAMED<-subset(THEDATAMED,DaysAtSea>0)

unique(THEDATA$Ecoregion)
medregion<-unique(THEDATAMED$Ecoregion)[!(unique(THEDATAMED$Ecoregion)%in%unique(THEDATA$Ecoregion))]

MED<-THEDATAMED[THEDATAMED$Ecoregion%in%medregion,]

#load the beautiful calc_bpue function




ecoregions<-unique(MED$Ecoregion)
ML4<-unique(MED$MetierL4)
species<-unique(MED$Species)

MED.df<-aggregate(n_ind~Ecoregion+MetierL4+Species,data=MED,sum)
MED.df[c("bpue","lwr","upr","model","replicates","base_model_heterogeneity","bpue.cond.name","bpue.cond","bpue.cond.lwr","bpue.cond.upr")]<-NA

#5:14

tic<-Sys.time()
for (i in 158:nrow(MED.df)) {
if(MED.df$n_ind[i]>0) {
MED.df[i,5:14]<-calc_bpue(MED.df$Ecoregion[i],MED.df$MetierL4[i],MED.df$Species[i],min_re_obs = 2,THEDATA=MED)
}
print(i)
flush.console()
}
Sys.time()-tic

save(MED.df,file="meddf.Rdata")

BPUE.est<-MED.df



#############################
#let's start to weed out those we can't say anything about
BPUE.est$colour<-"red"

BPUE.retain<-BPUE.est[!is.na(BPUE.est$model),]
BPUE.retain<-BPUE.retain[-which(BPUE.retain$model=="n_ind ~ 1" & BPUE.retain$base_model_heterogeneity!=FALSE),]
#BPUE.retain<-BPUE.retain[-which(is.na(BPUE.retain$bpue.cond)&((BPUE.retain$model!="n_ind ~ 1")|(BPUE.retain$model!="only one"))),]
# we have values for all

BPUE.retain$colour<-"green"
BPUE.retain$colour[BPUE.retain$replicates<5]<-"yellow"
BPUE.retain$QC1<-"green"
BPUE.retain$QC2<-"green"
BPUE.retain$QC3<-"green"

effort<-read.csv("Ag_effort_by_EcoR_Country_ML4_VS.csv")

BPUE.retain$QC3[20]<-"red"
#remove mismatch observatoin species

BPUE.final<- BPUE.retain[!(BPUE.retain$Ecoregion=="Aegean-Levantine Sea"&BPUE.retain$MetierL4=="OTB"&BPUE.retain$Species=="Caretta caretta"),]

#variables retained in models
#Country Year (for year we just want the 2022 random intercept),VesselLength_group
#monitoring methods and sampling method - we need to go back to the monitored effort - estimate the ratio fo method and apply this ratio to total effort
BPUE.final$mean_TB<-NA
BPUE.final$lwr_TB<-NA
BPUE.final$upr_TB<-NA
BPUE.final<-subset(BPUE.final,model!="only one")

for (i in 1:dim(BPUE.final)[1]) {
dat<-MED[MED$Ecoregion==BPUE.final$Ecoregion[i]& MED$MetierL4==BPUE.final$MetierL4[i] & MED$Species==BPUE.final$Species[i],]

if (BPUE.final$model[i]=="n_ind ~ 1") {
form<-as.formula(BPUE.final$model[i])
best<-tryCatch(glmmTMB(formula = form, offset=log(DaysAtSea), family = nbinom2, data = dat),warning=function(w) w$message)
effort.sub<-effort[effort$Ecoregion==BPUE.final$Ecoregion[i] & effort$MetierL4==BPUE.final$MetierL4[i] ,]
if (is.character(best)==FALSE) {

TE<-sum(effort.sub$FishingEffort)
bpue.re <-  as.data.frame(emmeans(best,~1,offset=log(TE),type="response"))
BPUE.final$mean_TB[i]<-(bpue.re$response) 
BPUE.final$lwr_TB[i]<- bpue.re$lower.CL#(bpue.re$asymp.LCL)
BPUE.final$upr_TB[i]<-  bpue.re$upper.CL#(bpue.re$asymp.UCL)

} 
}

print(i)
flush.console()
} #end i

write.csv(BPUE.final,file="BPUE_final_w_TB_MED.csv")

BPUE.final$CV.rel<-(BPUE.final$upr_TB-BPUE.final$lwr_TB)/BPUE.final$mean_TB # a measure of relative uncertainty (CI/mean)
BPUE.final$monitor.effort<-NA
BPUE.final$total.effort<-NA

for (i in 1:dim(BPUE.final)[1]) {
dat<-MED[MED$Ecoregion==BPUE.final$Ecoregion[i]& MED$MetierL4==BPUE.final$MetierL4[i] & MED$Species==BPUE.final$Species[i],]
effort.sub<-effort[effort$Ecoregion==BPUE.final$Ecoregion[i] & effort$MetierL4==BPUE.final$MetierL4[i] ,]
BPUE.final$monitor.effort[i]<-sum(dat$DaysAtSea)
BPUE.final$total.effort[i]<-sum(effort.sub$FishingEffort)

}


write.csv(BPUE.final,file="BPUE_final_w_TB_MED.csv")

#######################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################
####### visualise analytical summary

#integrate both BPUE.final
BPUE.final.MED<-BPUE.final

BPUE.final.nonMED<-read.csv("BPUE_final_w_TB.csv")
BPUE.final.MED<-BPUE.final.MED[,names(BPUE.final.nonMED)[2:21]]
BPUE.final.nonMED<-BPUE.final.nonMED[,-1]

BPUE.final<-rbind(BPUE.final.nonMED,BPUE.final.MED)

tried<-aggregate(DaysAtSea~Ecoregion+MetierL4+Species,data=subset(THEDATAMED,n_ind>0),"sum")

achieved<-subset(BPUE.final,is.na(upr_TB)==F)

tried$achieved<-0
for (i in 1:nrow(achieved)) {
tried$achieved[which(tried$Ecoregion==achieved$Ecoregion[i] & tried$MetierL4==achieved$MetierL4[i] & tried$Species==achieved$Species[i])]<-1

}

tried$Ecoregion<-factor(tried$Ecoregion)
tried$MetierL4<-factor(tried$MetierL4)
tried$Species<-factor(tried$Species)

tab.E.M<-table(tried[tried$achieved==1,]$Ecoregion,tried[tried$achieved==1,]$MetierL4)/table(tried$Ecoregion,tried$MetierL4)
tab.E.M[which(table(tried$Ecoregion,tried$MetierL4)==0)]<-NA

library(reshape2)
po=melt(tab.E.M)
head(po)
names(po)<-c("Ecoregion","MetierL4","value")

library(ggplot2)
library(scales) # for muted function
gg.E.M<-ggplot(po, aes(MetierL4,Ecoregion)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes(fill = po$value, label = round(po$value, 2))) + # write the values
  scale_fill_gradient2(name="Species proportion",
					   low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0.5) + 
	xlab("Metier level 4")+
	theme_minimal()
	
tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/summary_table_prop_Eco_metier_inclMED.tiff",width=30,height=10,units="cm",res=200)
gg.E.M	
dev.off()
##################################################################


tab.S.M<-table(tried[tried$achieved==1,]$Species,tried[tried$achieved==1,]$MetierL4)/table(tried$Species,tried$MetierL4)
tab.S.M[which(table(tried$Species,tried$MetierL4)==0)]<-NA

library(reshape2)
po=melt(tab.S.M)
head(po)
names(po)<-c("Species","MetierL4","value")

library(ggplot2)
library(scales) # for muted function
gg.S.M<-ggplot(po, aes(MetierL4,Species)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes(fill = po$value, label = round(po$value, 2))) + # write the values
  scale_fill_gradient2(name="Ecoregion proportion",
					   low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0.5) + 
	xlab("Metier level 4")+
	theme_minimal()
	
tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/summary_table_prop_Species_metier_inclMED.tiff",width=25,height=32,units="cm",res=200)
gg.S.M	
dev.off()


########################
###### BPUE beam table
final.taxon<-read.csv("BPUE.final_w_taxon.csv")

species.df<-data.frame(species=unique(final.taxon$Species),taxon=NA,common=NA)

species.df$taxon<-final.taxon$taxon[match(species.df$species,final.taxon$Species)]
species.df$common<-final.taxon$common[match(species.df$species,final.taxon$Species)]

unique(BPUE.final$Species)[!(unique(BPUE.final$Species)%in%species.df$species)]

species.df<-rbind(species.df,data.frame(species="Monachus monachus",taxon="Mammals",common="Mediterranean monk seal"))
species.df<-rbind(species.df,data.frame(species="Chelonia mydas",taxon="Reptiles",common="Green sea turtle"))
species.df<-rbind(species.df,data.frame(species="Gymnura altavela",taxon="Fish",common="Spiny butterfly ray"))
species.df<-rbind(species.df,data.frame(species="Huso huso",taxon="Fish",common="beluga sturgeon"))

BPUE.final$taxon<-species.df$taxon[match(BPUE.final$Species,species.df$species)]
BPUE.final$common<-species.df$common[match(BPUE.final$Species,species.df$species)]

write.csv(BPUE.final,file="BPUE_final_inclMED_TB_taxon.csv")

achieved<-subset(BPUE.final,is.na(upr_TB)==F)
achieved$label<-paste(achieved$Ecoregion,achieved$MetierL4,achieved$common,sep=" and ")


mammals<-ggplot(subset(achieved,taxon=="Mammals"& label!="Norwegian Sea and GNS and Harbor Porpoise"& label!="Azores and LHM and Short-beaked Common Dolphin" & label!="Aegean-Levantine Sea and LLS and Mediterranean monk seal"), #removing the few with unsensible CI
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()+
  ylim(0,5)

birds<-ggplot(subset(achieved,taxon=="Bird"), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()


reptiles<-ggplot(subset(achieved,taxon=="Reptiles"), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()+
  ylim(0,5)
 
 fish<-ggplot(subset(achieved,taxon=="Fish"&log10(upr_TB)<12&log10(lwr_TB)>(-1)), 
  aes(x = label, y=log10(lwr_TB))) +
  geom_crossbar(aes(ymin = log10(lwr_TB), ymax = log10(upr_TB)), width = 0.5, fill = "blue",fatten=0) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (log10 individuals)")+
  theme_minimal()
 

	
tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_fish_inclMed.tiff",width=20,height=35,units="cm",res=200)
fish
dev.off()


tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_mammals_inclMed.tiff",width=20,height=10,units="cm",res=200)
mammals
dev.off()


tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_birds_inclMed.tiff",width=20,height=10,units="cm",res=200)
birds
dev.off()


tiff("C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/totalbycatch_reptiles_inclMed.tiff",width=20,height=10,units="cm",res=200)
reptiles
dev.off()

######

table.print<-achieved[,c("taxon","Ecoregion","MetierL4","common","monitor.effort","total.effort","model","bpue","lwr","upr","lwr_TB","upr_TB")]

table.print$monitor.effort<-round(table.print$monitor.effort,0)
table.print$total.effort<-round(table.print$total.effort,0)
table.print$bpue<-round(table.print$bpue,6)
table.print$lwr<-round(table.print$lwr,7)
table.print$upr<-round(table.print$upr,7)
table.print$lwr_TB<-round(log10(table.print$lwr_TB),2)
table.print$upr_TB<-round(log10(table.print$upr_TB),2)


table.print$BPUE<-paste0(formatC(table.print$bpue,digits=5,format="f")," [",table.print$lwr," ; ",table.print$upr,"]")
#table.print$upr_TB<-formatC(as.numeric(table.print$upr_TB),format="d")

table.print$representability<-"a constant BPUE appears to be representative"
table.print$representability[table.print$model=="n_ind ~ 1 + (1 | Year)"]<-"there is between-year variability in BPUE"
table.print$representability[table.print$model=="n_ind ~ 1 + (1 | VesselLength_group)"]<-"there is between-vessel length category variability in BPUE"

table.final<-table.print[,c(1:6,13:14,11:12)]

write.csv(table.final,file="C:/Users/David/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/retainedBPUEinclMed.csv")

write.csv(table.final$lwr_TB,file="stupidXL_lower.csv",row.names =FALSE)
write.csv(table.final$upr_TB,file="stupidXL_upper.csv",row.names =FALSE)

#################################################################################################################
###### let's table up the QC process

library(taxize)
library(rlang)

THEDATA<-fread("monitor_effort_bycatch_ecoregion_species.csv", dec = ".")
THEDATA<-subset(THEDATA,DaysAtSea>0)



table.all<-aggregate(n_ind ~ Ecoregion + MetierL4 + Species,data=THEDATA,sum)
table.all<-aggregate(cbind(n_ind,DaysAtSea) ~ Ecoregion + MetierL4 + Species,data=THEDATA,sum)

no.bycatch<-subset(table.all,n_ind==0) # here we have 584 non-zero count for which we can engage with bpue estimation (at least for now the zeros are not dealt with, however in the future we will concentrate on finding a ways to extract information for some of theme_minimal
nob.species.df<-data.frame(species=unique(no.bycatch$Species),common=NA,class=NA)
nob.species.df$class<-tax_name(nob.species.df$species,get="class",db="itis")
#nob.species.df$common<-sci2comm(nob.species.df$species,db="ncbi")
no.bycatch$taxon<-"Fish"

no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Mammalia"&is.na(nob.species.df$class$class)==FALSE]]<-"Mammals"
no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Aves"&is.na(nob.species.df$class$class)==FALSE]]<-"Birds"
no.bycatch$taxon[no.bycatch$Species%in%nob.species.df$species[nob.species.df$class$class=="Reptilia"&is.na(nob.species.df$class$class)==FALSE]]<-"Reptilia"

no.bycatch$DaysAtSea<-ceiling(no.bycatch$DaysAtSea)

write.csv(no.bycatch,file="C:/Users/davlu/OneDrive - Danmarks Tekniske Universitet/WGBYC/wgbyc 2023/data/no_bycatch_detected.csv")

matchingmonitor<-read.csv("bpue_retain_taxa_issue.csv")

#QC1 can we fit a model and have enough monitoring
bpue.first<-read.csv("BPUE_est_2023.csv")
bpue.first<-subset(bpue.first,is.na(bpue)==FALSE)

df.QC1<-bpue.first[bpue.first$model=="n_ind ~ 1" & bpue.first$ base_model_heterogeneity==T,]
write.csv(df.QC1,"QC1_interceptonlymodel_heterogeneous.csv")

df.QC1replicates<-bpue.first[bpue.first$model=="n_ind ~ 1" & bpue.first$base_model_heterogeneity==F&bpue.first$replicates<5,]
write.csv(df.QC1replicates,"QC1_interceptonlymodel_homogeneous_lowreplicates.csv")

df.QC1replicates.table<-df.QC1replicates[,c(2,3,4,6,7,8)]
df.QC1replicates.table$bpue<-round(df.QC1replicates.table$bpue,3)
df.QC1replicates.table$lwr<-round(df.QC1replicates.table$lwr,4)
df.QC1replicates.table$upr<-round(df.QC1replicates.table$upr,4)
write.table(df.QC1replicates.table,"QC1_interceptonlymodel_homogeneous_lowreplicates_table.txt",dec = ".",sep=";")


df.QC2.QC3<-read.csv("BPUE_retain_QC2_QC3.csv")

########
#finally the taxonomic mismatch between species and observatoin method
df.taxa<-read.csv("bpue_retain_taxa_issue.csv")
df.retained<-df.taxa[df.taxa$taxa_bycatch_monitor_ok==FALSE&is.na(df.taxa$taxa_bycatch_monitor_ok)==FALSE&df.taxa$QC2=="green"&df.taxa$QC3=="green",]
write.csv(df.retained[,c(2:4)]y,file="bycatch_rejected_mismatch_species_method.csv")
