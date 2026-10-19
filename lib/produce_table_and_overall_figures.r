###produce table to print and main total bycatch estimate figures

table2print<-total_bycatch
names(table2print)[5:7]<-c("tbfinal.mean","tbfinal.lwr","tbfinal.upr")

table2print<-merge(table2print,bpue1,by=c("ecoregion","metierl4","species"),all=FALSE)
table2print<-merge(table2print,summaryyear,by=c("ecoregion","metierl4","species"),all.x=TRUE)
table2print<-merge(table2print,summaryall,by=c("ecoregion","metierl4","species"),all.x=TRUE)

###MP 07/01/26: Added a few steps to recover taxa information that is not in total_bycatch but required in lines below.
common_names <- fread("data/commonnames2024.csv")
common_names$common <- tolower(common_names$common)
common_names$taxon <- tolower(common_names$taxon)
table2print[obs3,on = c("species"), taxa := i.taxon] #Retreiving taxa info from obs3 instead of common_names table so the elasmobranchs are distinct from other fish
table2print[common_names,on = c("species"), common := i.common]
###

#dat$n_ind / dat$daysatsea

table.print<-table2print[,c("ecoregion","metierl4","taxa","species","common","n_ind","daysatsea","fishing_effort","n_ind_all","daysatsea_all","model.x","bpue","lwr","upr","tbfinal.mean","tbfinal.lwr","tbfinal.upr")]

table.print$interannual<-"none apparent"
table.print$interannual[grep("year",table.print$model)]<-"there is between-year variability in BPUE"

#07/01/2026: Removed capital letters in variable names so grep() works properly
table.print$key.representability<-"a constant BPUE appears to be representative"
table.print$key.representability[grep("metierl5",table.print$model)]<-"there is between-metier level 5 variability in BPUE"
table.print$key.representability[grep("vessellength_group",table.print$model)]<-"there is between-vessel length category variability in BPUE"
table.print$key.representability[grep("areacode",table.print$model)]<-"there is spatial variability in BPUE"
table.print$key.representability[grep("country",table.print$model)]<-"there is spatial variability in BPUE"
table.print$key.representability[grep("monitoringmethod",table.print$model)]<-"there is variability in BPUE depending on monitoring methods"
table.print$key.representability[grep("samplingprotocol",table.print$model)]<-"there is variability in BPUE depending on sampling protocols"

table.print$n_ind<-round(table.print$n_ind,0)
table.print$daysatsea<-round(table.print$daysatsea,0)
table.print$n_ind_all<-round(table.print$n_ind_all,0)
table.print$daysatsea_all<-round(table.print$daysatsea_all,0)
table.print$fishing_effort<-round(table.print$fishing_effort,0)
table.print$bpue<-round(table.print$bpue,6)
table.print$lwr<-round(table.print$lwr,7)
table.print$upr<-round(table.print$upr,7)
table.print$tbfinal.mean<-round(table.print$tbfinal.mean,0)
table.print$tbfinal.lwr<-round(table.print$tbfinal.lwr,1)
table.print$tbfinal.upr<-round(table.print$tbfinal.upr,1)
table.print$taxa<-factor(table.print$taxa,levels=c("mammals","seabirds","turtles","elasmobranchs","fish"))

colnames(table.print)<-c("Ecoregion", "metier L4", "Taxon", "Species", "Common name", 
							"# individuals 2025", "monitoring effort (DaS) 2025", 
							"Fishing effort (DaS) 2025","# individuals 2017-2025", "monitoring effort (DaS) 2017-2025", 
							"BPUE model", "BPUE", "lower", "upper",
							"total bycatch 2025", "TB lower", "TB upper", "interannual", "key variability in BPUE")
							
table.print$"# individuals 2017-2025"<-round(table.print$"# individuals 2017-2025",0)
table.print$"# individuals 2025"<-round(table.print$"# individuals 2025",0)

table.print$interannual[table.print$interannual=="none apparent"]<-"none"
table.print$interannual[table.print$interannual!="none"]<-"present"
table.print$BPUE<-round(table.print$BPUE,4)
table.print$lower<-round(table.print$lower,4)
table.print$upper<-round(table.print$upper,4)


all3<-aggregate(daysatseaf~ecoregion+metierl4,data=all2,function(x) sum(x,na.rm=TRUE))

temp<-merge(table.print,all3,by.x=c("Ecoregion","metier L4"),by.y=c("ecoregion","metierl4"),all.x=TRUE)
table.print$"Fishing effort (DaS) 2025"<-round(temp$daysatseaf,0)

temp<-merge(table.print,bpue2,by.x=c("Ecoregion","metier L4","Species"),by.y=c("ecoregion","metierl4","species"),all.x=TRUE)
table.print$"Monitoring coverage"<-round(temp$monitoring_coverage,5)



fwrite(table.print,file="results/bpue_table_print.csv",na="NA")


##### render the table for report
ft <- flextable(table.print[order(table.print$Ecoregion,table.print$Taxon,table.print$"metier L4",table.print$Species),])
ft <-theme_vanilla(ft)
set_table_properties(ft, width = 1, layout = "autofit")
# save_as_docx(
  # ft, 
  # path = "results/BPUE_table1.docx")

save_as_html(
  ft, 
  path = "results/BPUE_table1.html")

complete_tb2print<-complete_tb
names(complete_tb2print)[5:7]<-c("tbfinal.mean","tbfinal.lwr","tbfinal.upr")

complete_tb2print<-merge(complete_tb2print,bpue1,c("ecoregion","metierl4","species"),all.x=TRUE)
complete_tb2print<-merge(complete_tb2print,summaryyear,c("ecoregion","metierl4","species"),all.x=TRUE)
complete_tb2print<-merge(complete_tb2print,summaryall,c("ecoregion","metierl4","species"),all.x=TRUE)

complete_tb2print[common_names,on = c("species"), taxa := i.taxon]
complete_tb2print[common_names,on = c("species"), common := i.common]

CTB_sum = complete_tb2print[
  ,
  .(n_ind = sum(n_ind),
	monitoring_DaS=sum(daysatsea),
	n_ind_all = sum(n_ind_all),
	monitoring_DaS_all=sum(daysatsea_all),
	fishing_DaS=sum(fishing_effort),
	TB=sum(tbfinal.mean),
	lower=sum(tbfinal.lwr),
	upper=sum(tbfinal.upr),
	taxa=unique(taxa),
	metier=length(common),
	common=unique(common)
	),
  keyby = .(ecoregion, species)]


CTB_sum$n_ind<-round(CTB_sum$n_ind,0)
CTB_sum$monitoring_DaS<-round(CTB_sum$monitoring_DaS,0)
CTB_sum$n_ind_all<-round(CTB_sum$n_ind_all,0)
CTB_sum$monitoring_DaS_all<-round(CTB_sum$monitoring_DaS_all,0)
CTB_sum$fishing_DaS<-round(CTB_sum$fishing_DaS,0)
CTB_sum$TB<-round(CTB_sum$TB,0)
CTB_sum$lower<-round(CTB_sum$lower,1)
CTB_sum$upper<-round(CTB_sum$upper,1)
CTB_sum$taxa<-factor(CTB_sum$taxa,levels=c("mammals","seabirds","turtles","elasmobranchs","fish"))


colnames(CTB_sum)<-c("Ecoregion", "species","# individuals 2025","monitoring effort (DaS) 2025",
"# individuals 2017-2025","monitoring effort (DaS) 2017-2025", 
"Fishing effort (DaS) 2025","total bycatch 2025", "TB lower", "TB upper","Taxon", "# metiers","Common name")

CTB_sum<-CTB_sum[,c("Ecoregion","Taxon","Common name","# individuals 2025","monitoring effort (DaS) 2025",
					"# individuals 2017-2025","monitoring effort (DaS) 2017-2025", 
					"Fishing effort (DaS) 2025","total bycatch 2025", "TB lower", "TB upper")]



ftc <- flextable(CTB_sum[order(CTB_sum$Ecoregion,CTB_sum$Taxon,CTB_sum$"Common name"),])
ftc <-theme_vanilla(ftc)
set_table_properties(ftc, width = 1, layout = "autofit")
save_as_docx(
  ftc, 
  path = "results/Total_bycatch_table1.docx")

save_as_html(
  ftc, 
  path = "results/Total_bycatch_table1.html")


#####

#MP 13/01/2026 Redefined priority1, 2 and 3 using code form previous versions.
priority1<-c("phocoena phocoena","delphinus delphis","monachus monachus")
priority2<-"turtles"
priority3<-c("angel shark", "common skate", "guitarfish", "maltese
skate", "great white shark", "sand tiger shark", "smalltooth sand tiger shark",
"spiny butterfly ray", "sturgeon", "balearic shearwater","sterlet")

priority3<-paste(priority3,collapse="|")

priorty1_tb<-table.print[table.print$Species%in%priority1&!is.na(table.print$"total bycatch 2025"),]
priorty2_tb<-table.print[table.print$Taxon=="turtles"&!is.na(table.print$"total bycatch 2025"),]
priorty3_tb<-table.print[grep(priority3,tolower(table.print$"Common name")),]
priorty3_tb<-priorty3_tb[!is.na(priorty3_tb$"total bycatch 2025"),]

priority_list<-rbind(priorty1_tb,priorty2_tb,priorty3_tb)
priority_list<-priority_list[,-c(9,16)]

ftp <- flextable(priority_list[order(priority_list$Ecoregion,priority_list$Taxon,priority_list$"metier L4",priority_list$Species),])
ftp <-theme_vanilla(ftp)
set_table_properties(ftp, width = 1, layout = "autofit")
save_as_docx(
  ftp, 
  path = "results/priority_table1.docx")



#################################################################################################################
#### rate of achievement

achieved<-subset(total_bycatch,!is.na(tot_mean))
achieved$label<-paste(str_to_title(achieved$ecoregion),toupper(achieved$metierl4),achieved$common,sep=" and ") #Capitalized metier acronyms and moved species at the end of the label to match order in the plot captions below

#MP 13/01/2026 
#Adding a line to grab taxa information again, as the column does not currently exist in total_bycatch
achieved[obs3,on = c("species"), taxa := i.taxon] #Retreiving taxa info from obs3 instead of common_names table so the elasmobranchs are distinct from other fish

mammals<-ggplot(subset(achieved,taxa=="mammals"
								#& 
								#label!="Short-beaked Common Dolphin and bay of biscay and the iberian coast and ptb" &
								#label!="Harbor Porpoise and bay of biscay and the iberian coast and gtr" &
								#label!="Gray Seal and celtic seas and otm" &
								#label!="Gray Seal and baltic sea and gtr"
								), #removing the few with unsensible CI
  aes(x = label, y=(tot_mean)),colour="white") +
  geom_crossbar(aes(ymin = (tot_lwr), ymax = (tot_upr)), width = 0.5, fill = "light blue",fatten=0.5) + #increased fatten argument to 0.5 so the mean line is visible on the plot bars.
    
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (individuals) - note that the axis is on a log scale")+
  theme_minimal()+
  scale_y_log10()

birds<-ggplot(subset(achieved,taxa=="seabirds"
							  #& 
							  #label!="European Shag and adriatic sea and otb" &
							  #label!="Northern Fulmar and celtic seas and gns" &
							  #label!="Herring Gull and bay of biscay and the iberian coast and otb" &
							  #label!="Northern Gannet and greater north sea and gns" &
							  #label!="European Shag and adriatic sea and otb" 
							  ), 
   aes(x = label, y=(tot_mean)),colour="white") +
  geom_crossbar(aes(ymin = (tot_lwr), ymax = (tot_upr)), width = 0.5, fill = "light blue",fatten=0.5) +
   
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (individuals) - note that the axis is on a log scale")+
  theme_minimal()+
  scale_y_continuous(trans = "log10")


reptiles<-ggplot(subset(achieved,taxa=="turtles"), 
    aes(x = label, y=(tot_mean)),colour="white") +
  geom_crossbar(aes(ymin = (tot_lwr), ymax = (tot_upr)), width = 0.5, fill = "light blue",fatten=0.5) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (individuals) - note that the axis is on a log scale")+
  theme_minimal()+
  scale_y_continuous(trans = "log10")



elasmo<-ggplot(subset(achieved,taxa=="elasmobranchs"&tot_upr<1000000&tot_lwr>0.001), 
    aes(x = label, y=(tot_mean)),colour="white") +
  geom_crossbar(aes(ymin = (tot_lwr), ymax = (tot_upr)), width = 0.5, fill = "light blue",fatten=0.5) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (individuals) - note that the axis is on a log scale")+
  theme_minimal()+
  scale_y_continuous(trans = "log10")


 fish<-ggplot(subset(achieved,taxa=="fish"&tot_upr<100000), 
   aes(x = label, y=(tot_mean)),colour="white") +
  geom_crossbar(aes(ymin = (tot_lwr), ymax = (tot_upr)), width = 0.5, fill = "light blue",fatten=0.5) +
  coord_flip() +
  xlab("Ecoregion, Metier level 4, Species") + 
  ylab("Total Bycatch (individuals) - note that the axis is on a log scale")+
  theme_minimal()+
  scale_y_continuous(trans = "log10")


	
png("results/totalbycatch_mammals.png",width=20,height=35,units="cm",res=200)
mammals
dev.off()


png("results/totalbycatch_turtles.png",width=20,height=35,units="cm",res=200)
reptiles
dev.off()


png("results/totalbycatch_birds.png",width=20,height=35,units="cm",res=200)
birds
dev.off()


png("results/totalbycatch_elasmo.png",width=20,height=35,units="cm",res=200)
elasmo
dev.off()


png("results/totalbycatch_fish.png",width=20,height=35,units="cm",res=200)
fish
dev.off()




