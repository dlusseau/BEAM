
zeros_prob_vect<-function(n_ind=1L,daysatsea=daysatsea,rarity=veryrare,bootstrap=1000L) {
  mean(rbinom(bootstrap, size = daysatsea * n_ind, prob = rarity) == 0L) #speed up the initial approach 
}


zero_assessments<-function(bpue2,bpue1) {

### thresholds decided at ICES WKBBEAM Dec 2025
veryrare<-0.001
p<-0.01
################################################


zero.candidates<-bpue2[bpue2$model=="none",]

zero.candidates[, daysatsea := round(daysatsea, 0)]
zero.candidates <- zero.candidates[daysatsea != 0]

zero.candidates$p_veryrare<-unlist(lapply(zero.candidates$daysatsea,function(x) zeros_prob_vect(n_ind=1L,daysatsea=x,rarity=veryrare,bootstrap=1000L)))
zero.candidates$upper_parametric<-1.96 * sqrt(2 / (zero.candidates$daysatsea+2)^2)

keys <- c("ecoregion", "metierl4", "species")
zc_ok <- unique(zero.candidates[p_veryrare <= p, ..keys])
zc_upper <- zero.candidates[, .(upper_parametric = upper_parametric), by = keys]

bpue2[zc_ok, on = keys, bpue := fifelse(is.na(bpue), 0, bpue)]
bpue2[zc_ok, on = keys, lwr := fifelse(is.na(bpue), 0, bpue)]
bpue2[zc_upper, on = keys, upr := i.upper_parametric]

bpue1[zc_ok, on = keys, bpue := fifelse(is.na(bpue), 0, bpue)]
bpue1[zc_ok, on = keys, lwr := fifelse(is.na(bpue), 0, bpue)]
bpue1[zc_upper, on = keys, upr := i.upper_parametric]


fwrite(bpue2, file = "data/bpue2.csv", sep = ";",na="NA")
fwrite(bpue1, file = "data/bpue1.csv", sep = ";",na="NA")
}
