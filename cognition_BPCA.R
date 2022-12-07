# ASK 11/18/22
#+++ ADAPTED BY AP 12/7/22
# Note: this script will run BPCA on concatenated time 0 (baseline) and time 2 (year 2 followup) cognition data from ABCD. 
# The tasks included will be those that are measured at both time points: picvocab, flanker, pattern, reading, picture, RAVLT, LMT

#+++ use personal library instead of base sherlock
.libPaths('/home/users/apines/R/x86_64-pc-linux-gnu-library/4.1')

# Load in libraries
library(mvtnorm)
library(tableone)
library(parallel)
library(rstan)
library(loo)
library(gamm4)
library(Hmisc)
library(FactoMineR)
library(nFactors)
library(reshape2)
library(psych)
library(data.table)
library(mice)
library(abind)
library(cvTools)
library(modEvA)

# Load data:
nihtb<-read.table("/oak/stanford/groups/leanew1/users/apines/data/gp/abcd_tbss01.txt",header=TRUE)
lmt<-read.table("/oak/stanford/groups/leanew1/users/apines/data/gp/lmtp201.txt",header=TRUE)
ravlt<-read.table("/oak/stanford/groups/leanew1/users/apines/data/gp/abcd_ps01.txt",header=TRUE)
site_info <- readRDS("/oak/stanford/groups/leanew1/users/apines/data/gp/DEAP-siteID.rds")
colnames(site_info)[1]<-"subjectkey"

# Convert factors to numeric values
nihtb$nihtbx_flanker_uncorrected <- as.numeric(nihtb$nihtbx_flanker_uncorrected)
nihtb$nihtbx_picvocab_uncorrected <- as.numeric(nihtb$nihtbx_picvocab_uncorrected)
nihtb$nihtbx_pattern_uncorrected <- as.numeric(nihtb$nihtbx_pattern_uncorrected)
nihtb$nihtbx_picture_uncorrected <- as.numeric(nihtb$nihtbx_picture_uncorrected)
nihtb$nihtbx_reading_uncorrected <- as.numeric(nihtb$nihtbx_reading_uncorrected)
lmt$lmt_scr_perc_correct <- as.numeric(lmt$lmt_scr_perc_correct)
ravlt$pea_ravlt_sd_trial_i_tc <- as.numeric(ravlt$pea_ravlt_sd_trial_i_tc)
ravlt$pea_ravlt_sd_trial_ii_tc <- as.numeric(ravlt$pea_ravlt_sd_trial_ii_tc)
ravlt$pea_ravlt_sd_trial_iii_tc <- as.numeric(ravlt$pea_ravlt_sd_trial_iii_tc)
ravlt$pea_ravlt_sd_trial_iv_tc <- as.numeric(ravlt$pea_ravlt_sd_trial_iv_tc)
ravlt$pea_ravlt_sd_trial_v_tc <- as.numeric(ravlt$pea_ravlt_sd_trial_v_tc)


# Extract necessary cognition columns
nihtb_keep <- nihtb[,c("subjectkey","eventname","nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_reading_uncorrected","nihtbx_picture_uncorrected")]
lmt_keep <- lmt[,c("subjectkey","eventname","lmt_scr_perc_correct")]
ravlt_keep <- ravlt[,c("subjectkey","eventname","pea_ravlt_sd_trial_i_tc","pea_ravlt_sd_trial_ii_tc","pea_ravlt_sd_trial_iii_tc","pea_ravlt_sd_trial_iv_tc","pea_ravlt_sd_trial_v_tc")]


# for RAVLT we need to generate a summed score (code borrowed from Thompson et al. 2019):
ind_pea_ravlt = c(which(names(ravlt_keep)=="pea_ravlt_sd_trial_i_tc"),which(names(ravlt_keep)=="pea_ravlt_sd_trial_ii_tc"),
                  which(names(ravlt_keep)=="pea_ravlt_sd_trial_iii_tc"),which(names(ravlt_keep)=="pea_ravlt_sd_trial_iv_tc"),
                  which(names(ravlt_keep)=="pea_ravlt_sd_trial_v_tc")); names(ravlt_keep)[ind_pea_ravlt]; summary(ravlt_keep[,ind_pea_ravlt])
ravlt_keep$RAVLT <- apply(ravlt_keep[,ind_pea_ravlt],1,sum)


# Combine cognitive variables into a dataframe:
merged_cognition1<-merge(nihtb_keep,lmt_keep,by = c("subjectkey","eventname"))
merged_cognition_final<-merge(merged_cognition1,ravlt_keep,by = c("subjectkey","eventname"))

# Add in the family and site info:
site_info_keep <- site_info[site_info$event_name=="baseline_year_1_arm_1",c("subjectkey","abcd_site","rel_family_id")]
data_for_bpca <- merge(merged_cognition_final,site_info_keep,by="subjectkey")

# Subset just the complete cases
data_for_bpca <- data_for_bpca[complete.cases(data_for_bpca),]

# Subset participants with both T1 and T2 data
subIDlist<-unique(data_for_bpca$subjectkey)
subj_to_keep <- c()
for(subj in 1:length(unique(data_for_bpca$subjectkey))){
  if(length(data_for_bpca[data_for_bpca$subjectkey==subIDlist[subj],"eventname"])==2){
    subj_to_keep<-c(subj_to_keep,subIDlist[subj])
  }
}
data_for_bpca <- data_for_bpca[data_for_bpca$subjectkey %in% subj_to_keep,]

# Sanity check
dim(data_for_bpca[data_for_bpca$eventname=="baseline_year_1_arm_1",])[1]==dim(data_for_bpca[data_for_bpca$eventname=="2_year_follow_up_y_arm_1",])[1]

# Drop unneeded columns and adjust column names
data_for_bpca<-data_for_bpca[,c("subjectkey","nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_reading_uncorrected","nihtbx_picture_uncorrected","lmt_scr_perc_correct","RAVLT","abcd_site","rel_family_id")]
colnames(data_for_bpca)<-c("subjectkey","PicVocab","Flanker","Pattern","Reading","Picture","LMT","RAVLT","abcd_site","rel_family_id")

# Get indices for the nesting variables and cognitive variables, and check that they're correct
ind_nest = c(which(names(data_for_bpca)=="abcd_site"),which(names(data_for_bpca)=="rel_family_id"));names(data_for_bpca[,ind_nest])
ind_np = c(which(names(data_for_bpca)=="PicVocab"),which(names(data_for_bpca)=="Flanker"),
           which(names(data_for_bpca)=="Pattern"),which(names(data_for_bpca)=="Picture"),
           which(names(data_for_bpca)=="Reading"),which(names(data_for_bpca)=="RAVLT"),
           which(names(data_for_bpca)=="LMT")); names(data_for_bpca)[ind_np]


# Get numeric site and family variables (borrowed from Thompson et al. 2019 code):
data_for_bpca$site_num = as.numeric(substr(data_for_bpca$abcd_site,5,6))
data_for_bpca$fam_num = 0
ind=0
for(i in sort(unique(data_for_bpca$rel_family_id))){
  ind = ind+1
  data_for_bpca$fam_num[data_for_bpca$rel_family_id==i & !is.na(data_for_bpca$rel_family_id)] = ind
}


# Run a simple PCA 
ind_Y = ind_np; names(data_for_bpca)[ind_Y]
Y = as.matrix(scale(data_for_bpca[complete.cases(data_for_bpca[,c(ind_Y)]),ind_Y]))
ev = eigen(cor(Y))
ap = parallel(subject=nrow(Y),var=ncol(Y),rep=100,cent=.05)
nS = nScree(x=ev$values,aparallel=ap$eigen$qevpea)
plotnScree(nS)
ncomp = 3
#y.pca = psych::principal(Y, rotate="promax", nfactors=ncomp, scores=TRUE)
#y.pca$loadings
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
y.pca$loadings
# ASK checked the initial version of this on 11/18/22 and verified that the top 3 tasks loading onto the top 3 PC's are the same as for the Thompson et al. 2019 factor scores



# recode family and site variables according to code from Thompson et al. 2019 
data1 = data_for_bpca[complete.cases(data_for_bpca[,c(ind_Y)]),]; names(data1); dim(data1)
site_num = rep(NA,length(data1$site_num))
for(s in 1:length(unique(data1$site_num))){
  site_num[data1$site_num == unique(data1$site_num)[s]] = s
}
data1$site_num = site_num; rm(site_num)

#+++ AP - ask for clarification on this
data1$id_fam = 0
data1$fam_size = 0
ind=0
for(s in 1:length(unique(data1$site_num))){
  data_s = data1[data1$site_num == s, ]
  for(f in 1:length(unique(data_s$rel_family_id))){
    data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = 
      sum(data_s$rel_family_id == unique(data_s$rel_family_id)[f])
    if(sum(data_s$fam_size[data_s$rel_family_id == unique(data_s$rel_family_id)[f]])>1){	
      ind=ind+1	
      data_s$id_fam[data_s$rel_family_id == unique(data_s$rel_family_id)[f]] = ind
    }	
  }
  data1[data1$site_num == s, ] = data_s
}
data1 = data1[order(data1$site_num,data1$id_fam),]

Site = data1$site_num
Fam = data1$id_fam
ns = length(unique(Site))
ns_s = rep(0,ns)
for(s in 1:ns){
  ns_s[s] = sum(Site==s)
}
nf_s = rep(0,ns)
for(s in 1:ns){
  nf_s[s] = length(unique(Fam[Site==s]))-1
}
nf = sum(nf_s)


# Setup for BPCA (code and model file from Thompson et al. 2019)
Y = (as.matrix(scale(data1[,ind_Y]))); summary(Y)
N = nrow(Y)
P = ncol(Y)
Nsamples = 1000
Nchains = 3
model_file = "/oak/stanford/groups/leanew1/users/apines/data/gp/bppca.stan"
smod = stan_model(model_file)


nf=4829 # Arielle added this line -- necessary to make sure the number of families with more than 1 individual is correct (note that since every individual is repeated twice here, indivdiuals are considered "siblings" with themselves)


# Run BPCA
D_max = 5
sa.list = list()
log_lik.list = list()
looic.list = list()
for(d in 1:D_max){
  print(d)
  pca_data <- list(Y = Y, N = N, P = P, D = d, Fam =  Fam, Site = Site, ns = ns, nf = nf)
  set.seed(314)
  sa.list[[d]] = sampling(smod, data= pca_data, iter=Nsamples, chains=Nchains,init="random")
  #log_lik.list[[d]] <- extract_log_lik(sa.list[[d]])
  log_lik.list[[d]] <- extract(sa.list[[d]],"log_lik_marg")[[1]]
  looic.list[[d]] = loo(log_lik.list[[d]])
  save(sa.list,log_lik.list,looic.list,file="/oak/stanford/groups/leanew1/users/apines/data/gp/bppca_results.RData")
  print("###############################")
}	




