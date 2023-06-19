library(ppcor)
library(mgcv)
library(pammtools)
library(cowplot)
library(patchwork)
library(gratia)
library(scales)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(visreg)
library(ggExtra)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age','sex','income')]
# get length of df for later
lenDF=dim(masterdf)[1]
# will need to get full and reduced models for each boot, as well as a null distribution
# in addition to derivatives and fits, save F values of interaction + null distribution F values
# interactions to be tested: p*sex, p*poverty, p*sex*poverty
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r=as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r=as.numeric(masterdf$cbcl_scr_syn_external_r)
# sex (g in for x to get co-pilot around its NSFW filter) and poverty to factors
masterdf$seg<-as.ordered(masterdf$sex)
masterdf$poverty=0
masterdf$income<-as.numeric(masterdf$income)
# note that poverty is defined as income < 5: https://collection3165.readthedocs.io/en/stable/recommendations/#2-the-bids-participants-files-and-matched-groups
masterdf$poverty[masterdf$income<5]=1
masterdf$poverty=as.ordered(masterdf$poverty)
### initialize cross-boot vectors
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
iMaxVal=max(masterdf$cbcl_scr_syn_internal_r)
emaxVal=max(masterdf$cbcl_scr_syn_external_r)
F_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
M_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
P_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
R_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
PF_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
PM_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
RF_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
RM_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
F_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
M_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
P_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
R_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
PF_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
PM_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
RF_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
RM_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
F_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
M_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
P_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
R_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
PF_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
PM_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
RF_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
RM_extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
# predicted values: set to maximum value for ncol
F_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
M_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
P_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
R_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
PF_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
PM_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
RF_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
RM_pFit=matrix(0,nrow=10000,ncol=pMaxVal)
F_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
M_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
P_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
R_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
PF_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
PM_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
RF_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
RM_intFit=matrix(0,nrow=10000,ncol=iMaxVal)
F_extFit=matrix(0,nrow=10000,ncol=emaxVal)
M_extFit=matrix(0,nrow=10000,ncol=emaxVal)
P_extFit=matrix(0,nrow=10000,ncol=emaxVal)
R_extFit=matrix(0,nrow=10000,ncol=emaxVal)
PF_extFit=matrix(0,nrow=10000,ncol=emaxVal)
PM_extFit=matrix(0,nrow=10000,ncol=emaxVal)
RF_extFit=matrix(0,nrow=10000,ncol=emaxVal)
RM_extFit=matrix(0,nrow=10000,ncol=emaxVal)
# F statistic vector for each interaction
Fseg_p=rep(0,10000)
Fpov_p=rep(0,10000)
Fsegpov_p=rep(0,10000)
Fseg_int=rep(0,10000)
Fpov_int=rep(0,10000)
Fsegpov_int=rep(0,10000)
Fseg_ext=rep(0,10000)
Fpov_ext=rep(0,10000)
Fsegpov_ext=rep(0,10000)
# AIC for each model
AICreduced_p=rep(0,10000)
AICseg_p=rep(0,10000)
AICpov_p=rep(0,10000)
AICsegpov_p=rep(0,10000)
AICreduced_int=rep(0,10000)
AICseg_int=rep(0,10000)
AICpov_int=rep(0,10000)
AICsegpov_int=rep(0,10000)
AICreduced_ext=rep(0,10000)
AICseg_ext=rep(0,10000)
AICpov_ext=rep(0,10000)
AICsegpov_ext=rep(0,10000)
# comparison AICs: includes main effect terms but not interaction terms
AICseg_p_noIntrxn=rep(0,10000)
AICpov_p_noIntrxn=rep(0,10000)
AICsegpov_p_noIntrxn=rep(0,10000)
AICseg_int_noIntrxn=rep(0,10000)
AICpov_int_noIntrxn=rep(0,10000)
AICsegpov_int_noIntrxn=rep(0,10000)
AICseg_ext_noIntrxn=rep(0,10000)
AICpov_ext_noIntrxn=rep(0,10000)
AICsegpov_ext_noIntrxn=rep(0,10000)
# now equivalent vectors for null models
Fseg_null_p=rep(0,10000)
Fpov_null_p=rep(0,10000)
Fsegpov_null_p=rep(0,10000)
Fseg_null_int=rep(0,10000)
Fpov_null_int=rep(0,10000)
Fsegpov_null_int=rep(0,10000)
Fseg_null_ext=rep(0,10000)
Fpov_null_ext=rep(0,10000)
Fsegpov_null_ext=rep(0,10000)
# AIC
AICseg_null_p=rep(0,10000)
AICpov_null_p=rep(0,10000)
AICsegpov_null_p=rep(0,10000)
AICseg_null_int=rep(0,10000)
AICpov_null_int=rep(0,10000)
AICsegpov_null_int=rep(0,10000)
AICseg_null_ext=rep(0,10000)
AICpov_null_ext=rep(0,10000)
AICsegpov_null_ext=rep(0,10000)
# and three lil' vectors just to track max p int and ext over iterations
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
# and modular fit for pov vs. psuedopov 2 (matched #) fits
povFit=matrix(0,nrow=10000,ncol=pMaxVal)
pseudopovFit=matrix(0,nrow=10000,ncol=pMaxVal)
#################################################
##################### loop over manual bootstraps
#################################################
for (b in 1:10000){
	print(b)
	# get subjects to include in this bootstrap
	BootSubjs=sample(subjs,numSubjs,replace=T)
	### inefficient but interpretable loop
	# Create an empty dataframe to store the resampled observations
	bootSamp <- data.frame()
	for (j in 1:length(BootSubjs)){
		subject_obs <- masterdf[masterdf$subjectkey == BootSubjs[j], ]
		bootSamp <- rbind(bootSamp, subject_obs)
	}
	##############################
	# get max for this iteration
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
	bimax=max(bootSamp$cbcl_scr_syn_internal_r)
	bemax=max(bootSamp$cbcl_scr_syn_external_r)
	# last pre-step: need to create NULL variables. Use count of F and count of Poverty as count of membership to psuedo-groups, but randomly distribute membership
	# get count F in this boot
	Fcount=length(which(bootSamp$sex=="F"))
	# get count poverty in this boot
	Povcount=length(which(bootSamp$poverty==1))
	# randomly assign Fcount people to pseudoseg
	bootSamp$psuedoseg=0
	bootSamp$psuedoseg[sample(1:dim(bootSamp)[1],Fcount)]=1
	# randomly assign Povcount people to pseudopoverty
	bootSamp$psuedopoverty=0
	bootSamp$psuedopoverty[sample(1:dim(bootSamp)[1],Povcount)]=1
	# make a df that is just the same number of kids in poverty, but actually from the non-poverty group. The point is to see if we can recover g~p slope in a an equivalent number of non-poverty kids
	bootSamp$pseudopoverty2=0
	bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
	
	#
	######## I FIT MODELS
	#
	
	#### g as response variable: REDUCED MODELS
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	#### FULL MODELS
	pgAge_seg<-bam(g~s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg+s(interview_age),data=bootSamp)
	intgAge_seg<-bam(g~s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_internal_r)+seg+s(interview_age),data=bootSamp)
	extgAge_seg<-bam(g~s(cbcl_scr_syn_external_r,by=seg)+s(cbcl_scr_syn_external_r)+seg+s(interview_age),data=bootSamp)
	pgAge_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r)+poverty+s(interview_age),data=bootSamp)
	intgAge_pov<-bam(g~s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_internal_r)+poverty+s(interview_age),data=bootSamp)
	extgAge_pov<-bam(g~s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_external_r)+poverty+s(interview_age),data=bootSamp)
	# comparison models to evaluate AIC gain from interactions specifically. Should include main effects but minus interaction of interest
	pgAge_seg_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r)+seg+s(interview_age),data=bootSamp)
	intgAge_seg_noIntrxn=bam(g~s(cbcl_scr_syn_internal_r)+seg+s(interview_age),data=bootSamp)
	extgAge_seg_noIntrxn=bam(g~s(cbcl_scr_syn_external_r)+seg+s(interview_age),data=bootSamp)
	pgAge_pov_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r)+poverty+s(interview_age),data=bootSamp)
	intgAge_pov_noIntrxn=bam(g~s(cbcl_scr_syn_internal_r)+poverty+s(interview_age),data=bootSamp)
	extgAge_pov_noIntrxn=bam(g~s(cbcl_scr_syn_external_r)+poverty+s(interview_age),data=bootSamp)
	pgAge_seg_pov_noIntrxn<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg*poverty+s(interview_age),data=bootSamp)
	intgAge_seg_pov_noIntrxn<-bam(g~s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_internal_r)+seg*poverty+s(interview_age),data=bootSamp)
	extgAge_seg_pov_noIntrxn<-bam(g~s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_external_r,by=seg)+s(cbcl_scr_syn_external_r)+seg*poverty+s(interview_age),data=bootSamp)
	#### VERY FULL MODELS
	pgAge_seg_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_totprob_r, by = interaction(seg, poverty)),data=bootSamp)
	intgAge_seg_pov<-bam(g~s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_internal_r)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_internal_r, by = interaction(seg, poverty)),data=bootSamp)
	extgAge_seg_pov<-bam(g~s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_external_r,by=seg)+s(cbcl_scr_syn_external_r)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_external_r, by = interaction(seg, poverty)),data=bootSamp)
	# fit null models for use later (no reduced, reduced is the same as real reduced)
	#### FULL NULL MODELS
	pgAge_seg_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedoseg)+s(cbcl_scr_syn_totprob_r)+psuedoseg+s(interview_age),data=bootSamp)
	intgAge_seg_n<-bam(g~s(cbcl_scr_syn_internal_r,by=psuedoseg)+s(cbcl_scr_syn_internal_r)+psuedoseg+s(interview_age),data=bootSamp)
	extgAge_seg_n<-bam(g~s(cbcl_scr_syn_external_r,by=psuedoseg)+s(cbcl_scr_syn_external_r)+psuedoseg+s(interview_age),data=bootSamp)
	pgAge_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty)+s(cbcl_scr_syn_totprob_r)+psuedopoverty+s(interview_age),data=bootSamp)
	intgAge_pov_n<-bam(g~s(cbcl_scr_syn_internal_r,by=psuedopoverty)+s(cbcl_scr_syn_internal_r)+psuedopoverty+s(interview_age),data=bootSamp)
	extgAge_pov_n<-bam(g~s(cbcl_scr_syn_external_r,by=psuedopoverty)+s(cbcl_scr_syn_external_r)+psuedopoverty+s(interview_age),data=bootSamp)
	#### VERY FULL NULL MODELS
	pgAge_seg_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty)+s(cbcl_scr_syn_totprob_r,by=psuedoseg)+s(cbcl_scr_syn_totprob_r)+psuedoseg*psuedopoverty+s(interview_age)+s(cbcl_scr_syn_totprob_r, by = interaction(psuedoseg, psuedopoverty)),data=bootSamp)
	intgAge_seg_pov_n<-bam(g~s(cbcl_scr_syn_internal_r,by=psuedopoverty)+s(cbcl_scr_syn_internal_r,by=psuedoseg)+s(cbcl_scr_syn_internal_r)+psuedoseg*psuedopoverty+s(interview_age)+s(cbcl_scr_syn_internal_r, by = interaction(psuedoseg, psuedopoverty)),data=bootSamp)
	extgAge_seg_pov_n<-bam(g~s(cbcl_scr_syn_external_r,by=psuedopoverty)+s(cbcl_scr_syn_external_r,by=psuedoseg)+s(cbcl_scr_syn_external_r)+psuedoseg*psuedopoverty+s(interview_age)+s(cbcl_scr_syn_external_r, by = interaction(psuedoseg, psuedopoverty)),data=bootSamp)
	
	#
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#
	
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(1:bpmax)
	eachIntcount=seq(1:bimax)
	eachExtcount=seq(1:bemax)
	# set age to to median for predict df
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax))
	predictDFint=data.frame(eachIntcount,rep(median(bootSamp$interview_age),bimax))
	predictDFext=data.frame(eachExtcount,rep(median(bootSamp$interview_age),bemax))
	# set predict df to interacting factors of interest (seg, poverty)
	predictDF_segp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax),rep("F",bpmax))
	predictDF_povp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax),rep("1",bpmax))
	predictDF_segpovp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax),rep("F",bpmax),rep("1",bpmax))
	predictDF_segint=data.frame(eachIntcount,rep(median(bootSamp$interview_age),bimax),rep("F",bimax))
	predictDF_povint=data.frame(eachIntcount,rep(median(bootSamp$interview_age),bimax),rep("1",bimax))
	predictDF_segpovint=data.frame(eachIntcount,rep(median(bootSamp$interview_age),bimax),rep("F",bimax),rep("1",bimax))
	predictDF_segext=data.frame(eachExtcount,rep(median(bootSamp$interview_age),bemax),rep("F",bemax))
	predictDF_povext=data.frame(eachExtcount,rep(median(bootSamp$interview_age),bemax),rep("1",bemax))
	predictDF_segpovext=data.frame(eachExtcount,rep(median(bootSamp$interview_age),bemax),rep("F",bemax),rep("1",bemax))
	# set colnames so predict can work
	colnames(predictDF_segp)=c('cbcl_scr_syn_totprob_r','interview_age','seg')
	colnames(predictDF_povp)=c('cbcl_scr_syn_totprob_r','interview_age','poverty')
	colnames(predictDF_segpovp)=c('cbcl_scr_syn_totprob_r','interview_age','seg','poverty')
	colnames(predictDF_segint)=c('cbcl_scr_syn_internal_r','interview_age','seg')
	colnames(predictDF_povint)=c('cbcl_scr_syn_internal_r','interview_age','poverty')
	colnames(predictDF_segpovint)=c('cbcl_scr_syn_internal_r','interview_age','seg','poverty')
	colnames(predictDF_segext)=c('cbcl_scr_syn_external_r','interview_age','seg')
	colnames(predictDF_povext)=c('cbcl_scr_syn_external_r','interview_age','poverty')
	colnames(predictDF_segpovext)=c('cbcl_scr_syn_external_r','interview_age','seg','poverty')
	# predict girl
	forfitp_F=predict(pgAge_seg,predictDF_segp)
	forfitint_F=predict(intgAge_seg,predictDF_segint)
	forfitext_F=predict(extgAge_seg,predictDF_segext)
	forDerivp_F=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=bpmax)
	forDerivint_F=derivatives(intgAge_seg,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segint,n=bimax)
	forDerivext_F=derivatives(extgAge_seg,term="s(cbcl_scr_syn_external_r)",data=predictDF_segext,n=bemax)
	# predict boy
	predictDF_segp$seg=rep("M",bpmax)
	forfitp_M=predict(pgAge_seg,predictDF_segp)
	forDerivp_M=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=bpmax)
	predictDF_segint$seg=rep("M",bimax)
	forfitint_M=predict(intgAge_seg,predictDF_segint)
	forDerivint_M=derivatives(intgAge_seg,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segint,n=bimax)
	predictDF_segext$seg=rep("M",bemax)
	forfitext_M=predict(extgAge_seg,predictDF_segext)
	forDerivext_M=derivatives(extgAge_seg,term="s(cbcl_scr_syn_external_r)",data=predictDF_segext,n=bemax)
	# predict poverty
	forfitp_P=predict(pgAge_pov,predictDF_povp)
	forDerivp_P=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=bpmax)
	forfitint_P=predict(intgAge_pov,predictDF_povint)
	forDerivint_P=derivatives(intgAge_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_povint,n=bimax)
	forfitext_P=predict(extgAge_pov,predictDF_povext)
	forDerivext_P=derivatives(extgAge_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_povext,n=bemax)
	# predict rich
	predictDF_povp$poverty=rep("0",bpmax)
	forfitp_R=predict(pgAge_pov,predictDF_povp)
	forDerivp_R=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=bpmax)
	predictDF_povint$poverty=rep("0",bimax)
	forfitint_R=predict(intgAge_pov,predictDF_povint)
	forDerivint_R=derivatives(intgAge_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_povint,n=bimax)
	predictDF_povext$poverty=rep("0",bemax)
	forfitext_R=predict(extgAge_pov,predictDF_povext)
	forDerivext_R=derivatives(extgAge_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_povext,n=bemax)
	## predict across triple interaction, funtime
	# poverty female
	forfitp_PF=predict(pgAge_seg_pov,predictDF_segpovp)
	forDerivp_PF=derivatives(pgAge_seg_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segpovp,n=bpmax)
	forfitint_PF=predict(intgAge_seg_pov,predictDF_segpovint)
	forDerivint_PF=derivatives(intgAge_seg_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segpovint,n=bimax)
	forfitext_PF=predict(extgAge_seg_pov,predictDF_segpovext)
	forDerivext_PF=derivatives(extgAge_seg_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_segpovext,n=bemax)
	# poverty male
	predictDF_segpovp$seg=rep("M",bpmax)
	forfitp_PM=predict(pgAge_seg_pov,predictDF_segpovp)
	forDerivp_PM=derivatives(pgAge_seg_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segpovp,n=bpmax)
	predictDF_segpovint$seg=rep("M",bimax)
	forfitint_PM=predict(intgAge_seg_pov,predictDF_segpovint)
	forDerivint_PM=derivatives(intgAge_seg_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segpovint,n=bimax)
	predictDF_segpovext$seg=rep("M",bemax)
	forfitext_PM=predict(extgAge_seg_pov,predictDF_segpovext)
	forDerivext_PM=derivatives(extgAge_seg_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_segpovext,n=bemax)
	# rich male
	predictDF_segpovp$poverty=rep("0",bpmax)
	forfitp_RM=predict(pgAge_seg_pov,predictDF_segpovp)
	forDerivp_RM=derivatives(pgAge_seg_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segpovp,n=bpmax)
	predictDF_segpovint$poverty=rep("0",bimax)
	forfitint_RM=predict(intgAge_seg_pov,predictDF_segpovint)
	forDerivint_RM=derivatives(intgAge_seg_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segpovint,n=bimax)
	predictDF_segpovext$poverty=rep("0",bemax)
	forfitext_RM=predict(extgAge_seg_pov,predictDF_segpovext)
	forDerivext_RM=derivatives(extgAge_seg_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_segpovext,n=bemax)
	# rich female
	predictDF_segpovp$seg=rep("F",bpmax)
	forfitp_RF=predict(pgAge_seg_pov,predictDF_segpovp)
	forDerivp_RF=derivatives(pgAge_seg_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segpovp,n=bpmax)
	predictDF_segpovint$seg=rep("F",bimax)
	forfitint_RF=predict(intgAge_seg_pov,predictDF_segpovint)
	forDerivint_RF=derivatives(intgAge_seg_pov,term="s(cbcl_scr_syn_internal_r)",data=predictDF_segpovint,n=bimax)
	predictDF_segpovext$seg=rep("F",bemax)
	forfitext_RF=predict(extgAge_seg_pov,predictDF_segpovext)
	forDerivext_RF=derivatives(extgAge_seg_pov,term="s(cbcl_scr_syn_external_r)",data=predictDF_segpovext,n=bemax)
	# print out fit
	F_pFit[b,1:bpmax]=forfitp_F
	F_intFit[b,1:bimax]=forfitint_F
	F_extFit[b,1:bemax]=forfitext_F
	M_pFit[b,1:bpmax]=forfitp_M
	M_intFit[b,1:bimax]=forfitint_M
	M_extFit[b,1:bemax]=forfitext_M
	P_pFit[b,1:bpmax]=forfitp_P
	P_intFit[b,1:bimax]=forfitint_P
	P_extFit[b,1:bemax]=forfitext_P
	R_pFit[b,1:bpmax]=forfitp_R
	R_intFit[b,1:bimax]=forfitint_R
	R_extFit[b,1:bemax]=forfitext_R
	PF_pFit[b,1:bpmax]=forfitp_PF
	PF_intFit[b,1:bimax]=forfitint_PF
	PF_extFit[b,1:bemax]=forfitext_PF
	PM_pFit[b,1:bpmax]=forfitp_PM
	PM_intFit[b,1:bimax]=forfitint_PM
	PM_extFit[b,1:bemax]=forfitext_PM
	RM_pFit[b,1:bpmax]=forfitp_RM
	RM_intFit[b,1:bimax]=forfitint_RM
	RM_extFit[b,1:bemax]=forfitext_RM
	RF_pFit[b,1:bpmax]=forfitp_RF
	RF_intFit[b,1:bimax]=forfitint_RF
	RF_extFit[b,1:bemax]=forfitext_RF
	# print out derivatives
	F_pDeriv[b,1:bpmax]=forDerivp_F$derivative
	F_intDeriv[b,1:bimax]=forDerivint_F$derivative
	F_extDeriv[b,1:bemax]=forDerivext_F$derivative
	M_pDeriv[b,1:bpmax]=forDerivp_M$derivative
	M_intDeriv[b,1:bimax]=forDerivint_M$derivative
	M_extDeriv[b,1:bemax]=forDerivext_M$derivative
	P_pDeriv[b,1:bpmax]=forDerivp_P$derivative
	P_intDeriv[b,1:bimax]=forDerivint_P$derivative
	P_extDeriv[b,1:bemax]=forDerivext_P$derivative
	R_pDeriv[b,1:bpmax]=forDerivp_R$derivative
	R_intDeriv[b,1:bimax]=forDerivint_R$derivative
	R_extDeriv[b,1:bemax]=forDerivext_R$derivative
	PF_pDeriv[b,1:bpmax]=forDerivp_PF$derivative
	PF_intDeriv[b,1:bimax]=forDerivint_PF$derivative
	PF_extDeriv[b,1:bemax]=forDerivext_PF$derivative
	PM_pDeriv[b,1:bpmax]=forDerivp_PM$derivative
	PM_intDeriv[b,1:bimax]=forDerivint_PM$derivative
	PM_extDeriv[b,1:bemax]=forDerivext_PM$derivative
	RM_pDeriv[b,1:bpmax]=forDerivp_RM$derivative
	RM_intDeriv[b,1:bimax]=forDerivint_RM$derivative
	RM_extDeriv[b,1:bemax]=forDerivext_RM$derivative
	RF_pDeriv[b,1:bpmax]=forDerivp_RF$derivative
	RF_intDeriv[b,1:bimax]=forDerivint_RF$derivative
	RF_extDeriv[b,1:bemax]=forDerivext_RF$derivative
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
	intMax[b]=bimax
	extMax[b]=bemax
	
	#
	########## III NOW PRINT OUT F STATISTICS AND AIC, INCLUDING FOR NULLS (psuedopoverty and psuedoseg)
	#
	
	# model summaries
	reducedSum_p=summary(pgAge)
	segSum_p=summary(pgAge_seg)
	povSum_p=summary(pgAge_pov)
	segpovSum_p=summary(pgAge_seg_pov)
	reducedSum_int=summary(intgAge)
	segSum_int=summary(intgAge_seg)
	povSum_int=summary(intgAge_pov)
	segpovSum_int=summary(intgAge_seg_pov)
	reducedSum_ext=summary(extgAge)
	segSum_ext=summary(extgAge_seg)
	povSum_ext=summary(extgAge_pov)
	segpovSum_ext=summary(extgAge_seg_pov)
	# null model summaries
	seg_n_Sum_p=summary(pgAge_seg_n)
	pov_n_Sum_p=summary(pgAge_pov_n)
	segpov_n_Sum_p=summary(pgAge_seg_pov_n)
	seg_n_Sum_int=summary(intgAge_seg_n)
	pov_n_Sum_int=summary(intgAge_pov_n)
	segpov_n_Sum_int=summary(intgAge_seg_pov_n)
	seg_n_Sum_ext=summary(extgAge_seg_n)
	pov_n_Sum_ext=summary(extgAge_pov_n)
	segpov_n_Sum_ext=summary(extgAge_seg_pov_n)
	# extract smooths table
	reducedSmooths_p=reducedSum_p$s.table
	segSmooths_p=segSum_p$s.table
	povSmooths_p=povSum_p$s.table
	segpovSmooths_p=segpovSum_p$s.table
	reducedSmooths_int=reducedSum_int$s.table
	segSmooths_int=segSum_int$s.table
	povSmooths_int=povSum_int$s.table
	segpovSmooths_int=segpovSum_int$s.table
	reducedSmooths_ext=reducedSum_ext$s.table
	segSmooths_ext=segSum_ext$s.table
	povSmooths_ext=povSum_ext$s.table
	segpovSmooths_ext=segpovSum_ext$s.table
	# null smooths table
	seg_n_Smooths_p=seg_n_Sum_p$s.table
	pov_n_Smooths_p=pov_n_Sum_p$s.table
	segpov_n_Smooths_p=segpov_n_Sum_p$s.table
	seg_n_Smooths_int=seg_n_Sum_int$s.table
	pov_n_Smooths_int=pov_n_Sum_int$s.table
	segpov_n_Smooths_int=segpov_n_Sum_int$s.table
	seg_n_Smooths_ext=seg_n_Sum_ext$s.table
	pov_n_Smooths_ext=pov_n_Sum_ext$s.table
	segpov_n_Smooths_ext=segpov_n_Sum_ext$s.table
	### F stats
	# full
	Fseg_p[b]=segSmooths_p[1,'F']
	Fseg_int[b]=segSmooths_int[1,'F']
	Fseg_ext[b]=segSmooths_ext[1,'F']
	Fpov_p[b]=povSmooths_p[1,'F']
	Fpov_int[b]=povSmooths_int[1,'F']
	Fpov_ext[b]=povSmooths_ext[1,'F']
	# very full - note summing will likely render AIC comparison more useful
	Fsegpov_p[b]=sum(segpovSmooths_p[5:8,'F'])
	Fsegpov_int[b]=sum(segpovSmooths_int[5:8,'F'])
	Fsegpov_ext[b]=sum(segpovSmooths_ext[5:8,'F'])
	### AIC
	# reduced
	AICreduced_p[b]=AIC(pgAge)
	AICreduced_int[b]=AIC(intgAge)
	AICreduced_ext[b]=AIC(extgAge)
	# full
	AICseg_p[b]=AIC(pgAge_seg)
	AICseg_int[b]=AIC(intgAge_seg)
	AICseg_ext[b]=AIC(extgAge_seg)
	AICpov_p[b]=AIC(pgAge_pov)
	AICpov_int[b]=AIC(intgAge_pov)
	AICpov_ext[b]=AIC(extgAge_pov)
	# very full
	AICsegpov_p[b]=AIC(pgAge_seg_pov)
	AICsegpov_int[b]=AIC(intgAge_seg_pov)
	AICsegpov_ext[b]=AIC(extgAge_seg_pov)
	# AIC of direct comparison models: includes main effect for covariates but without interactions of interest
	AICseg_p_noIntrxn[b]=AIC(pgAge_seg_noIntrxn)
	AICseg_int_noIntrxn[b]=AIC(intgAge_seg_noIntrxn)
	AICseg_ext_noIntrxn[b]=AIC(extgAge_seg_noIntrxn)
	AICpov_p_noIntrxn[b]=AIC(pgAge_pov_noIntrxn)
	AICpov_int_noIntrxn[b]=AIC(intgAge_pov_noIntrxn)
	AICpov_ext_noIntrxn[b]=AIC(extgAge_pov_noIntrxn)
	AICsegpov_p_noIntrxn[b]=AIC(pgAge_seg_pov_noIntrxn)
	AICsegpov_int_noIntrxn[b]=AIC(intgAge_seg_pov_noIntrxn)
	AICsegpov_ext_noIntrxn[b]=AIC(extgAge_seg_pov_noIntrxn)
	### null F stats
	# full
	Fseg_null_p[b]=seg_n_Smooths_p[1,'F']
	Fseg_null_int[b]=seg_n_Smooths_int[1,'F']
	Fseg_null_ext[b]=seg_n_Smooths_ext[1,'F']
	Fpov_null_p[b]=pov_n_Smooths_p[1,'F']
	Fpov_null_int[b]=pov_n_Smooths_int[1,'F']
	Fpov_null_ext[b]=pov_n_Smooths_ext[1,'F']
	# very full
	Fsegpov_null_p[b]=sum(segpov_n_Smooths_p[5:8,'F'])
	Fsegpov_null_int[b]=sum(segpov_n_Smooths_int[5:8,'F'])
	Fsegpov_null_ext[b]=sum(segpov_n_Smooths_ext[5:8,'F'])
	### null AIC
	# full
	AICseg_null_p[b]=AIC(pgAge_seg_n)
	AICseg_null_int[b]=AIC(intgAge_seg_n)
	AICseg_null_ext[b]=AIC(extgAge_seg_n)
	AICpov_null_p[b]=AIC(pgAge_pov_n)
	AICpov_null_int[b]=AIC(intgAge_pov_n)
	AICpov_null_ext[b]=AIC(extgAge_pov_n)
	# very full
	AICsegpov_null_p[b]=AIC(pgAge_seg_pov_n)
	AICsegpov_null_int[b]=AIC(intgAge_seg_pov_n)
	AICsegpov_null_ext[b]=AIC(extgAge_seg_pov_n)

	# modular portion for poverty and equivalently-sized group from non-poverty:
	# extract poverty and pseudopov2 group
	povertyGroup=subset(bootSamp,poverty==1)
	pseudoPov2Group=subset(bootSamp,pseudopoverty2==1)
	# fit g~p to poverty group
	pov_gp=gam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=povertyGroup)
	# fit g~p to same-size nonpoverty group
	nonpov_gp=gam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=pseudoPov2Group)
	# make prediction df to use for both models
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax))
	colnames(predictDFp)=c('cbcl_scr_syn_totprob_r','interview_age')
	# extract fit (to become derivatives in postproc)
	povFit[b]=predict(pov_gp,newdata=predictDFp)
	pseudopovFit[b]=predict(nonpov_gp,newdata=predictDFp)
}
# SAVEOUT
# save out version with all F stats and AICs, include max values for all iterations as well
outdf=data.frame(Fseg_p,Fseg_int,Fseg_ext,Fpov_p,Fpov_int,Fpov_ext,Fsegpov_p,Fsegpov_int,Fsegpov_ext,AICreduced_p,AICreduced_int,AICreduced_ext,AICseg_p,AICseg_int,AICseg_ext,AICpov_p,AICpov_int,AICpov_ext,AICsegpov_p,AICsegpov_int,AICsegpov_ext,Fseg_null_p,Fseg_null_int,Fseg_null_ext,Fpov_null_p,Fpov_null_int,Fpov_null_ext,Fsegpov_null_p,Fsegpov_null_int,Fsegpov_null_ext,AICseg_null_p,AICseg_null_int,AICseg_null_ext,AICpov_null_p,AICpov_null_int,AICpov_null_ext,AICsegpov_null_p,AICsegpov_null_int,AICsegpov_null_ext,AICseg_p_noIntrxn,AICseg_int_noIntrxn,AICseg_ext_noIntrxn,AICpov_p_noIntrxn,AICpov_int_noIntrxn,AICpov_ext_noIntrxn,AICsegpov_p_noIntrxn,AICsegpov_int_noIntrxn,AICsegpov_ext_noIntrxn,pMax,intMax,extMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpFandAIC.rds')
# save out fits
outdf=data.frame(F_pFit,F_intFit,F_extFit,M_pFit,M_intFit,M_extFit,P_pFit,P_intFit,P_extFit,R_pFit,R_intFit,R_extFit,PF_pFit,PF_intFit,PF_extFit,PM_pFit,PM_intFit,PM_extFit,RM_pFit,RM_intFit,RM_extFit,RF_pFit,RF_intFit,RF_extFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpFits.rds')
# save out derivatives
outdf=data.frame(F_pDeriv,F_intDeriv,F_extDeriv,M_pDeriv,M_intDeriv,M_extDeriv,P_pDeriv,P_intDeriv,P_extDeriv,R_pDeriv,R_intDeriv,R_extDeriv,PF_pDeriv,PF_intDeriv,PF_extDeriv,PM_pDeriv,PM_intDeriv,PM_extDeriv,RM_pDeriv,RM_intDeriv,RM_extDeriv,RF_pDeriv,RF_intDeriv,RF_extDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpDerivs.rds')
print('done with g~p fit bootstrapping!')
# save out modular poverty and equivlanetly-sized poverty fits
outdf=data.frame(povFit,pseudopovFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpPovNonPov.rds')
