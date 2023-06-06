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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age')]
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
Fseg=rep(0,10000)
Fpoverty=rep(0,10000)
Fsegpoverty=rep(0,10000)
# AIC for each model
AICreduced=rep(0,10000)
AICseg=rep(0,10000)
AICpoverty=rep(0,10000)
AICsegpoverty=rep(0,10000)
# now equivalent vectors for null models
Fseg_null=rep(0,10000)
Fpoverty_null=rep(0,10000)
Fsegpoverty_null=rep(0,10000)
# AIC
AICseg_null=rep(0,10000)
AICpoverty_null=rep(0,10000)
AICsexpoverty_null=rep(0,10000)
# loop over manual bootstrap
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
	bootSamp$psuedoseg=rep(0,numSubjs)
	bootSamp$psuedoseg[sample(1:numSubjs,Fcount)]=1
	# randomly assign Povcount people to pseudopoverty
	bootSamp$psuedopoverty=rep(0,numSubjs)
	bootSamp$psuedopoverty[sample(1:numSubjs,Povcount)]=1
	######## I PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
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
	#### VERY FULL MODELS
	pgAge_seg_pov<-gam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=seg)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_totprob_r, by = interaction(seg, poverty)),data=masterdf)
	intgAge_seg_pov<-gam(g~s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=seg)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_internal_r, by = interaction(seg, poverty)),data=masterdf)
	extgAge_seg_pov<-gam(g~s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_external_r,by=seg)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_external_r, by = interaction(seg, poverty)),data=masterdf)
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
	# predict boy
	predictDF_segp$seg=rep("M",bpmax)
	forfitp_M=predict(pgAge_seg,predictDF_segp)
	predictDF_segint$seg=rep("M",bimax)
	forfitint_M=predict(intgAge_seg,predictDF_segint)
	predictDF_segext$seg=rep("M",bemax)
	forfitext_M=predict(extgAge_seg,predictDF_segext)
	# predict poverty
	forfitp_P=predict(pgAge_seg_pov,predictDF_povp)
	forfitint_P=predict(intgAge_seg_pov,predictDF_povint)
	forfitext_P=predict(extgAge_seg_pov,predictDF_povext)
	# predict rich
	predictDF_povp$poverty=rep("0",bpmax)
	forfitp_R=predict(pgAge_seg_pov,predictDF_povp)
	predictDF_povint$poverty=rep("0",bimax)
	forfitint_R=predict(intgAge_seg_pov,predictDF_povint)
	predictDF_povext$poverty=rep("0",bemax)
	fotfitext_R=predict(extgAge_seg_pov,predictDF_povext)
	## predict across triple interaction, funtime
	# poverty female
	forfitp_PF=predict(pgAge_seg_pov,predictDF_segpovp)
	forfitint_PF=predict(intgAge_seg_pov,predictDF_segpovint)
	forfitext_PF=predict(extgAge_seg_pov,predictDF_segpovext)
	# poverty male
	predictDF_segpovp$seg=rep("M",bpmax)
	forfitp_PM=predict(pgAge_seg_pov,predictDF_segpovp)
	predictDF_segpovint$seg=rep("M",bimax)
	forfitint_PM=predict(intgAge_seg_pov,predictDF_segpovint)
	predictDF_segpovext$seg=rep("M",bemax)
	forfitext_PM=predict(extgAge_seg_pov,predictDF_segpovext)
	# rich male
	predictDF_segpovp$poverty=rep("0",bpmax)
	forfitp_RM=predict(pgAge_seg_pov,predictDF_segpovp)
	predictDF_segpovint$poverty=rep("0",bimax)
	forfitint_RM=predict(intgAge_seg_pov,predictDF_segpovint)
	predictDF_segpovext$poverty=rep("0",bemax)
	forfitext_RM=predict(extgAge_seg_pov,predictDF_segpovext)
	# rich female
	predictDF_segpovp$seg=rep("F",bpmax)
	forfitp_RF=predict(pgAge_seg_pov,predictDF_segpovp)
	predictDF_segpovint$seg=rep("F",bimax)
	forfitint_RF=predict(intgAge_seg_pov,predictDF_segpovint)
	predictDF_segpovext$seg=rep("F",bemax)
	forfitext_RF=predict(extgAge_seg_pov,predictDF_segpovext)
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
	# use same procedure to print out deriatives of fit ADAPT DERIVATIVES TO FIT ON BOTH LEVELS OF FACTORS
	F_pDeriv[b,1:bpmax]=derivatives(pgAge_seg
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE,n=bpmax)
	forSplineint=derivatives(intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE,n=bimax)
	forSplineext=derivatives(extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE,n=bemax)
	forSplinesom=derivatives(somgAge,term='s(cbcl_scr_syn_somatic_r)',partial_match = TRUE,n=bsommax)
	forSplineanx=derivatives(anxgAge,term='s(cbcl_scr_syn_anxdep_r)',partial_match = TRUE,n=banxmax)
	forSplinetho=derivatives(thogAge,term='s(cbcl_scr_syn_thought_r)',partial_match = TRUE,n=bthomax)
	forSplinewit=derivatives(witgAge,term='s(cbcl_scr_syn_withdep_r)',partial_match = TRUE,n=bwitmax)
	forSplinesoc=derivatives(socgAge,term='s(cbcl_scr_syn_social_r)',partial_match = TRUE,n=bsocmax)
	forSplineatt=derivatives(attgAge,term='s(cbcl_scr_syn_attention_r)',partial_match = TRUE,n=battmax)
	forSplinerul=derivatives(rulgAge,term='s(cbcl_scr_syn_rulebreak_r)',partial_match = TRUE,n=brulmax)
	forSplineagg=derivatives(agggAge,term='s(cbcl_scr_syn_aggressive_r)',partial_match = TRUE,n=baggmax)
	# print out fit derivatives
	pDeriv[b,1:bpmax]=forSplinep$derivative
	intDeriv[b,1:bimax]=forSplineint$derivative
	extDeriv[b,1:bemax]=forSplineext$derivative
	somDeriv[b,1:bsommax]=forSplinesom$derivative
	anxDeriv[b,1:banxmax]=forSplineanx$derivative
	thoDeriv[b,1:bthomax]=forSplinetho$derivative
	witDeriv[b,1:bwitmax]=forSplinewit$derivative
	socDeriv[b,1:bsocmax]=forSplinesoc$derivative
	attDeriv[b,1:battmax]=forSplineatt$derivative
	rulDeriv[b,1:brulmax]=forSplinerul$derivative
	aggDeriv[b,1:baggmax]=forSplineagg$derivative
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
	intMax[b]=bimax
	extMax[b]=bemax
	somMax[b]=bsommax
	anxMax[b]=banxmax
	thoMax[b]=bthomax
	witMax[b]=bwitmax
	socMax[b]=bsocmax
	attMax[b]=battmax
	rulMax[b]=brulmax
	aggMax[b]=baggmax
	########## II NOW PRINT OUT F STATISTICS AND AIC, INCLUDING FOR NULLS
}
# SAVEOUT
# save out version with all cbcl factors
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,somlinBoots,anxlinBoots,tholinBoots,witlinBoots,soclinBoots,attlinBoots,rullinBoots,agglinBoots,pMax,intMax,extMax,somMax,anxMax,thoMax,witMax,socMax,attMax,rulMax,aggMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,socDeriv,attDeriv,rulDeriv,aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots.rds')
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots.rds')

print('done with g~p fit bootstrapping!')
