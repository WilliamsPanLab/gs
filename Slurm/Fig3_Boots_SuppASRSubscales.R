library(mgcv)
library(pammtools)
library(patchwork)
library(gratia)
library(scales)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','ASRAnxDepr','ASRWithdrawn','ASRSomatic','ASRThought','ASRAttn','ASRAggr','ASRRulB','ASRInt','ASRExt','g','subjectkey','interview_age')]
masterdf$ASR_anxdep=as.numeric(masterdf$ASRAnxDepr)
masterdf$ASR_withdep=as.numeric(masterdf$ASRWithdrawn)
masterdf$ASR_somatic=as.numeric(masterdf$ASRSomatic)
masterdf$ASR_thought=as.numeric(masterdf$ASRThought)
masterdf$ASR_attention=as.numeric(masterdf$ASRAttn)
masterdf$ASR_aggressive=as.numeric(masterdf$ASRAggr)
masterdf$ASR_rulebreak=as.numeric(masterdf$ASRRulB)
masterdf$ASRInt=as.numeric(masterdf$ASRInt)
masterdf$ASRExt=as.numeric(masterdf$ASRExt)
# get length of df for later
lenDF=dim(masterdf)[1]
# will need to get full and reduced models for each boot, as well as a null distribution
# in addition to derivatives and fits, save F values of interaction + null distribution F values
# interactions to be tested: p*sex, p*poverty, p*sex*poverty
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r=as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r=as.numeric(masterdf$cbcl_scr_syn_external_r)
masterdf$cbcl_scr_syn_somatic_r=as.numeric(masterdf$cbcl_scr_syn_somatic_r)
masterdf$cbcl_scr_syn_anxdep_r=as.numeric(masterdf$cbcl_scr_syn_anxdep_r)
masterdf$cbcl_scr_syn_thought_r=as.numeric(masterdf$cbcl_scr_syn_thought_r)
masterdf$cbcl_scr_syn_withdep_r=as.numeric(masterdf$cbcl_scr_syn_withdep_r)
masterdf$cbcl_scr_syn_social_r=as.numeric(masterdf$cbcl_scr_syn_social_r)
masterdf$cbcl_scr_syn_attention_r=as.numeric(masterdf$cbcl_scr_syn_attention_r)
masterdf$cbcl_scr_syn_rulebreak_r=as.numeric(masterdf$cbcl_scr_syn_rulebreak_r)
masterdf$cbcl_scr_syn_aggressive_r=as.numeric(masterdf$cbcl_scr_syn_aggressive_r)
# sex (g in for x to get co-pilot around its NSFW filter) and poverty to factors
masterdf$seg<-as.ordered(masterdf$sex)
masterdf$poverty=0
masterdf$income<-as.numeric(masterdf$income)
# note that poverty is defined as income < 5: https://collection3165.readthedocs.io/en/stable/recommendations/#2-the-bids-participants-files-and-matched-groups
masterdf$poverty[masterdf$income<5]=1
masterdf$poverty=as.ordered(masterdf$poverty)
### initialize cross-boot vectors
# predicted derivatives: set to maximum value for ncol +1, as 0-maxvalue is 1 longer than maxvalue
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
F_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
M_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
P_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
R_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
# predicted values: set to maximum value for ncol
F_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
M_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
P_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
R_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
# make matrices for all the subscale fits as well (sep. matrices for girl, boy, nonpoverty, poverty)
intMaxVal=max(masterdf$cbcl_scr_syn_internal_r)
extMaxVal=max(masterdf$cbcl_scr_syn_external_r)
somMaxVal=max(masterdf$cbcl_scr_syn_somatic_r)
anxMaxVal=max(masterdf$cbcl_scr_syn_anxdep_r)
thoMaxVal=max(masterdf$cbcl_scr_syn_thought_r)
depMaxVal=max(masterdf$cbcl_scr_syn_withdep_r)
socMaxVal=max(masterdf$cbcl_scr_syn_social_r)
attMaxVal=max(masterdf$cbcl_scr_syn_attention_r)
rulMaxVal=max(masterdf$cbcl_scr_syn_rulebreak_r)
aggMaxVal=max(masterdf$cbcl_scr_syn_aggressive_r)
# internalizing
F_intFit=matrix(0,nrow=10000,ncol=(intMaxVal+1))
M_intFit=matrix(0,nrow=10000,ncol=(intMaxVal+1))
P_intFit=matrix(0,nrow=10000,ncol=(intMaxVal+1))
R_intFit=matrix(0,nrow=10000,ncol=(intMaxVal+1))
# externalizing
F_extFit=matrix(0,nrow=10000,ncol=(extMaxVal+1))
M_extFit=matrix(0,nrow=10000,ncol=(extMaxVal+1))
P_extFit=matrix(0,nrow=10000,ncol=(extMaxVal+1))
R_extFit=matrix(0,nrow=10000,ncol=(extMaxVal+1))
# somatic
F_somFit=matrix(0,nrow=10000,ncol=(somMaxVal+1))
M_somFit=matrix(0,nrow=10000,ncol=(somMaxVal+1))
P_somFit=matrix(0,nrow=10000,ncol=(somMaxVal+1))
R_somFit=matrix(0,nrow=10000,ncol=(somMaxVal+1))
# anxiety/depression
F_anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal+1))
M_anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal+1))
P_anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal+1))
R_anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal+1))
# thought
F_thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal+1))
M_thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal+1))
P_thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal+1))
R_thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal+1))
# withdrawn depression
F_witdepFit=matrix(0,nrow=10000,ncol=(depMaxVal+1))
M_witdepFit=matrix(0,nrow=10000,ncol=(depMaxVal+1))
P_witdepFit=matrix(0,nrow=10000,ncol=(depMaxVal+1))
R_witdepFit=matrix(0,nrow=10000,ncol=(depMaxVal+1))
# social
F_socFit=matrix(0,nrow=10000,ncol=(socMaxVal+1))
M_socFit=matrix(0,nrow=10000,ncol=(socMaxVal+1))
P_socFit=matrix(0,nrow=10000,ncol=(socMaxVal+1))
R_socFit=matrix(0,nrow=10000,ncol=(socMaxVal+1))
# attention
F_attFit=matrix(0,nrow=10000,ncol=(attMaxVal+1))
M_attFit=matrix(0,nrow=10000,ncol=(attMaxVal+1))
P_attFit=matrix(0,nrow=10000,ncol=(attMaxVal+1))
R_attFit=matrix(0,nrow=10000,ncol=(attMaxVal+1))
# rule-breaking
F_rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal+1))
M_rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal+1))
P_rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal+1))
R_rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal+1))
# aggression
F_aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal+1))
M_aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal+1))
P_aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal+1))
R_aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal+1))

# a lil' vector just to track max p int and ext over iterations
pMax=rep(0,10000)
# and modular fit for pov vs. psuedopov 2 (matched #) fits
povFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
pseudopovFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
FullNonpovFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
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
	bintmax=max(bootSamp$cbcl_scr_syn_internal_r)
	bextmax=max(bootSamp$cbcl_scr_syn_external_r)
	bsommax=max(bootSamp$cbcl_scr_syn_somatic_r)
	banxmax=max(bootSamp$cbcl_scr_syn_anxdep_r)
	bthomax=max(bootSamp$cbcl_scr_syn_thought_r)
	bwitdepmax=max(bootSamp$cbcl_scr_syn_withdep_r)
	bsocmax=max(bootSamp$cbcl_scr_syn_social_r)
	battmax=max(bootSamp$cbcl_scr_syn_attention_r)
	brulmax=max(bootSamp$cbcl_scr_syn_rulebreak_r)
	baggmax=max(bootSamp$cbcl_scr_syn_aggressive_r)
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
	#### FULL MODELS
	pgAge_seg<-bam(g~s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg+s(interview_age),data=bootSamp)
	pgAge_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r)+poverty+s(interview_age),data=bootSamp)
	IntgAge_seg<-bam(g~s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_internal_r)+seg+s(interview_age),data=bootSamp)
	ExtgAge_seg<-bam(g~s(cbcl_scr_syn_external_r,by=seg)+s(cbcl_scr_syn_external_r)+seg+s(interview_age),data=bootSamp)
	SomgAge_seg<-bam(g~s(cbcl_scr_syn_somatic_r,by=seg)+s(cbcl_scr_syn_somatic_r)+seg+s(interview_age),data=bootSamp)
	AnxgAge_seg<-bam(g~s(cbcl_scr_syn_anxdep_r,by=seg)+s(cbcl_scr_syn_anxdep_r)+seg+s(interview_age),data=bootSamp)
	ThogAge_seg<-bam(g~s(cbcl_scr_syn_thought_r,by=seg)+s(cbcl_scr_syn_thought_r)+seg+s(interview_age),data=bootSamp)
	WitDepgAge_seg<-bam(g~s(cbcl_scr_syn_withdep_r,by=seg)+s(cbcl_scr_syn_withdep_r)+seg+s(interview_age),data=bootSamp)
	SocgAge_seg<-bam(g~s(cbcl_scr_syn_social_r,by=seg)+s(cbcl_scr_syn_social_r)+seg+s(interview_age),data=bootSamp)
	AttgAge_seg<-bam(g~s(cbcl_scr_syn_attention_r,by=seg)+s(cbcl_scr_syn_attention_r)+seg+s(interview_age),data=bootSamp)
	RulgAge_seg<-bam(g~s(cbcl_scr_syn_rulebreak_r,by=seg)+s(cbcl_scr_syn_rulebreak_r)+seg+s(interview_age),data=bootSamp)
	AgggAge_seg<-bam(g~s(cbcl_scr_syn_aggressive_r,by=seg)+s(cbcl_scr_syn_aggressive_r)+seg+s(interview_age),data=bootSamp)
	IntgAge_pov<-bam(g~s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_internal_r)+poverty+s(interview_age),data=bootSamp)
	ExtgAge_pov<-bam(g~s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_external_r)+poverty+s(interview_age),data=bootSamp)
	SomgAge_pov<-bam(g~s(cbcl_scr_syn_somatic_r,by=poverty)+s(cbcl_scr_syn_somatic_r)+poverty+s(interview_age),data=bootSamp)
	AnxgAge_pov<-bam(g~s(cbcl_scr_syn_anxdep_r,by=poverty)+s(cbcl_scr_syn_anxdep_r)+poverty+s(interview_age),data=bootSamp)
	ThogAge_pov<-bam(g~s(cbcl_scr_syn_thought_r,by=poverty)+s(cbcl_scr_syn_thought_r)+poverty+s(interview_age),data=bootSamp)
	WitDepgAge_pov<-bam(g~s(cbcl_scr_syn_withdep_r,by=poverty)+s(cbcl_scr_syn_withdep_r)+poverty+s(interview_age),data=bootSamp)
	SocgAge_pov<-bam(g~s(cbcl_scr_syn_social_r,by=poverty)+s(cbcl_scr_syn_social_r)+poverty+s(interview_age),data=bootSamp)
	AttgAge_pov<-bam(g~s(cbcl_scr_syn_attention_r,by=poverty)+s(cbcl_scr_syn_attention_r)+poverty+s(interview_age),data=bootSamp)
	RulgAge_pov<-bam(g~s(cbcl_scr_syn_rulebreak_r,by=poverty)+s(cbcl_scr_syn_rulebreak_r)+poverty+s(interview_age),data=bootSamp)
	AgggAge_pov<-bam(g~s(cbcl_scr_syn_aggressive_r,by=poverty)+s(cbcl_scr_syn_aggressive_r)+poverty+s(interview_age),data=bootSamp)

	# comparison models to evaluate AIC gain from interactions specifically. Should include main effects but minus interaction of interest
	pgAge_seg_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r)+seg+s(interview_age),data=bootSamp)
	pgAge_pov_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r)+poverty+s(interview_age),data=bootSamp)
	pgAge_seg_pov_noIntrxn<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg*poverty+s(interview_age),data=bootSamp)
	#### VERY FULL MODELS
	pgAge_seg_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_totprob_r)+seg*poverty+s(interview_age)+s(cbcl_scr_syn_totprob_r, by = interaction(seg, poverty)),data=bootSamp)
	# fit null models for use later (no reduced, reduced is the same as real reduced)
	#### FULL NULL MODELS
	pgAge_seg_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedoseg)+s(cbcl_scr_syn_totprob_r)+psuedoseg+s(interview_age),data=bootSamp)
	pgAge_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty)+s(cbcl_scr_syn_totprob_r)+psuedopoverty+s(interview_age),data=bootSamp)
	#### VERY FULL NULL MODELS
	pgAge_seg_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty)+s(cbcl_scr_syn_totprob_r,by=psuedoseg)+s(cbcl_scr_syn_totprob_r)+psuedoseg*psuedopoverty+s(interview_age)+s(cbcl_scr_syn_totprob_r, by = interaction(psuedoseg, psuedopoverty)),data=bootSamp)
	
	#
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#
	
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(0:bpmax)
	eachIntcount=seq(0:bintmax)
	eachExtcount=seq(0:bextmax)
	eachSomcount=seq(0:bsommax)
	eachAnxcount=seq(0:banxmax)
	eachThocount=seq(0:bthomax)
	eachWitDepcount=seq(0:bwitdepmax)
	eachSoccount=seq(0:bsocmax)
	eachAttcount=seq(0:battmax)
	eachRulcount=seq(0:brulmax)
	eachAggcount=seq(0:baggmax)

	# set age to to median for predict df, +1 because 0:max sequence is one lnger than max
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)))
	predictDFi=data.frame(eachIntcount,rep(median(bootSamp$interview_age),(bintmax+1)))
	predictDFe=data.frame(eachExtcount,rep(median(bootSamp$interview_age),(bextmax+1)))
	predictDFs=data.frame(eachSomcount,rep(median(bootSamp$interview_age),(bsommax+1)))
	predictDFa=data.frame(eachAnxcount,rep(median(bootSamp$interview_age),(banxmax+1)))
	predictDFt=data.frame(eachThocount,rep(median(bootSamp$interview_age),(bthomax+1)))
	predictDFwd=data.frame(eachWitDepcount,rep(median(bootSamp$interview_age),(bwitdepmax+1)))
	predictDFsoc=data.frame(eachSoccount,rep(median(bootSamp$interview_age),(bsocmax+1)))
	predictDFatt=data.frame(eachAttcount,rep(median(bootSamp$interview_age),(battmax+1)))
	predictDFrul=data.frame(eachRulcount,rep(median(bootSamp$interview_age),(brulmax+1)))
	predictDFagg=data.frame(eachAggcount,rep(median(bootSamp$interview_age),(baggmax+1)))
	# set predict df to interacting factors of interest (seg, poverty)
	predictDF_segp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("F",(bpmax+1)))
	predictDF_povp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("1",(bpmax+1)))
	predictDF_segpovp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("F",(bpmax+1)),rep("1",(bpmax+1)))
	predictDF_segInt=data.frame(eachIntcount,rep(median(bootSamp$interview_age),(bintmax+1)),rep("F",(bintmax+1)))
	predictDF_povInt=data.frame(eachIntcount,rep(median(bootSamp$interview_age),(bintmax+1)),rep("1",(bintmax+1)))
	predictDF_segExt=data.frame(eachExtcount,rep(median(bootSamp$interview_age),(bextmax+1)),rep("F",(bextmax+1)))
	predictDF_povExt=data.frame(eachExtcount,rep(median(bootSamp$interview_age),(bextmax+1)),rep("1",(bextmax+1)))
	predictDF_segSom=data.frame(eachSomcount,rep(median(bootSamp$interview_age),(bsommax+1)),rep("F",(bsommax+1)))
	predictDF_povSom=data.frame(eachSomcount,rep(median(bootSamp$interview_age),(bsommax+1)),rep("1",(bsommax+1)))
	predictDF_segAnx=data.frame(eachAnxcount,rep(median(bootSamp$interview_age),(banxmax+1)),rep("F",(banxmax+1)))
	predictDF_povAnx=data.frame(eachAnxcount,rep(median(bootSamp$interview_age),(banxmax+1)),rep("1",(banxmax+1)))
	predictDF_segTho=data.frame(eachThocount,rep(median(bootSamp$interview_age),(bthomax+1)),rep("F",(bthomax+1)))
	predictDF_povTho=data.frame(eachThocount,rep(median(bootSamp$interview_age),(bthomax+1)),rep("1",(bthomax+1)))
	predictDF_segWitDep=data.frame(eachWitDepcount,rep(median(bootSamp$interview_age),(bwitdepmax+1)),rep("F",(bwitdepmax+1)))
	predictDF_povWitDep=data.frame(eachWitDepcount,rep(median(bootSamp$interview_age),(bwitdepmax+1)),rep("1",(bwitdepmax+1)))
	predictDF_segSoc=data.frame(eachSoccount,rep(median(bootSamp$interview_age),(bsocmax+1)),rep("F",(bsocmax+1)))
	predictDF_povSoc=data.frame(eachSoccount,rep(median(bootSamp$interview_age),(bsocmax+1)),rep("1",(bsocmax+1)))
	predictDF_segAtt=data.frame(eachAttcount,rep(median(bootSamp$interview_age),(battmax+1)),rep("F",(battmax+1)))
	predictDF_povAtt=data.frame(eachAttcount,rep(median(bootSamp$interview_age),(battmax+1)),rep("1",(battmax+1)))
	predictDF_segRul=data.frame(eachRulcount,rep(median(bootSamp$interview_age),(brulmax+1)),rep("F",(brulmax+1)))
	predictDF_povRul=data.frame(eachRulcount,rep(median(bootSamp$interview_age),(brulmax+1)),rep("1",(brulmax+1)))
	predictDF_segAgg=data.frame(eachAggcount,rep(median(bootSamp$interview_age),(baggmax+1)),rep("F",(baggmax+1)))
	predictDF_povAgg=data.frame(eachAggcount,rep(median(bootSamp$interview_age),(baggmax+1)),rep("1",(baggmax+1)))
	# set colnames so predict can work
	colnames(predictDF_segp)=c('cbcl_scr_syn_totprob_r','interview_age','seg')
	colnames(predictDF_povp)=c('cbcl_scr_syn_totprob_r','interview_age','poverty')
	colnames(predictDF_segpovp)=c('cbcl_scr_syn_totprob_r','interview_age','seg','poverty')
	colnames(predictDF_segInt)=c('cbcl_scr_syn_internal_r','interview_age','seg')
	colnames(predictDF_povInt)=c('cbcl_scr_syn_internal_r','interview_age','poverty')
	colnames(predictDF_segExt)=c('cbcl_scr_syn_external_r','interview_age','seg')
	colnames(predictDF_povExt)=c('cbcl_scr_syn_external_r','interview_age','poverty')
	colnames(predictDF_segSom)=c('cbcl_scr_syn_somatic_r','interview_age','seg')
	colnames(predictDF_povSom)=c('cbcl_scr_syn_somatic_r','interview_age','poverty')
	colnames(predictDF_segAnx)=c('cbcl_scr_syn_anxdep_r','interview_age','seg')
	colnames(predictDF_povAnx)=c('cbcl_scr_syn_anxdep_r','interview_age','poverty')
	colnames(predictDF_segTho)=c('cbcl_scr_syn_thought_r','interview_age','seg')
	colnames(predictDF_povTho)=c('cbcl_scr_syn_thought_r','interview_age','poverty')
	colnames(predictDF_segWitDep)=c('cbcl_scr_syn_withdep_r','interview_age','seg')
	colnames(predictDF_povWitDep)=c('cbcl_scr_syn_withdep_r','interview_age','poverty')
	colnames(predictDF_segSoc)=c('cbcl_scr_syn_social_r','interview_age','seg')
	colnames(predictDF_povSoc)=c('cbcl_scr_syn_social_r','interview_age','poverty')
	colnames(predictDF_segAtt)=c('cbcl_scr_syn_attention_r','interview_age','seg')
	colnames(predictDF_povAtt)=c('cbcl_scr_syn_attention_r','interview_age','poverty')
	colnames(predictDF_segRul)=c('cbcl_scr_syn_rulebreak_r','interview_age','seg')
	colnames(predictDF_povRul)=c('cbcl_scr_syn_rulebreak_r','interview_age','poverty')
	colnames(predictDF_segAgg)=c('cbcl_scr_syn_aggressive_r','interview_age','seg')
	colnames(predictDF_povAgg)=c('cbcl_scr_syn_aggressive_r','interview_age','poverty')

	# predict girl
	forfitp_F=predict(pgAge_seg,predictDF_segp)
	forDerivp_F=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=(bpmax+1))
	forfitint_F=predict(IntgAge_seg,predictDF_segInt)
	forfitext_F=predict(ExtgAge_seg,predictDF_segExt)
	forfitsom_F=predict(SomgAge_seg,predictDF_segSom)
	forfitanx_F=predict(AnxgAge_seg,predictDF_segAnx)
	forfittho_F=predict(ThogAge_seg,predictDF_segTho)
	forfitwitdep_F=predict(WitDepgAge_seg,predictDF_segWitDep)
	forfitsoc_F=predict(SocgAge_seg,predictDF_segSoc)
	forfitatt_F=predict(AttgAge_seg,predictDF_segAtt)
	forfitrul_F=predict(RulgAge_seg,predictDF_segRul)
	forfitagg_F=predict(AgggAge_seg,predictDF_segAgg)
	# predict boy
	predictDF_segp$seg=rep("M",(bpmax+1))
	predictDF_segInt$seg=rep("M",(bintmax+1))
	predictDF_segExt$seg=rep("M",(bextmax+1))
	predictDF_segSom$seg=rep("M",(bsommax+1))
	predictDF_segAnx$seg=rep("M",(banxmax+1))
	predictDF_segTho$seg=rep("M",(bthomax+1))
	predictDF_segWitDep$seg=rep("M",(bwitdepmax+1))
	predictDF_segSoc$seg=rep("M",(bsocmax+1))
	predictDF_segAtt$seg=rep("M",(battmax+1))
	predictDF_segRul$seg=rep("M",(brulmax+1))
	predictDF_segAgg$seg=rep("M",(baggmax+1))
	forfitp_M=predict(pgAge_seg,predictDF_segp)
	forfitint_M=predict(IntgAge_seg,predictDF_segInt)
	forfitext_M=predict(ExtgAge_seg,predictDF_segExt)
	forfitsom_M=predict(SomgAge_seg,predictDF_segSom)
	forfitanx_M=predict(AnxgAge_seg,predictDF_segAnx)
	forfittho_M=predict(ThogAge_seg,predictDF_segTho)
	forfitwitdep_M=predict(WitDepgAge_seg,predictDF_segWitDep)
	forfitsoc_M=predict(SocgAge_seg,predictDF_segSoc)
	forfitatt_M=predict(AttgAge_seg,predictDF_segAtt)
	forfitrul_M=predict(RulgAge_seg,predictDF_segRul)
	forfitagg_M=predict(AgggAge_seg,predictDF_segAgg)
	forDerivp_M=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=(bpmax+1))
	# predict poverty
	forfitp_P=predict(pgAge_pov,predictDF_povp)
	forDerivp_P=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=(bpmax+1))
	forfitint_P=predict(IntgAge_pov,predictDF_povInt)
	forfitext_P=predict(ExtgAge_pov,predictDF_povExt)
	forfitsom_P=predict(SomgAge_pov,predictDF_povSom)
	forfitanx_P=predict(AnxgAge_pov,predictDF_povAnx)
	forfittho_P=predict(ThogAge_pov,predictDF_povTho)
	forfitwitdep_P=predict(WitDepgAge_pov,predictDF_povWitDep)
	forfitsoc_P=predict(SocgAge_pov,predictDF_povSoc)
	forfitatt_P=predict(AttgAge_pov,predictDF_povAtt)
	forfitrul_P=predict(RulgAge_pov,predictDF_povRul)
	forfitagg_P=predict(AgggAge_pov,predictDF_povAgg)
	# predict rich
	predictDF_povp$poverty=rep("0",(bpmax+1))
	predictDF_povInt$poverty=rep("0",(bintmax+1))
	predictDF_povExt$poverty=rep("0",(bextmax+1))
	predictDF_povSom$poverty=rep("0",(bsommax+1))
	predictDF_povAnx$poverty=rep("0",(banxmax+1))
	predictDF_povTho$poverty=rep("0",(bthomax+1))
	predictDF_povWitDep$poverty=rep("0",(bwitdepmax+1))
	predictDF_povSoc$poverty=rep("0",(bsocmax+1))
	predictDF_povAtt$poverty=rep("0",(battmax+1))
	predictDF_povRul$poverty=rep("0",(brulmax+1))
	predictDF_povAgg$poverty=rep("0",(baggmax+1))
	forfitp_R=predict(pgAge_pov,predictDF_povp)
	forfitint_R=predict(IntgAge_pov,predictDF_povInt)
	forfitext_R=predict(ExtgAge_pov,predictDF_povExt)
	forfitsom_R=predict(SomgAge_pov,predictDF_povSom)
	forfitanx_R=predict(AnxgAge_pov,predictDF_povAnx)
	forfittho_R=predict(ThogAge_pov,predictDF_povTho)
	forfitwitdep_R=predict(WitDepgAge_pov,predictDF_povWitDep)
	forfitsoc_R=predict(SocgAge_pov,predictDF_povSoc)
	forfitatt_R=predict(AttgAge_pov,predictDF_povAtt)
	forfitrul_R=predict(RulgAge_pov,predictDF_povRul)
	forfitagg_R=predict(AgggAge_pov,predictDF_povAgg)
	forDerivp_R=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=(bpmax+1))
	# print out fit
	F_pFit[b,1:(bpmax+1)]=forfitp_F
	M_pFit[b,1:(bpmax+1)]=forfitp_M
	P_pFit[b,1:(bpmax+1)]=forfitp_P
	R_pFit[b,1:(bpmax+1)]=forfitp_R
	F_intFit[b,1:(bintmax+1)]=forfitint_F
	M_intFit[b,1:(bintmax+1)]=forfitint_M
	P_intFit[b,1:(bintmax+1)]=forfitint_P
	R_intFit[b,1:(bintmax+1)]=forfitint_R
	F_extFit[b,1:(bextmax+1)]=forfitext_F
	M_extFit[b,1:(bextmax+1)]=forfitext_M
	P_extFit[b,1:(bextmax+1)]=forfitext_P
	R_extFit[b,1:(bextmax+1)]=forfitext_R
	F_somFit[b,1:(bsommax+1)]=forfitsom_F
	M_somFit[b,1:(bsommax+1)]=forfitsom_M
	P_somFit[b,1:(bsommax+1)]=forfitsom_P
	R_somFit[b,1:(bsommax+1)]=forfitsom_R
	F_anxFit[b,1:(banxmax+1)]=forfitanx_F
	M_anxFit[b,1:(banxmax+1)]=forfitanx_M
	P_anxFit[b,1:(banxmax+1)]=forfitanx_P
	R_anxFit[b,1:(banxmax+1)]=forfitanx_R
	F_thoFit[b,1:(bthomax+1)]=forfittho_F
	M_thoFit[b,1:(bthomax+1)]=forfittho_M
	P_thoFit[b,1:(bthomax+1)]=forfittho_P
	R_thoFit[b,1:(bthomax+1)]=forfittho_R
	F_witdepFit[b,1:(bwitdepmax+1)]=forfitwitdep_F
	M_witdepFit[b,1:(bwitdepmax+1)]=forfitwitdep_M
	P_witdepFit[b,1:(bwitdepmax+1)]=forfitwitdep_P
	R_witdepFit[b,1:(bwitdepmax+1)]=forfitwitdep_R
	F_socFit[b,1:(bsocmax+1)]=forfitsoc_F
	M_socFit[b,1:(bsocmax+1)]=forfitsoc_M
	P_socFit[b,1:(bsocmax+1)]=forfitsoc_P
	R_socFit[b,1:(bsocmax+1)]=forfitsoc_R
	F_attFit[b,1:(battmax+1)]=forfitatt_F
	M_attFit[b,1:(battmax+1)]=forfitatt_M
	P_attFit[b,1:(battmax+1)]=forfitatt_P
	R_attFit[b,1:(battmax+1)]=forfitatt_R
	F_rulFit[b,1:(brulmax+1)]=forfitrul_F
	M_rulFit[b,1:(brulmax+1)]=forfitrul_M
	P_rulFit[b,1:(brulmax+1)]=forfitrul_P
	R_rulFit[b,1:(brulmax+1)]=forfitrul_R
	F_aggFit[b,1:(baggmax+1)]=forfitagg_F
	M_aggFit[b,1:(baggmax+1)]=forfitagg_M
	P_aggFit[b,1:(baggmax+1)]=forfitagg_P
	R_aggFit[b,1:(baggmax+1)]=forfitagg_R
	# print out derivatives
	F_pDeriv[b,1:(bpmax+1)]=forDerivp_F$derivative
	M_pDeriv[b,1:(bpmax+1)]=forDerivp_M$derivative
	P_pDeriv[b,1:(bpmax+1)]=forDerivp_P$derivative
	R_pDeriv[b,1:(bpmax+1)]=forDerivp_R$derivative
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
	#
	# III modular portion for poverty and equivalently-sized group from non-poverty:
	# extract poverty and pseudopov2 group
	povertyGroup=subset(bootSamp,poverty==1)
	pseudoPov2Group=subset(bootSamp,pseudopoverty2==1)
	nonpovGroup=subset(bootSamp,poverty==0)
	# find maximum cbcl_scr_syn_totprob_r value present across both subgroups, make prediction df with that range of cbcl_scr_syn_totprob_r values
	bpmax_2=min(c(max(povertyGroup$cbcl_scr_syn_totprob_r),max(pseudoPov2Group$cbcl_scr_syn_totprob_r)))
	eachPcount2=seq(0,bpmax_2)
	predictDF_2=data.frame(eachPcount2,rep(median(bootSamp$interview_age),(bpmax_2+1)))
	# fit g~p to poverty group
	pov_gp=gam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=povertyGroup)
	# fit g~p to same-size nonpoverty group
	nonpov_gp=gam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=pseudoPov2Group)
	# fit g~p to the full nonpoverty group
	fnonpov_gp=gam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=nonpovGroup)
	# make prediction df to use for both models
	colnames(predictDF_2)=c('cbcl_scr_syn_totprob_r','interview_age')
	# extract fit derivatives
	povFit[b,1:(bpmax_2+1)]=derivatives(pov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative
	pseudopovFit[b,1:(bpmax_2+1)]=derivatives(nonpov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative
	FullNonpovFit[b,1:(bpmax_2+1)]=derivatives(fnonpov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative

}
# SAVEOUT
# save out version with all F stats and AICs, include max values for all iterations as well
# save out fits
outdf=data.frame(F_pFit,M_pFit,P_pFit,R_pFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpFits.rds')
# save out fits of all 10 subscales (int, ext, som, anx, tho, withdep, soc, rulbreak, attention, agg) and all 4 (F,M,P,R) versions
outdf=data.frame(F_intFit,M_intFit,P_intFit,R_intFit,F_extFit,M_extFit,P_extFit,R_extFit,F_somFit,M_somFit,P_somFit,R_somFit,F_anxFit,M_anxFit,P_anxFit,R_anxFit,F_thoFit,M_thoFit,P_thoFit,R_thoFit,F_witdepFit,M_witdepFit,P_witdepFit,R_witdepFit,F_socFit,M_socFit,P_socFit,R_socFit,F_rulFit,M_rulFit,P_rulFit,R_rulFit,F_attFit,M_attFit,P_attFit,R_attFit,F_aggFit,M_aggFit,P_aggFit,R_aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpSubscaleFits.rds')
# save out derivatives
outdf=data.frame(F_pDeriv,M_pDeriv,P_pDeriv,R_pDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpDerivs.rds')
print('done with g~p fit bootstrapping!')
# save out modular poverty and equivlanetly-sized poverty fits
outdf=data.frame(povFit,pseudopovFit,FullNonpovFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpPovNonPov.rds')
