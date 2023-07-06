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
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','ASRAnxDepr','ASRWithdrawn','ASRSomatic','ASRThought','ASRAttn','ASRAggr','ASRRulB','ASRInt','ASRExt','g','subjectkey','interview_age')]
# get length of df for later
lenDF=dim(masterdf)[1]
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
masterdf$ASR_anxdep=as.numeric(masterdf$ASRAnxDepr)
masterdf$ASR_withdep=as.numeric(masterdf$ASRWithdrawn)
masterdf$ASR_somatic=as.numeric(masterdf$ASRSomatic)
masterdf$ASR_thought=as.numeric(masterdf$ASRThought)
masterdf$ASR_attention=as.numeric(masterdf$ASRAttn)
masterdf$ASR_aggressive=as.numeric(masterdf$ASRAggr)
masterdf$ASR_rulebreak=as.numeric(masterdf$ASRRulB)
masterdf$ASRInt=as.numeric(masterdf$ASRInt)
masterdf$ASRExt=as.numeric(masterdf$ASRExt)
# Note social for children is removed due to lack of equivalent and intrusive for adults
### initialize cross-boot vectors
# linear? 1xbootlength
plinBoots=rep(0,10000)
intlinBoots=rep(0,10000)
extlinBoots=rep(0,10000)
somLinBoots=rep(0,10000)
anxLinBoots=rep(0,10000)
thoLinBoots=rep(0,10000)
witLinBoots=rep(0,10000)
attLinBoots=rep(0,10000)
rulLinBoots=rep(0,10000)
aggLinBoots=rep(0,10000)
asrpLinBoots=rep(0,10000)
asrintLinBoots=rep(0,10000)
asrextLinBoots=rep(0,10000)
asrsomLinBoots=rep(0,10000)
asranxLinBoots=rep(0,10000)
asrthoLinBoots=rep(0,10000)
asrwitLinBoots=rep(0,10000)
asrattLinBoots=rep(0,10000)
asrrulLinBoots=rep(0,10000)
asraggLinBoots=rep(0,10000)
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$parentPcount)
iMaxVal=max(masterdf$ASRInt)
emaxVal=max(masterdf$ASRExt)
somMaxVal=max(masterdf$ASR_somatic)
anxMaxVal=max(masterdf$ASR_anxdep)
thoMaxVal=max(masterdf$ASR_thought)
witMaxVal=max(masterdf$ASR_withdep)
attMaxVal=max(masterdf$ASR_attention)
rulMaxVal=max(masterdf$ASR_rulebreak)
aggMaxVal=max(masterdf$ASR_aggressive)
asrpMaxVal=max(masterdf$parentPcount)
asriMaxVal=max(masterdf$ASRInt)
asreMaxVal=max(masterdf$ASRExt)
asrsomMaxVal=max(masterdf$ASR_somatic)
asranxMaxVal=max(masterdf$ASR_anxdep)
asrthoMaxVal=max(masterdf$ASR_thought)
asrwitMaxVal=max(masterdf$ASR_withdep)
asrattMaxVal=max(masterdf$ASR_attention)
asrrulMaxVal=max(masterdf$ASR_rulebreak)
asraggMaxVal=max(masterdf$ASR_aggressive)
# predicted derivatives: 10000xncol (+1 because 0-max is 1 longer than max)
pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
intDeriv=matrix(0,nrow=10000,ncol=(iMaxVal+1))
extDeriv=matrix(0,nrow=10000,ncol=(emaxVal+1))
somDeriv=matrix(0,nrow=10000,ncol=(somMaxVal+1))
anxDeriv=matrix(0,nrow=10000,ncol=(anxMaxVal+1))
thoDeriv=matrix(0,nrow=10000,ncol=(thoMaxVal+1))
witDeriv=matrix(0,nrow=10000,ncol=(witMaxVal+1))
attDeriv=matrix(0,nrow=10000,ncol=(attMaxVal+1))
rulDeriv=matrix(0,nrow=10000,ncol=(rulMaxVal+1))
aggDeriv=matrix(0,nrow=10000,ncol=(aggMaxVal+1))
asrpDeriv=matrix(0,nrow=10000,ncol=(asrpMaxVal+1))
asrintDeriv=matrix(0,nrow=10000,ncol=(asriMaxVal+1))
asrextDeriv=matrix(0,nrow=10000,ncol=(asreMaxVal+1))
asrsomDeriv=matrix(0,nrow=10000,ncol=(asrsomMaxVal+1))
asranxDeriv=matrix(0,nrow=10000,ncol=(asranxMaxVal+1))
asrthoDeriv=matrix(0,nrow=10000,ncol=(asrthoMaxVal+1))
asrwitDeriv=matrix(0,nrow=10000,ncol=(asrwitMaxVal+1))
asrattDeriv=matrix(0,nrow=10000,ncol=(asrattMaxVal+1))
asrrulDeriv=matrix(0,nrow=10000,ncol=(asrrulMaxVal+1))
asraggDeriv=matrix(0,nrow=10000,ncol=(asraggMaxVal+1))
# predicted values: set to maximum value for ncol
pFit=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFit=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFit=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFit=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFit=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
attFit=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
asrPFit=matrix(0,nrow=10000,ncol=(asrpMaxVal)+1)
asrintFit=matrix(0,nrow=10000,ncol=(asriMaxVal)+1)
asrextFit=matrix(0,nrow=10000,ncol=(asreMaxVal)+1)
asrsomFit=matrix(0,nrow=10000,ncol=(asrsomMaxVal)+1)
asranxFit=matrix(0,nrow=10000,ncol=(asranxMaxVal)+1)
asrthoFit=matrix(0,nrow=10000,ncol=(asrthoMaxVal)+1)
asrwitFit=matrix(0,nrow=10000,ncol=(asrwitMaxVal)+1)
asrattFit=matrix(0,nrow=10000,ncol=(asrattMaxVal)+1)
asrrulFit=matrix(0,nrow=10000,ncol=(asrrulMaxVal)+1)
asraggFit=matrix(0,nrow=10000,ncol=(asraggMaxVal)+1)
# maximum value in each iteration
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
somMax=rep(0,10000)
anxMax=rep(0,10000)
thoMax=rep(0,10000)
witMax=rep(0,10000)
attMax=rep(0,10000)
rulMax=rep(0,10000)
aggMax=rep(0,10000)
# loop over manual bootstrap
for (b in 1:2000){
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
	bsommax=max(bootSamp$cbcl_scr_syn_somatic_r)
	banxmax=max(bootSamp$cbcl_scr_syn_anxdep_r)
	bthomax=max(bootSamp$cbcl_scr_syn_thought_r)
	bwitmax=max(bootSamp$cbcl_scr_syn_withdep_r)
	battmax=max(bootSamp$cbcl_scr_syn_attention_r)
	brulmax=max(bootSamp$cbcl_scr_syn_rulebreak_r)
	baggmax=max(bootSamp$cbcl_scr_syn_aggressive_r)
	basrpmax=max(bootSamp$parentPcount)
	basrintmax=max(bootSamp$ASRInt)
	basrextmax=max(bootSamp$ASRExt)
	basrsommax=max(bootSamp$ASR_somatic)	
	basranxmax=max(bootSamp$ASR_anxdep)
	basrthomax=max(bootSamp$ASR_thought)
	basrwitmax=max(bootSamp$ASR_withdep)
	basrattmax=max(bootSamp$ASR_attention)
	basrrulmax=max(bootSamp$ASR_rulebreak)
	basraggmax=max(bootSamp$ASR_aggressive)
	######## I FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(g~cbcl_scr_syn_totprob_r+s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age),data=bootSamp)
	plinBoots[b]=summary(pgAge)$s.pv[1]
	intgAge<-bam(g~cbcl_scr_syn_internal_r+s(cbcl_scr_syn_internal_r,m=c(2,0))+s(interview_age),data=bootSamp)
	intlinBoots[b]=summary(intgAge)$s.pv[1]
	extgAge<-bam(g~cbcl_scr_syn_external_r+s(cbcl_scr_syn_external_r,m=c(2,0))+s(interview_age),data=bootSamp)
	extlinBoots[b]=summary(extgAge)$s.pv[1]
	somgAge<-bam(g~cbcl_scr_syn_somatic_r+s(cbcl_scr_syn_somatic_r,m=c(2,0))+s(interview_age),data=bootSamp)
	somLinBoots[b]=summary(somgAge)$s.pv[1]
	anxgAge<-bam(g~cbcl_scr_syn_anxdep_r+s(cbcl_scr_syn_anxdep_r,m=c(2,0))+s(interview_age),data=bootSamp)
	anxLinBoots[b]=summary(anxgAge)$s.pv[1]
	thogAge<-bam(g~cbcl_scr_syn_thought_r+s(cbcl_scr_syn_thought_r,m=c(2,0))+s(interview_age),data=bootSamp)
	thoLinBoots[b]=summary(thogAge)$s.pv[1]
	witgAge<-bam(g~cbcl_scr_syn_withdep_r+s(cbcl_scr_syn_withdep_r,m=c(2,0))+s(interview_age),data=bootSamp)
	witLinBoots[b]=summary(witgAge)$s.pv[1]
	attgAge<-bam(g~cbcl_scr_syn_attention_r+s(cbcl_scr_syn_attention_r,m=c(2,0))+s(interview_age),data=bootSamp)
	attLinBoots[b]=summary(attgAge)$s.pv[1]
	rulgAge<-bam(g~cbcl_scr_syn_rulebreak_r+s(cbcl_scr_syn_rulebreak_r,m=c(2,0))+s(interview_age),data=bootSamp)
	rulLinBoots[b]=summary(rulgAge)$s.pv[1]
	agggAge<-bam(g~cbcl_scr_syn_aggressive_r+s(cbcl_scr_syn_aggressive_r,m=c(2,0))+s(interview_age),data=bootSamp)
	aggLinBoots[b]=summary(agggAge)$s.pv[1]
	#### ASR
	asrpgAge<-bam(g~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age),data=bootSamp)
	asrpLinBoots[b]=summary(asrpgAge)$s.pv[1]
	asrintgAge<-bam(g~ASRInt+s(ASRInt,m=c(2,0))+s(interview_age),data=bootSamp)
	asrintLinBoots[b]=summary(asrintgAge)$s.pv[1]
	asrextgAge<-bam(g~ASRExt+s(ASRExt,m=c(2,0))+s(interview_age),data=bootSamp)
	asrextLinBoots[b]=summary(asrextgAge)$s.pv[1]
	asrsomgAge<-bam(g~ASR_somatic+s(ASR_somatic,m=c(2,0))+s(interview_age),data=bootSamp)
	asrsomLinBoots[b]=summary(asrsomgAge)$s.pv[1]
	asranxgAge<-bam(g~ASR_anxdep+s(ASR_anxdep,m=c(2,0))+s(interview_age),data=bootSamp)
	asranxLinBoots[b]=summary(asranxgAge)$s.pv[1]
	asrthogAge<-bam(g~ASR_thought+s(ASR_thought,m=c(2,0))+s(interview_age),data=bootSamp)
	asrthoLinBoots[b]=summary(asrthogAge)$s.pv[1]
	asrwitgAge<-bam(g~ASR_withdep+s(ASR_withdep,m=c(2,0))+s(interview_age),data=bootSamp)
	asrwitLinBoots[b]=summary(asrwitgAge)$s.pv[1]
	asrattgAge<-bam(g~ASR_attention+s(ASR_attention,m=c(2,0))+s(interview_age),data=bootSamp)
	asrattLinBoots[b]=summary(asrattgAge)$s.pv[1]
	asrrulgAge<-bam(g~ASR_rulebreak+s(ASR_rulebreak,m=c(2,0))+s(interview_age),data=bootSamp)
	asrrulLinBoots[b]=summary(asrrulgAge)$s.pv[1]
	asragggAge<-bam(g~ASR_aggressive+s(ASR_aggressive,m=c(2,0))+s(interview_age),data=bootSamp)
	asraggLinBoots[b]=summary(asragggAge)$s.pv[1]
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable, add asr and predict on asr
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age)+s(parentPcount),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age)+s(ASRInt),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age)+s(ASRExt),data=bootSamp)
	somgAge<-bam(g~s(cbcl_scr_syn_somatic_r)+s(interview_age)+s(ASR_somatic),data=bootSamp)
	anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r)+s(interview_age)+s(ASR_anxdep),data=bootSamp)
	thogAge<-bam(g~s(cbcl_scr_syn_thought_r)+s(interview_age)+s(ASR_thought),data=bootSamp)
	witgAge<-bam(g~s(cbcl_scr_syn_withdep_r)+s(interview_age)+s(ASR_withdep),data=bootSamp)
	attgAge<-bam(g~s(cbcl_scr_syn_attention_r)+s(interview_age)+s(ASR_attention),data=bootSamp)
	rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r)+s(interview_age)+s(ASR_rulebreak),data=bootSamp)
	agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r)+s(interview_age)+s(ASR_aggressive),data=bootSamp)
	#### ASR
	asrpgAge<-bam(g~s(parentPcount)+s(interview_age),data=bootSamp)
	asrintgAge<-bam(g~s(ASRInt)+s(interview_age),data=bootSamp)
	asrextgAge<-bam(g~s(ASRExt)+s(interview_age),data=bootSamp)
	asrsomgAge<-bam(g~s(ASR_somatic)+s(interview_age),data=bootSamp)
	asranxgAge<-bam(g~s(ASR_anxdep)+s(interview_age),data=bootSamp)
	asrthogAge<-bam(g~s(ASR_thought)+s(interview_age),data=bootSamp)
	asrwitgAge<-bam(g~s(ASR_withdep)+s(interview_age),data=bootSamp)
	asrattgAge<-bam(g~s(ASR_attention)+s(interview_age),data=bootSamp)
	asrrulgAge<-bam(g~s(ASR_rulebreak)+s(interview_age),data=bootSamp)
	asragggAge<-bam(g~s(ASR_aggressive)+s(interview_age),data=bootSamp)
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(0:bpmax)
	eachIntcount=seq(0:bimax)
	eachExtcount=seq(0:bemax)
	eachSomcount=seq(0:bsommax)
	eachAnxcount=seq(0:banxmax)
	eachThocount=seq(0:bthomax)
	eachWitcount=seq(0:bwitmax)
	eachAttcount=seq(0:battmax)
	eachRulcount=seq(0:brulmax)
	eachAggcount=seq(0:baggmax)
	eachasrPcount=seq(0:basrpmax)
	eachasrIntcount=seq(0:basrintmax)
	eachasrExtcount=seq(0:basrextmax)
	eachasrSomcount=seq(0:basrsommax)
	eachasrAnxcount=seq(0:basranxmax)
	eachasrThocount=seq(0:basrthomax)
	eachasrWitcount=seq(0:basrwitmax)
	eachasrAttcount=seq(0:basrattmax)
	eachasrRulcount=seq(0:basrrulmax)
	eachasrAggcount=seq(0:basraggmax)
	# set age to to median for predict df, also set child symptom score to median for predict df
	#####################################
	predictDFp=data.frame(eachasrPcount,rep(median(bootSamp$cbcl_scr_syn_totprob_r),(basrpmax+1)),rep(median(bootSamp$interview_age),(basrpmax+1)))
	predictDFint=data.frame(eachasrIntcount,rep(median(bootSamp$cbcl_scr_syn_internal_r),(basrintmax+1)),rep(median(bootSamp$interview_age),(basrintmax+1)))
	predictDFext=data.frame(eachasrExtcount,rep(median(bootSamp$cbcl_scr_syn_external_r),(basrextmax+1)),rep(median(bootSamp$interview_age),(basrextmax+1)))
	predictDFsom=data.frame(eachasrSomcount,rep(median(bootSamp$cbcl_scr_syn_somatic_r),(basrsommax+1)),rep(median(bootSamp$interview_age),(basrsommax+1)))
	predictDFanx=data.frame(eachasrAnxcount,rep(median(bootSamp$cbcl_scr_syn_anxdep_r),(basranxmax+1)),rep(median(bootSamp$interview_age),(basranxmax+1)))
	predictDFtho=data.frame(eachasrThocount,rep(median(bootSamp$cbcl_scr_syn_thought_r),(basrthomax+1)),rep(median(bootSamp$interview_age),(basrthomax+1)))
	predictDFwit=data.frame(eachasrWitcount,rep(median(bootSamp$cbcl_scr_syn_withdep_r),(basrwitmax+1)),rep(median(bootSamp$interview_age),(basrwitmax+1)))
	predictDFatt=data.frame(eachasrAttcount,rep(median(bootSamp$cbcl_scr_syn_attention_r),(basrattmax+1)),rep(median(bootSamp$interview_age),(basrattmax+1)))
	predictDFrul=data.frame(eachasrRulcount,rep(median(bootSamp$cbcl_scr_syn_rulebreak_r),(basrrulmax+1)),rep(median(bootSamp$interview_age),(basrrulmax+1)))
	predictDFagg=data.frame(eachasrAggcount,rep(median(bootSamp$cbcl_scr_syn_aggressive_r),(basraggmax+1)),rep(median(bootSamp$interview_age),(basraggmax+1)))
	predictDFasrp=data.frame(eachasrPcount,rep(median(bootSamp$interview_age),(basrpmax+1)))
	predictDFasrint=data.frame(eachasrIntcount,rep(median(bootSamp$interview_age),(basrintmax+1)))
	predictDFasrext=data.frame(eachasrExtcount,rep(median(bootSamp$interview_age),(basrextmax+1)))
	predictDFasrsom=data.frame(eachasrSomcount,rep(median(bootSamp$interview_age),(basrsommax+1)))
	predictDFasranx=data.frame(eachasrAnxcount,rep(median(bootSamp$interview_age),(basranxmax+1)))
	predictDFasrtho=data.frame(eachasrThocount,rep(median(bootSamp$interview_age),(basrthomax+1)))
	predictDFasrwit=data.frame(eachasrWitcount,rep(median(bootSamp$interview_age),(basrwitmax+1)))
	predictDFasratt=data.frame(eachasrAttcount,rep(median(bootSamp$interview_age),(basrattmax+1)))
	predictDFasrrul=data.frame(eachasrRulcount,rep(median(bootSamp$interview_age),(basrrulmax+1)))
	predictDFasragg=data.frame(eachasrAggcount,rep(median(bootSamp$interview_age),(basraggmax+1)))
	# set colnames so predict can work
	colnames(predictDFp)=c('parentPcount','cbcl_scr_syn_totprob_r','interview_age')
	colnames(predictDFint)=c('ASRInt','cbcl_scr_syn_internal_r','interview_age')
	colnames(predictDFext)=c('ASRExt','cbcl_scr_syn_external_r','interview_age')
	colnames(predictDFsom)=c('ASR_somatic','cbcl_scr_syn_somatic_r','interview_age')
	colnames(predictDFanx)=c('ASR_anxdep','cbcl_scr_syn_anxdep_r','interview_age')
	colnames(predictDFtho)=c('ASR_thought','cbcl_scr_syn_thought_r','interview_age')
	colnames(predictDFwit)=c('ASR_withdep','cbcl_scr_syn_withdep_r','interview_age')
	colnames(predictDFatt)=c('ASR_attention','cbcl_scr_syn_attention_r','interview_age')
	colnames(predictDFrul)=c('ASR_rulebreak','cbcl_scr_syn_rulebreak_r','interview_age')
	colnames(predictDFagg)=c('ASR_aggressive','cbcl_scr_syn_aggressive_r','interview_age')
	colnames(predictDFasrp)=c('parentPcount','interview_age')
	colnames(predictDFasrint)=c('ASRInt','interview_age')
	colnames(predictDFasrext)=c('ASRExt','interview_age')
	colnames(predictDFasrsom)=c('ASR_somatic','interview_age')
	colnames(predictDFasranx)=c('ASR_anxdep','interview_age')
	colnames(predictDFasrtho)=c('ASR_thought','interview_age')
	colnames(predictDFasrwit)=c('ASR_withdep','interview_age')
	colnames(predictDFasratt)=c('ASR_attention','interview_age')
	colnames(predictDFasrrul)=c('ASR_rulebreak','interview_age')
	colnames(predictDFasragg)=c('ASR_aggressive','interview_age')
	# predict
	forFitP=predict(pgAge,predictDFp)
	forFitInt=predict(intgAge,predictDFint)
	forFitExt=predict(extgAge,predictDFext)
	forFitSom=predict(somgAge,predictDFsom)
	forFitAnx=predict(anxgAge,predictDFanx)
	forFitTho=predict(thogAge,predictDFtho)
	forFitWit=predict(witgAge,predictDFwit)
	forFitAtt=predict(attgAge,predictDFatt)
	forFitRul=predict(rulgAge,predictDFrul)
	forFitAgg=predict(agggAge,predictDFagg)
	forFitasrP=predict(asrpgAge,predictDFasrp)
	forFitasrInt=predict(asrintgAge,predictDFasrint)
	forFitasrExt=predict(asrextgAge,predictDFasrext)
	forFitasrSom=predict(asrsomgAge,predictDFasrsom)
	forFitasrAnx=predict(asranxgAge,predictDFasranx)
	forFitasrTho=predict(asrthogAge,predictDFasrtho)
	forFitasrWit=predict(asrwitgAge,predictDFasrwit)
	forFitasrAtt=predict(asrattgAge,predictDFasratt)
	forFitasrRul=predict(asrrulgAge,predictDFasrrul)
	forFitasrAgg=predict(asragggAge,predictDFasragg)
	# print out fit
	pFit[b,1:(basrpmax+1)]=forFitP
	intFit[b,1:(basrintmax+1)]=forFitInt
	extFit[b,1:(basrextmax+1)]=forFitExt
	somFit[b,1:(basrsommax+1)]=forFitSom
	anxFit[b,1:(basranxmax+1)]=forFitAnx
	thoFit[b,1:(basrthomax+1)]=forFitTho
	witFit[b,1:(basrwitmax+1)]=forFitWit
	attFit[b,1:(basrattmax+1)]=forFitAtt
	rulFit[b,1:(basrrulmax+1)]=forFitRul
	aggFit[b,1:(basraggmax+1)]=forFitAgg
	asrPFit[b,1:(basrpmax+1)]=forFitasrP
	asrintFit[b,1:(basrintmax+1)]=forFitasrInt
	asrextFit[b,1:(basrextmax+1)]=forFitasrExt
	asrsomFit[b,1:(basrsommax+1)]=forFitasrSom
	asranxFit[b,1:(basranxmax+1)]=forFitasrAnx
	asrthoFit[b,1:(basrthomax+1)]=forFitasrTho
	asrwitFit[b,1:(basrwitmax+1)]=forFitasrWit
	asrattFit[b,1:(basrattmax+1)]=forFitasrAtt
	asrrulFit[b,1:(basrrulmax+1)]=forFitasrRul
	asraggFit[b,1:(basraggmax+1)]=forFitasrAgg
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(parentPcount)',partial_match = TRUE,n=(basrpmax+1))
	forSplineint=derivatives(intgAge,term='s(ASRInt)',partial_match = TRUE,n=(basrintmax+1))
	forSplineext=derivatives(extgAge,term='s(ASRExt)',partial_match = TRUE,n=(basrextmax+1))
	forSplinesom=derivatives(somgAge,term='s(ASR_somatic)',partial_match = TRUE,n=(basrsommax+1))
	forSplineanx=derivatives(anxgAge,term='s(ASR_anxdep)',partial_match = TRUE,n=(basranxmax+1))
	forSplinetho=derivatives(thogAge,term='s(ASR_thought)',partial_match = TRUE,n=(basrthomax+1))
	forSplinewit=derivatives(witgAge,term='s(ASR_withdep)',partial_match = TRUE,n=(basrwitmax+1))
	forSplineatt=derivatives(attgAge,term='s(ASR_attention)',partial_match = TRUE,n=(basrattmax+1))
	forSplinerul=derivatives(rulgAge,term='s(ASR_rulebreak)',partial_match = TRUE,n=(basrrulmax+1))
	forSplineagg=derivatives(agggAge,term='s(ASR_aggressive)',partial_match = TRUE,n=(basraggmax+1))
	# asr
	forSplineasrp=derivatives(asrpgAge,term='s(parentPcount)',partial_match = TRUE,n=(basrpmax+1))
	forSplineasrint=derivatives(asrintgAge,term='s(ASRInt)',partial_match = TRUE,n=(basrintmax+1))
	forSplineasrext=derivatives(asrextgAge,term='s(ASRExt)',partial_match = TRUE,n=(basrextmax+1))
	forSplineasrsom=derivatives(asrsomgAge,term='s(ASR_somatic)',partial_match = TRUE,n=(basrsommax+1))
	forSplineasranx=derivatives(asranxgAge,term='s(ASR_anxdep)',partial_match = TRUE,n=(basranxmax+1))
	forSplineasrtho=derivatives(asrthogAge,term='s(ASR_thought)',partial_match = TRUE,n=(basrthomax+1))
	forSplineasrwit=derivatives(asrwitgAge,term='s(ASR_withdep)',partial_match = TRUE,n=(basrwitmax+1))
	forSplineasratt=derivatives(asrattgAge,term='s(ASR_attention)',partial_match = TRUE,n=(basrattmax+1))
	forSplineasrrul=derivatives(asrrulgAge,term='s(ASR_rulebreak)',partial_match = TRUE,n=(basrrulmax+1))
	forSplineasragg=derivatives(asragggAge,term='s(ASR_aggressive)',partial_match = TRUE,n=(basraggmax+1))
	# print out fit derivatives
	pDeriv[b,1:(basrpmax+1)]=forSplinep$derivative
	intDeriv[b,1:(basrintmax+1)]=forSplineint$derivative
	extDeriv[b,1:(basrextmax+1)]=forSplineext$derivative
	somDeriv[b,1:(basrsommax+1)]=forSplinesom$derivative
	anxDeriv[b,1:(basranxmax+1)]=forSplineanx$derivative
	thoDeriv[b,1:(basrthomax+1)]=forSplinetho$derivative
	witDeriv[b,1:(basrwitmax+1)]=forSplinewit$derivative
	attDeriv[b,1:(basrattmax+1)]=forSplineatt$derivative
	rulDeriv[b,1:(basrrulmax+1)]=forSplinerul$derivative
	aggDeriv[b,1:(basraggmax+1)]=forSplineagg$derivative
	# asr
	asrpDeriv[b,1:(basrpmax+1)]=forSplineasrp$derivative
	asrintDeriv[b,1:(basrintmax+1)]=forSplineasrint$derivative
	asrextDeriv[b,1:(basrextmax+1)]=forSplineasrext$derivative
	asrsomDeriv[b,1:(basrsommax+1)]=forSplineasrsom$derivative
	asranxDeriv[b,1:(basranxmax+1)]=forSplineasranx$derivative
	asrthoDeriv[b,1:(basrthomax+1)]=forSplineasrtho$derivative
	asrwitDeriv[b,1:(basrwitmax+1)]=forSplineasrwit$derivative
	asrattDeriv[b,1:(basrattmax+1)]=forSplineasratt$derivative
	asrrulDeriv[b,1:(basrrulmax+1)]=forSplineasrrul$derivative
	asraggDeriv[b,1:(basraggmax+1)]=forSplineasragg$derivative
}
# SAVEOUT
# save out version with all cbcl and asr linboots
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,somLinBoots,anxLinBoots,thoLinBoots,witLinBoots,attLinBoots,rulLinBoots,aggLinBoots,asrpLinBoots,asrintLinBoots,asrextLinBoots,asrsomLinBoots,asranxLinBoots,asrthoLinBoots,asrwitLinBoots,asrattLinBoots,asrrulLinBoots,asraggLinBoots)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots_asr2k.rds')
# save out version with all cbcl and asr derivs
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,attDeriv,rulDeriv,aggDeriv,asrpDeriv,asrintDeriv,asrextDeriv,asrsomDeriv,asranxDeriv,asrthoDeriv,asrwitDeriv,asrattDeriv,asrrulDeriv,asraggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_asr2k.rds')
# save out version with all cbcl and asr fits
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,attFit,rulFit,aggFit,asrPFit,asrintFit,asrextFit,asrsomFit,asranxFit,asrthoFit,asrwitFit,asrattFit,asrrulFit,asraggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_asr2k.rds')
print('done with g~p fit bootstrapping!')
