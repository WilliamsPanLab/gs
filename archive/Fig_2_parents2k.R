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
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','ASRAnxDepr','ASRWithdrawn','ASRSomatic','ASRThought','ASRAttn','ASRAggr','ASRIntrusive','ASRRulB','ASRInt','ASRExt','g','subjectkey','interview_age')]
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
masterdf$ASR_intrusive=as.numeric(masterdf$ASRIntrusive)
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
asrintrLinBoots=rep(0,10000)
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
asrintrMaxVal=max(masterdf$ASR_intrusive)
# predicted derivatives: 10000xncol
pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
somDeriv=matrix(0,nrow=10000,ncol=somMaxVal)
anxDeriv=matrix(0,nrow=10000,ncol=anxMaxVal)
thoDeriv=matrix(0,nrow=10000,ncol=thoMaxVal)
witDeriv=matrix(0,nrow=10000,ncol=witMaxVal)
attDeriv=matrix(0,nrow=10000,ncol=attMaxVal)
rulDeriv=matrix(0,nrow=10000,ncol=rulMaxVal)
aggDeriv=matrix(0,nrow=10000,ncol=aggMaxVal)
asrpDeriv=matrix(0,nrow=10000,ncol=asrpMaxVal)
asrintDeriv=matrix(0,nrow=10000,ncol=asriMaxVal)
asrextDeriv=matrix(0,nrow=10000,ncol=asreMaxVal)
asrsomDeriv=matrix(0,nrow=10000,ncol=asrsomMaxVal)
asranxDeriv=matrix(0,nrow=10000,ncol=asranxMaxVal)
asrthoDeriv=matrix(0,nrow=10000,ncol=asrthoMaxVal)
asrwitDeriv=matrix(0,nrow=10000,ncol=asrwitMaxVal)
asrattDeriv=matrix(0,nrow=10000,ncol=asrattMaxVal)
asrrulDeriv=matrix(0,nrow=10000,ncol=asrrulMaxVal)
asraggDeriv=matrix(0,nrow=10000,ncol=asraggMaxVal)
asrintrDeriv=matrix(0,nrow=10000,ncol=asrintrMaxVal)
# predicted values: set to maximum value for ncol
pFit=matrix(0,nrow=10000,ncol=pMaxVal)
intFit=matrix(0,nrow=10000,ncol=iMaxVal)
extFit=matrix(0,nrow=10000,ncol=emaxVal)
somFit=matrix(0,nrow=10000,ncol=somMaxVal)
anxFit=matrix(0,nrow=10000,ncol=anxMaxVal)
thoFit=matrix(0,nrow=10000,ncol=thoMaxVal)
witFit=matrix(0,nrow=10000,ncol=witMaxVal)
attFit=matrix(0,nrow=10000,ncol=attMaxVal)
rulFit=matrix(0,nrow=10000,ncol=rulMaxVal)
aggFit=matrix(0,nrow=10000,ncol=aggMaxVal)
asrPFit=matrix(0,nrow=10000,ncol=asrpMaxVal)
asrintFit=matrix(0,nrow=10000,ncol=asriMaxVal)
asrextFit=matrix(0,nrow=10000,ncol=asreMaxVal)
asrsomFit=matrix(0,nrow=10000,ncol=asrsomMaxVal)
asranxFit=matrix(0,nrow=10000,ncol=asranxMaxVal)
asrthoFit=matrix(0,nrow=10000,ncol=asrthoMaxVal)
asrwitFit=matrix(0,nrow=10000,ncol=asrwitMaxVal)
asrattFit=matrix(0,nrow=10000,ncol=asrattMaxVal)
asrrulFit=matrix(0,nrow=10000,ncol=asrrulMaxVal)
asraggFit=matrix(0,nrow=10000,ncol=asraggMaxVal)
asrintrFit=matrix(0,nrow=10000,ncol=asrintrMaxVal)
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
intrMax=rep(0,10000)
set.seed(1)
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
	basrintrmax=max(bootSamp$ASR_intrusive)
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
	asrintrgAge<-bam(g~ASR_intrusive+s(ASR_intrusive,m=c(2,0))+s(interview_age),data=bootSamp)
	asrintrLinBoots[b]=summary(asrintrgAge)$s.pv[1]
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
	asrintrgAge<-bam(g~s(ASR_intrusive)+s(interview_age),data=bootSamp)
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(1:bpmax)
	eachIntcount=seq(1:bimax)
	eachExtcount=seq(1:bemax)
	eachSomcount=seq(1:bsommax)
	eachAnxcount=seq(1:banxmax)
	eachThocount=seq(1:bthomax)
	eachWitcount=seq(1:bwitmax)
	eachAttcount=seq(1:battmax)
	eachRulcount=seq(1:brulmax)
	eachAggcount=seq(1:baggmax)
	eachasrPcount=seq(1:basrpmax)
	eachasrIntcount=seq(1:basrintmax)
	eachasrExtcount=seq(1:basrextmax)
	eachasrSomcount=seq(1:basrsommax)
	eachasrAnxcount=seq(1:basranxmax)
	eachasrThocount=seq(1:basrthomax)
	eachasrWitcount=seq(1:basrwitmax)
	eachasrAttcount=seq(1:basrattmax)
	eachasrRulcount=seq(1:basrrulmax)
	eachasrAggcount=seq(1:basraggmax)
	eachasrIntrcount=seq(1:basrintrmax)
	# set age to to median for predict df, also set child symptom score to median for predict df
	#####################################
	predictDFp=data.frame(eachasrPcount,rep(median(bootSamp$cbcl_scr_syn_totprob_r),basrpmax),rep(median(bootSamp$interview_age),basrpmax))
	predictDFint=data.frame(eachasrIntcount,rep(median(bootSamp$cbcl_scr_syn_internal_r),basrintmax),rep(median(bootSamp$interview_age),basrintmax))
	predictDFext=data.frame(eachasrExtcount,rep(median(bootSamp$cbcl_scr_syn_external_r),basrextmax),rep(median(bootSamp$interview_age),basrextmax))
	predictDFsom=data.frame(eachasrSomcount,rep(median(bootSamp$cbcl_scr_syn_somatic_r),basrsommax),rep(median(bootSamp$interview_age),basrsommax))
	predictDFanx=data.frame(eachasrAnxcount,rep(median(bootSamp$cbcl_scr_syn_anxdep_r),basranxmax),rep(median(bootSamp$interview_age),basranxmax))
	predictDFtho=data.frame(eachasrThocount,rep(median(bootSamp$cbcl_scr_syn_thought_r),basrthomax),rep(median(bootSamp$interview_age),basrthomax))
	predictDFwit=data.frame(eachasrWitcount,rep(median(bootSamp$cbcl_scr_syn_withdep_r),basrwitmax),rep(median(bootSamp$interview_age),basrwitmax))
	predictDFatt=data.frame(eachasrAttcount,rep(median(bootSamp$cbcl_scr_syn_attention_r),basrattmax),rep(median(bootSamp$interview_age),basrattmax))
	predictDFrul=data.frame(eachasrRulcount,rep(median(bootSamp$cbcl_scr_syn_rulebreak_r),basrrulmax),rep(median(bootSamp$interview_age),basrrulmax))
	predictDFagg=data.frame(eachasrAggcount,rep(median(bootSamp$cbcl_scr_syn_aggressive_r),basraggmax),rep(median(bootSamp$interview_age),basraggmax))
	predictDFasrp=data.frame(eachasrPcount,rep(median(bootSamp$interview_age),basrpmax))
	predictDFasrint=data.frame(eachasrIntcount,rep(median(bootSamp$interview_age),basrintmax))
	predictDFasrext=data.frame(eachasrExtcount,rep(median(bootSamp$interview_age),basrextmax))
	predictDFasrsom=data.frame(eachasrSomcount,rep(median(bootSamp$interview_age),basrsommax))
	predictDFasranx=data.frame(eachasrAnxcount,rep(median(bootSamp$interview_age),basranxmax))
	predictDFasrtho=data.frame(eachasrThocount,rep(median(bootSamp$interview_age),basrthomax))
	predictDFasrwit=data.frame(eachasrWitcount,rep(median(bootSamp$interview_age),basrwitmax))
	predictDFasratt=data.frame(eachasrAttcount,rep(median(bootSamp$interview_age),basrattmax))
	predictDFasrrul=data.frame(eachasrRulcount,rep(median(bootSamp$interview_age),basrrulmax))
	predictDFasragg=data.frame(eachasrAggcount,rep(median(bootSamp$interview_age),basraggmax))
	predictDFasrintr=data.frame(eachasrIntrcount,rep(median(bootSamp$interview_age),basrintrmax))
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
	colnames(predictDFasrintr)=c('ASR_intrusive','interview_age')
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
	forFitasrIntr=predict(asrintrgAge,predictDFasrintr)
	# print out fit
	pFit[b,1:basrpmax]=forFitP
	intFit[b,1:basrintmax]=forFitInt
	extFit[b,1:basrextmax]=forFitExt
	somFit[b,1:basrsommax]=forFitSom
	anxFit[b,1:basranxmax]=forFitAnx
	thoFit[b,1:basrthomax]=forFitTho
	witFit[b,1:basrwitmax]=forFitWit
	attFit[b,1:basrattmax]=forFitAtt
	rulFit[b,1:basrrulmax]=forFitRul
	aggFit[b,1:basraggmax]=forFitAgg
	asrPFit[b,1:basrpmax]=forFitasrP
	asrintFit[b,1:basrintmax]=forFitasrInt
	asrextFit[b,1:basrextmax]=forFitasrExt
	asrsomFit[b,1:basrsommax]=forFitasrSom
	asranxFit[b,1:basranxmax]=forFitasrAnx
	asrthoFit[b,1:basrthomax]=forFitasrTho
	asrwitFit[b,1:basrwitmax]=forFitasrWit
	asrattFit[b,1:basrattmax]=forFitasrAtt
	asrrulFit[b,1:basrrulmax]=forFitasrRul
	asraggFit[b,1:basraggmax]=forFitasrAgg
	asrintrFit[b,1:basrintrmax]=forFitasrIntr
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(parentPcount)',partial_match = TRUE,n=basrpmax)
	forSplineint=derivatives(intgAge,term='s(ASRInt)',partial_match = TRUE,n=basrintmax)
	forSplineext=derivatives(extgAge,term='s(ASRExt)',partial_match = TRUE,n=basrextmax)
	forSplinesom=derivatives(somgAge,term='s(ASR_somatic)',partial_match = TRUE,n=basrsommax)
	forSplineanx=derivatives(anxgAge,term='s(ASR_anxdep)',partial_match = TRUE,n=basranxmax)
	forSplinetho=derivatives(thogAge,term='s(ASR_thought)',partial_match = TRUE,n=basrthomax)
	forSplinewit=derivatives(witgAge,term='s(ASR_withdep)',partial_match = TRUE,n=basrwitmax)
	forSplineatt=derivatives(attgAge,term='s(ASR_attention)',partial_match = TRUE,n=basrattmax)
	forSplinerul=derivatives(rulgAge,term='s(ASR_rulebreak)',partial_match = TRUE,n=basrrulmax)
	forSplineagg=derivatives(agggAge,term='s(ASR_aggressive)',partial_match = TRUE,n=basraggmax)
	# asr
	forSplineasrp=derivatives(asrpgAge,term='s(parentPcount)',partial_match = TRUE,n=basrpmax)
	forSplineasrint=derivatives(asrintgAge,term='s(ASRInt)',partial_match = TRUE,n=basrintmax)
	forSplineasrext=derivatives(asrextgAge,term='s(ASRExt)',partial_match = TRUE,n=basrextmax)
	forSplineasrsom=derivatives(asrsomgAge,term='s(ASR_somatic)',partial_match = TRUE,n=basrsommax)
	forSplineasranx=derivatives(asranxgAge,term='s(ASR_anxdep)',partial_match = TRUE,n=basranxmax)
	forSplineasrtho=derivatives(asrthogAge,term='s(ASR_thought)',partial_match = TRUE,n=basrthomax)
	forSplineasrwit=derivatives(asrwitgAge,term='s(ASR_withdep)',partial_match = TRUE,n=basrwitmax)
	forSplineasratt=derivatives(asrattgAge,term='s(ASR_attention)',partial_match = TRUE,n=basrattmax)
	forSplineasrrul=derivatives(asrrulgAge,term='s(ASR_rulebreak)',partial_match = TRUE,n=basrrulmax)
	forSplineasragg=derivatives(asragggAge,term='s(ASR_aggressive)',partial_match = TRUE,n=basraggmax)
	forSplineasrintr=derivatives(asrintrgAge,term='s(ASR_intrusive)',partial_match = TRUE,n=basrintrmax)
	# print out fit derivatives
	pDeriv[b,1:basrpmax]=forSplinep$derivative
	intDeriv[b,1:basrintmax]=forSplineint$derivative
	extDeriv[b,1:basrextmax]=forSplineext$derivative
	somDeriv[b,1:basrsommax]=forSplinesom$derivative
	anxDeriv[b,1:basranxmax]=forSplineanx$derivative
	thoDeriv[b,1:basrthomax]=forSplinetho$derivative
	witDeriv[b,1:basrwitmax]=forSplinewit$derivative
	attDeriv[b,1:basrattmax]=forSplineatt$derivative
	rulDeriv[b,1:basrrulmax]=forSplinerul$derivative
	aggDeriv[b,1:basraggmax]=forSplineagg$derivative
	# asr
	asrpDeriv[b,1:basrpmax]=forSplineasrp$derivative
	asrintDeriv[b,1:basrintmax]=forSplineasrint$derivative
	asrextDeriv[b,1:basrextmax]=forSplineasrext$derivative
	asrsomDeriv[b,1:basrsommax]=forSplineasrsom$derivative
	asranxDeriv[b,1:basranxmax]=forSplineasranx$derivative
	asrthoDeriv[b,1:basrthomax]=forSplineasrtho$derivative
	asrwitDeriv[b,1:basrwitmax]=forSplineasrwit$derivative
	asrattDeriv[b,1:basrattmax]=forSplineasratt$derivative
	asrrulDeriv[b,1:basrrulmax]=forSplineasrrul$derivative
	asraggDeriv[b,1:basraggmax]=forSplineasragg$derivative
	asrintrDeriv[b,1:basrintrmax]=forSplineasrintr$derivative
}
# SAVEOUT
# save out version with all cbcl and asr linboots
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,somLinBoots,anxLinBoots,thoLinBoots,witLinBoots,attLinBoots,rulLinBoots,aggLinBoots,asrpLinBoots,asrintLinBoots,asrextLinBoots,asrsomLinBoots,asranxLinBoots,asrthoLinBoots,asrwitLinBoots,asrattLinBoots,asrrulLinBoots,asraggLinBoots,asrintrLinBoots)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots_asr2k.rds')
# save out version with all cbcl and asr derivs
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,attDeriv,rulDeriv,aggDeriv,asrpDeriv,asrintDeriv,asrextDeriv,asrsomDeriv,asranxDeriv,asrthoDeriv,asrwitDeriv,asrattDeriv,asrrulDeriv,asraggDeriv,asrintrDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_asr2k.rds')
# save out version with all cbcl and asr fits
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,attFit,rulFit,aggFit,asrPFit,asrintFit,asrextFit,asrsomFit,asranxFit,asrthoFit,asrwitFit,asrattFit,asrrulFit,asraggFit,asrintrFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_asr2k.rds')
print('done with g~p fit bootstrapping!')
