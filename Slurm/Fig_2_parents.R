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
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','ASRAnxDepr','ASRWithdrawn','ASRSomatic','ASRThought','ASRAttn','ASRAggr','ASRRulB','ASRInt','ASRExt','g','subjectkey','interview_age','Pov_v2')]
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
# initialize version of each for poverty and nonpoverty groups
pFitPov=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
attFitPov=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
asrPFitPov=matrix(0,nrow=10000,ncol=(asrpMaxVal)+1)
asrintFitPov=matrix(0,nrow=10000,ncol=(asriMaxVal)+1)
asrextFitPov=matrix(0,nrow=10000,ncol=(asreMaxVal)+1)
asrsomFitPov=matrix(0,nrow=10000,ncol=(asrsomMaxVal)+1)
asranxFitPov=matrix(0,nrow=10000,ncol=(asranxMaxVal)+1)
asrthoFitPov=matrix(0,nrow=10000,ncol=(asrthoMaxVal)+1)
asrwitFitPov=matrix(0,nrow=10000,ncol=(asrwitMaxVal)+1)
asrattFitPov=matrix(0,nrow=10000,ncol=(asrattMaxVal)+1)
asrrulFitPov=matrix(0,nrow=10000,ncol=(asrrulMaxVal)+1)
asraggFitPov=matrix(0,nrow=10000,ncol=(asraggMaxVal)+1)
pFitNonPov=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitNonPov=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitNonPov=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitNonPov=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitNonPov=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitNonPov=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitNonPov=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
attFitNonPov=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitNonPov=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitNonPov=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
asrPFitNonPov=matrix(0,nrow=10000,ncol=(asrpMaxVal)+1)
asrintFitNonPov=matrix(0,nrow=10000,ncol=(asriMaxVal)+1)
asrextFitNonPov=matrix(0,nrow=10000,ncol=(asreMaxVal)+1)
asrsomFitNonPov=matrix(0,nrow=10000,ncol=(asrsomMaxVal)+1)
asranxFitNonPov=matrix(0,nrow=10000,ncol=(asranxMaxVal)+1)
asrthoFitNonPov=matrix(0,nrow=10000,ncol=(asrthoMaxVal)+1)
asrwitFitNonPov=matrix(0,nrow=10000,ncol=(asrwitMaxVal)+1)
asrattFitNonPov=matrix(0,nrow=10000,ncol=(asrattMaxVal)+1)
asrrulFitNonPov=matrix(0,nrow=10000,ncol=(asrrulMaxVal)+1)
asraggFitNonPov=matrix(0,nrow=10000,ncol=(asraggMaxVal)+1)
# initialize difference in AIC for each asr model with and without poverty interaction
asrPDiff=rep(0,10000)
asrintDiff=rep(0,10000)
asrextDiff=rep(0,10000)
asrsomDiff=rep(0,10000)
asranxDiff=rep(0,10000)
asrthoDiff=rep(0,10000)
asrwitDiff=rep(0,10000)
asrattDiff=rep(0,10000)
asrrulDiff=rep(0,10000)
asraggDiff=rep(0,10000)
# difference in AIC for each asr model with psuedo poverty interaction
asrPDiffPseudo=rep(0,10000)
asrintDiffPseudo=rep(0,10000)
asrextDiffPseudo=rep(0,10000)
asrsomDiffPseudo=rep(0,10000)
asranxDiffPseudo=rep(0,10000)
asrthoDiffPseudo=rep(0,10000)
asrwitDiffPseudo=rep(0,10000)
asrattDiffPseudo=rep(0,10000)
asrrulDiffPseudo=rep(0,10000)
asraggDiffPseudo=rep(0,10000)
# initialize difference in adjusted r^2 for each asr model with and without poverty interaction
asrPDiffAdj=rep(0,10000)
asrintDiffAdj=rep(0,10000)
asrextDiffAdj=rep(0,10000)
asrsomDiffAdj=rep(0,10000)
asranxDiffAdj=rep(0,10000)
asrthoDiffAdj=rep(0,10000)
asrwitDiffAdj=rep(0,10000)
asrattDiffAdj=rep(0,10000)
asrrulDiffAdj=rep(0,10000)
asraggDiffAdj=rep(0,10000)
# initialize difference in adjusted r^2 for each asr model with psuedo poverty interaction
asrPDiffAdjPseudo=rep(0,10000)
asrintDiffAdjPseudo=rep(0,10000)
asrextDiffAdjPseudo=rep(0,10000)
asrsomDiffAdjPseudo=rep(0,10000)
asranxDiffAdjPseudo=rep(0,10000)
asrthoDiffAdjPseudo=rep(0,10000)
asrwitDiffAdjPseudo=rep(0,10000)
asrattDiffAdjPseudo=rep(0,10000)
asrrulDiffAdjPseudo=rep(0,10000)
asraggDiffAdjPseudo=rep(0,10000)
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
# finally, add pov classification (v2)
masterdf$poverty=0
masterdf$income<-as.numeric(masterdf$income)
# poverty now defined in sample construction
masterdf$poverty[masterdf$Pov_v2==1]=1
masterdf$poverty=as.factor(masterdf$poverty)
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
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
        # last pre-step: need to create NULL variables. Use count of Poverty as count of membership to psuedo-groups, but randomly distribute membership
        # get count poverty in this boot
        Povcount=length(which(bootSamp$poverty==1))
        # randomly assign Povcount people to pseudopoverty
        bootSamp$psuedopoverty=0
        bootSamp$psuedopoverty[sample(1:dim(bootSamp)[1],Povcount)]=1
        # make a df that is just the same number of kids in poverty, but actually from the non-poverty group. The point is to see if we can recover g~p slope in a an equivalent number of non-poverty kids
        bootSamp$pseudopoverty2=0
        bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
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
	######## I PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable, add asr and predict on asr
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+s(parentPcount,k=4),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+s(ASRInt,k=4),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+s(ASRExt,k=4),data=bootSamp)
	somgAge<-bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+s(ASR_somatic,k=4),data=bootSamp)
	anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+s(ASR_anxdep,k=4),data=bootSamp)
	thogAge<-bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+s(ASR_thought,k=4),data=bootSamp)
	witgAge<-bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+s(ASR_withdep,k=4),data=bootSamp)
	attgAge<-bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+s(ASR_attention,k=4),data=bootSamp)
	rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+s(ASR_rulebreak,k=4),data=bootSamp)
	agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+s(ASR_aggressive,k=4),data=bootSamp)
	#### ASR
	asrpgAge<-bam(g~s(parentPcount,k=4)+s(interview_age,k=4),data=bootSamp)
	asrintgAge<-bam(g~s(ASRInt,k=4)+s(interview_age,k=4),data=bootSamp)
	asrextgAge<-bam(g~s(ASRExt,k=4)+s(interview_age,k=4),data=bootSamp)
	asrsomgAge<-bam(g~s(ASR_somatic,k=4)+s(interview_age,k=4),data=bootSamp)
	asranxgAge<-bam(g~s(ASR_anxdep,k=4)+s(interview_age,k=4),data=bootSamp)
	asrthogAge<-bam(g~s(ASR_thought,k=4)+s(interview_age,k=4),data=bootSamp)
	asrwitgAge<-bam(g~s(ASR_withdep,k=4)+s(interview_age,k=4),data=bootSamp)
	asrattgAge<-bam(g~s(ASR_attention,k=4)+s(interview_age,k=4),data=bootSamp)
	asrrulgAge<-bam(g~s(ASR_rulebreak,k=4)+s(interview_age,k=4),data=bootSamp)
	asragggAge<-bam(g~s(ASR_aggressive,k=4)+s(interview_age,k=4),data=bootSamp)
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
	####### II PREDICT POVERTY INTERACTIONS #######
	# fit models with standalone poverty term
	asrpgAge_pov=bam(g~s(parentPcount,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrintgAge_pov=bam(g~s(ASRInt,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrextgAge_pov=bam(g~s(ASRExt,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrsomgAge_pov=bam(g~s(ASR_somatic,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asranxgAge_pov=bam(g~s(ASR_anxdep,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrthogAge_pov=bam(g~s(ASR_thought,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrwitgAge_pov=bam(g~s(ASR_withdep,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrattgAge_pov=bam(g~s(ASR_attention,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrrulgAge_pov=bam(g~s(ASR_rulebreak,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asragggAge_pov=bam(g~s(ASR_aggressive,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	# fit versions with poverty interaction on symptoms
	asrpgAge_povint=bam(g~s(parentPcount,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrintgAge_povint=bam(g~s(ASRInt,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrextgAge_povint=bam(g~s(ASRExt,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrsomgAge_povint=bam(g~s(ASR_somatic,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asranxgAge_povint=bam(g~s(ASR_anxdep,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrthogAge_povint=bam(g~s(ASR_thought,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrwitgAge_povint=bam(g~s(ASR_withdep,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrattgAge_povint=bam(g~s(ASR_attention,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asrrulgAge_povint=bam(g~s(ASR_rulebreak,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	asragggAge_povint=bam(g~s(ASR_aggressive,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	# fit versions on psuedopoverty variable
	asrpgAge_povpseudo=bam(g~s(parentPcount,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrintgAge_povpseudo=bam(g~s(ASRInt,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrextgAge_povpseudo=bam(g~s(ASRExt,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrsomgAge_povpseudo=bam(g~s(ASR_somatic,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asranxgAge_povpseudo=bam(g~s(ASR_anxdep,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrthogAge_povpseudo=bam(g~s(ASR_thought,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrwitgAge_povpseudo=bam(g~s(ASR_withdep,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrattgAge_povpseudo=bam(g~s(ASR_attention,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrrulgAge_povpseudo=bam(g~s(ASR_rulebreak,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asragggAge_povpseudo=bam(g~s(ASR_aggressive,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	# fit versions with psuedopoverty interaction on symptoms
	asrpgAge_povpseudoint=bam(g~s(parentPcount,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrintgAge_povpseudoint=bam(g~s(ASRInt,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrextgAge_povpseudoint=bam(g~s(ASRExt,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrsomgAge_povpseudoint=bam(g~s(ASR_somatic,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asranxgAge_povpseudoint=bam(g~s(ASR_anxdep,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrthogAge_povpseudoint=bam(g~s(ASR_thought,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrwitgAge_povpseudoint=bam(g~s(ASR_withdep,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrattgAge_povpseudoint=bam(g~s(ASR_attention,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asrrulgAge_povpseudoint=bam(g~s(ASR_rulebreak,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	asragggAge_povpseudoint=bam(g~s(ASR_aggressive,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	# get AIC differences
	asrPDiff[b]=AIC(asrpgAge_pov)-AIC(asrpgAge_povint)
	asrintDiff[b]=AIC(asrintgAge_pov)-AIC(asrintgAge_povint)
	asrextDiff[b]=AIC(asrextgAge_pov)-AIC(asrextgAge_povint)
	asrsomDiff[b]=AIC(asrsomgAge_pov)-AIC(asrsomgAge_povint)
	asranxDiff[b]=AIC(asranxgAge_pov)-AIC(asranxgAge_povint)
	asrthoDiff[b]=AIC(asrthogAge_pov)-AIC(asrthogAge_povint)
	asrwitDiff[b]=AIC(asrwitgAge_pov)-AIC(asrwitgAge_povint)
	asrattDiff[b]=AIC(asrattgAge_pov)-AIC(asrattgAge_povint)
	asrrulDiff[b]=AIC(asrrulgAge_pov)-AIC(asrrulgAge_povint)
	asraggDiff[b]=AIC(asragggAge_pov)-AIC(asragggAge_povint)
	asrPDiff_psuedo[b]=AIC(asrpgAge_povpseudo)-AIC(asrpgAge_povpseudoint)
	asrintDiff_psuedo[b]=AIC(asrintgAge_povpseudo)-AIC(asrintgAge_povpseudoint)
	asrextDiff_psuedo[b]=AIC(asrextgAge_povpseudo)-AIC(asrextgAge_povpseudoint)
	asrsomDiff_psuedo[b]=AIC(asrsomgAge_povpseudo)-AIC(asrsomgAge_povpseudoint)
	asranxDiff_psuedo[b]=AIC(asranxgAge_povpseudo)-AIC(asranxgAge_povpseudoint)
	asrthoDiff_psuedo[b]=AIC(asrthogAge_povpseudo)-AIC(asrthogAge_povpseudoint)
	asrwitDiff_psuedo[b]=AIC(asrwitgAge_povpseudo)-AIC(asrwitgAge_povpseudoint)
	asrattDiff_psuedo[b]=AIC(asrattgAge_povpseudo)-AIC(asrattgAge_povpseudoint)
	asrrulDiff_psuedo[b]=AIC(asrrulgAge_povpseudo)-AIC(asrrulgAge_povpseudoint)
	asraggDiff_psuedo[b]=AIC(asragggAge_povpseudo)-AIC(asragggAge_povpseudoint)
	# get adjusted r^2 diffferences
	asrPDiffAdj[b]=summary(asrpgAge_povint)$adj.r.squared-summary(asrpgAge_pov)$adj.r.squared
	asrintDiffAdj[b]=summary(asrintgAge_povint)$adj.r.squared-summary(asrintgAge_pov)$adj.r.squared
	asrextDiffAdj[b]=summary(asrextgAge_povint)$adj.r.squared-summary(asrextgAge_pov)$adj.r.squared
	asrsomDiffAdj[b]=summary(asrsomgAge_povint)$adj.r.squared-summary(asrsomgAge_pov)$adj.r.squared
	asranxDiffAdj[b]=summary(asranxgAge_povint)$adj.r.squared-summary(asranxgAge_pov)$adj.r.squared
	asrthoDiffAdj[b]=summary(asrthogAge_povint)$adj.r.squared-summary(asrthogAge_pov)$adj.r.squared
	asrwitDiffAdj[b]=summary(asrwitgAge_povint)$adj.r.squared-summary(asrwitgAge_pov)$adj.r.squared
	asrattDiffAdj[b]=summary(asrattgAge_povint)$adj.r.squared-summary(asrattgAge_pov)$adj.r.squared
	asrrulDiffAdj[b]=summary(asrrulgAge_povint)$adj.r.squared-summary(asrrulgAge_pov)$adj.r.squared
	asraggDiffAdj[b]=summary(asragggAge_povint)$adj.r.squared-summary(asragggAge_pov)$adj.r.squared
	asrPDiffAdj_psuedo[b]=summary(asrpgAge_povpseudoint)$adj.r.squared-summary(asrpgAge_povpseudo)$adj.r.squared
	asrintDiffAdj_psuedo[b]=summary(asrintgAge_povpseudoint)$adj.r.squared-summary(asrintgAge_povpseudo)$adj.r.squared
	asrextDiffAdj_psuedo[b]=summary(asrextgAge_povpseudoint)$adj.r.squared-summary(asrextgAge_povpseudo)$adj.r.squared
	asrsomDiffAdj_psuedo[b]=summary(asrsomgAge_povpseudoint)$adj.r.squared-summary(asrsomgAge_povpseudo)$adj.r.squared
	asranxDiffAdj_psuedo[b]=summary(asranxgAge_povpseudoint)$adj.r.squared-summary(asranxgAge_povpseudo)$adj.r.squared
	asrthoDiffAdj_psuedo[b]=summary(asrthogAge_povpseudoint)$adj.r.squared-summary(asrthogAge_povpseudo)$adj.r.squared
	asrwitDiffAdj_psuedo[b]=summary(asrwitgAge_povpseudoint)$adj.r.squared-summary(asrwitgAge_povpseudo)$adj.r.squared
	asrattDiffAdj_psuedo[b]=summary(asrattgAge_povpseudoint)$adj.r.squared-summary(asrattgAge_povpseudo)$adj.r.squared
	asrrulDiffAdj_psuedo[b]=summary(asrrulgAge_povpseudoint)$adj.r.squared-summary(asrrulgAge_povpseudo)$adj.r.squared
	asraggDiffAdj_psuedo[b]=summary(asragggAge_povpseudoint)$adj.r.squared-summary(asragggAge_povpseudo)$adj.r.squared
	# make new predict dataframes with poverty variable
	predictDFasrp$poverty=1
	predictDFasrint$poverty=1
	predictDFasrext$poverty=1
	predictDFasrsom$poverty=1
	predictDFasranx$poverty=1
	predictDFasrtho$poverty=1
	predictDFasrwit$poverty=1
	predictDFasratt$poverty=1
	predictDFasrrul$poverty=1
	predictDFasragg$poverty=1
	# get poverty fits
	forFitasrP=predict(asrpgAge_povint,predictDFasrp)
	forFitasrInt=predict(asrintgAge_povint,predictDFasrint)
	forFitasrExt=predict(asrextgAge_povint,predictDFasrext)
	forFitasrSom=predict(asrsomgAge_povint,predictDFasrsom)
	forFitasrAnx=predict(asranxgAge_povint,predictDFasranx)
	forFitasrTho=predict(asrthogAge_povint,predictDFasrtho)
	forFitasrWit=predict(asrwitgAge_povint,predictDFasrwit)
	forFitasrAtt=predict(asrattgAge_povint,predictDFasratt)
	forFitasrRul=predict(asrrulgAge_povint,predictDFasrrul)
	forFitasrAgg=predict(asragggAge_povint,predictDFasragg)
	# print out fit
	asrPFitPov[b,1:(basrpmax+1)]=forFitasrP
	asrintFitPov[b,1:(basrintmax+1)]=forFitasrInt
	asrextFitPov[b,1:(basrextmax+1)]=forFitasrExt
	asrsomFitPov[b,1:(basrsommax+1)]=forFitasrSom
	asranxFitPov[b,1:(basranxmax+1)]=forFitasrAnx
	asrthoFitPov[b,1:(basrthomax+1)]=forFitasrTho
	asrwitFitPov[b,1:(basrwitmax+1)]=forFitasrWit
	asrattFitPov[b,1:(basrattmax+1)]=forFitasrAtt
	asrrulFitPov[b,1:(basrrulmax+1)]=forFitasrRul
	asraggFitPov[b,1:(basraggmax+1)]=forFitasrAgg
	# update predict df to nonpoverty, get nonpoverty fits
	predictDFasrp$poverty=0
	predictDFasrint$poverty=0
	predictDFasrext$poverty=0
	predictDFasrsom$poverty=0
	predictDFasranx$poverty=0
	predictDFasrtho$poverty=0
	predictDFasrwit$poverty=0
	predictDFasratt$poverty=0
	predictDFasrrul$poverty=0
	predictDFasragg$poverty=0
	# get nonpoverty fits
	forFitasrP=predict(asrpgAge_povint,predictDFasrp)
	forFitasrInt=predict(asrintgAge_povint,predictDFasrint)
	forFitasrExt=predict(asrextgAge_povint,predictDFasrext)
	forFitasrSom=predict(asrsomgAge_povint,predictDFasrsom)
	forFitasrAnx=predict(asranxgAge_povint,predictDFasranx)
	forFitasrTho=predict(asrthogAge_povint,predictDFasrtho)
	forFitasrWit=predict(asrwitgAge_povint,predictDFasrwit)
	forFitasrAtt=predict(asrattgAge_povint,predictDFasratt)
	forFitasrRul=predict(asrrulgAge_povint,predictDFasrrul)
	forFitasrAgg=predict(asragggAge_povint,predictDFasragg)
	# print out fit
	asrPFitNonPov[b,1:(basrpmax+1)]=forFitasrP
	asrintFitNonPov[b,1:(basrintmax+1)]=forFitasrInt
	asrextFitNonPov[b,1:(basrextmax+1)]=forFitasrExt
	asrsomFitNonPov[b,1:(basrsommax+1)]=forFitasrSom
	asranxFitNonPov[b,1:(basranxmax+1)]=forFitasrAnx
	asrthoFitNonPov[b,1:(basrthomax+1)]=forFitasrTho
	asrwitFitNonPov[b,1:(basrwitmax+1)]=forFitasrWit
	asrattFitNonPov[b,1:(basrattmax+1)]=forFitasrAtt
	asrrulFitNonPov[b,1:(basrrulmax+1)]=forFitasrRul
	asraggFitNonPov[b,1:(basraggmax+1)]=forFitasrAgg
}
# SAVEOUT
# save out all difference in AIC vectors, differences in adjusted R^2 vectors
outdf=data.frame(asrPDiff,asrintDiff,asrextDiff,asrsomDiff,asranxDiff,asrthoDiff,asrwitDiff,asrattDiff,asrrulDiff,asraggDiff,asrPDiffAdj,asrintDiffAdj,asrextDiffAdj,asrsomDiffAdj,asranxDiffAdj,asrthoDiffAdj,asrwitDiffAdj,asrattDiffAdj,asrrulDiffAdj,asraggDiffAdj)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_asr.rds')
# save out pseudo versions for comparison if needed
outdf=data.frame(asrPDiff_psuedo,asrintDiff_psuedo,asrextDiff_psuedo,asrsomDiff_psuedo,asranxDiff_psuedo,asrthoDiff_psuedo,asrwitDiff_psuedo,asrattDiff_psuedo,asrrulDiff_psuedo,asraggDiff_psuedo,asrPDiffAdj_psuedo,asrintDiffAdj_psuedo,asrextDiffAdj_psuedo,asrsomDiffAdj_psuedo,asranxDiffAdj_psuedo,asrthoDiffAdj_psuedo,asrwitDiffAdj_psuedo,asrattDiffAdj_psuedo,asrrulDiffAdj_psuedo,asraggDiffAdj_psuedo)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_asr_psuedo.rds')
# save out poverty and nonpoverty fits
outdf=data.frame(asrPFitPov,asrPFitNonPov,asrintFitPov,asrintFitNonPov,asrextFitPov,asrextFitNonPov,asrsomFitPov,asrsomFitNonPov,asranxFitPov,asranxFitNonPov,asrthoFitPov,asrthoFitNonPov,asrwitFitPov,asrwitFitNonPov,asrattFitPov,asrattFitNonPov,asrrulFitPov,asrrulFitNonPov,asraggFitPov,asraggFitNonPov)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_asr_pNp.rds')
# save out version with all cbcl and asr derivs
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,attDeriv,rulDeriv,aggDeriv,asrpDeriv,asrintDeriv,asrextDeriv,asrsomDeriv,asranxDeriv,asrthoDeriv,asrwitDeriv,asrattDeriv,asrrulDeriv,asraggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_asr.rds')
# save out version with all cbcl and asr fits
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,attFit,rulFit,aggFit,asrPFit,asrintFit,asrextFit,asrsomFit,asranxFit,asrthoFit,asrwitFit,asrattFit,asrrulFit,asraggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_asr.rds')
print('done with g~p fit bootstrapping!')
