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
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
iMaxVal=max(masterdf$cbcl_scr_syn_internal_r)
emaxVal=max(masterdf$cbcl_scr_syn_external_r)
somMaxVal=max(masterdf$cbcl_scr_syn_somatic_r)
anxMaxVal=max(masterdf$cbcl_scr_syn_anxdep_r)
thoMaxVal=max(masterdf$cbcl_scr_syn_thought_r)
witMaxVal=max(masterdf$cbcl_scr_syn_withdep_r)
socMaxVal=max(masterdf$cbcl_scr_syn_social_r)
attMaxVal=max(masterdf$cbcl_scr_syn_attention_r)
rulMaxVal=max(masterdf$cbcl_scr_syn_rulebreak_r)
aggMaxVal=max(masterdf$cbcl_scr_syn_aggressive_r)
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
# predicted values: set to maximum value for ncol
pFit=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFit=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFit=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFit=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFit=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFit=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
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
socFitPov=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
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
socFitNonPov=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
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
# initialize difference in AIC for each cbcl model with and without poverty interaction
pDiff=rep(0,10000)
intDiff=rep(0,10000)
extDiff=rep(0,10000)
somDiff=rep(0,10000)
anxDiff=rep(0,10000)
thoDiff=rep(0,10000)
witDiff=rep(0,10000)
socDiff=rep(0,10000)
attDiff=rep(0,10000)
rulDiff=rep(0,10000)
aggDiff=rep(0,10000)
# difference in AIC for each cbcl model with psuedo poverty interaction
pDiffPseudo=rep(0,10000)
intDiffPseudo=rep(0,10000)
extDiffPseudo=rep(0,10000)
somDiffPseudo=rep(0,10000)
anxDiffPseudo=rep(0,10000)
thoDiffPseudo=rep(0,10000)
witDiffPseudo=rep(0,10000)
socDiffPseudo=rep(0,10000)
attDiffPseudo=rep(0,10000)
rulDiffPseudo=rep(0,10000)
aggDiffPseudo=rep(0,10000)
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
# maximum value in each iteration
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
somMax=rep(0,10000)
anxMax=rep(0,10000)
thoMax=rep(0,10000)
witMax=rep(0,10000)
socMax=rep(0,10000)
attMax=rep(0,10000)
rulMax=rep(0,10000)
aggMax=rep(0,10000)
# finally, add pov classification (v2)
masterdf$poverty=0
# poverty now defined in sample construction
masterdf$poverty[masterdf$Pov_v2==1]=1
masterdf$poverty=as.factor(masterdf$poverty)
# loop over manual bootstrap
for (b in 8001:10000){
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
	# get unique instances of subjects, note duplcated and !duplicated are same because 2x instances of each subj
	unqSubjs=bootSamp[duplicated(bootSamp$subjectkey),]
        # get count poverty in this boot
	subjects_with_poverty = unqSubjs$subjectkey[unqSubjs$poverty == 1]
	subjects_with_pseudopoverty = sample(unqSubjs$subjectkey, length(subjects_with_poverty))
	# create a vector to assign pseudopoverty to entire subjects
	bootSamp$pseudopoverty = 0
        # randomly assign Povcount people to pseudopoverty
	bootSamp$pseudopoverty[bootSamp$subjectkey %in% subjects_with_pseudopoverty] = 1
	bootSamp$pseudopoverty=as.factor(bootSamp$pseudopoverty)
        # make a df that is just the same number of kids in poverty, but actually from the non-poverty group. The point is to see if we can recover g~p slope in a an equivalent number of non-poverty kids
        #bootSamp$pseudopoverty2=0
        #bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
	##############################
	# get max for this iteration
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
	bimax=max(bootSamp$cbcl_scr_syn_internal_r)
	bemax=max(bootSamp$cbcl_scr_syn_external_r)
	bsommax=max(bootSamp$cbcl_scr_syn_somatic_r)
	banxmax=max(bootSamp$cbcl_scr_syn_anxdep_r)
	bthomax=max(bootSamp$cbcl_scr_syn_thought_r)
	bwitmax=max(bootSamp$cbcl_scr_syn_withdep_r)
	bsocmax=max(bootSamp$cbcl_scr_syn_social_r)
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
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4),data=bootSamp)
	somgAge<-bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4),data=bootSamp)
	anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4),data=bootSamp)
	thogAge<-bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4),data=bootSamp)
	witgAge<-bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4),data=bootSamp)
	socgAge<-bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4),data=bootSamp)
	attgAge<-bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4),data=bootSamp)
	rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4),data=bootSamp)
	agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4),data=bootSamp)
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
	eachSoccount=seq(0:bsocmax)
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
	predictDFp=data.frame(eachPcount,rep(median(masterdf$interview_age),(bpmax+1)))
	predictDFint=data.frame(eachIntcount,rep(median(masterdf$interview_age),(bimax+1)))
	predictDFext=data.frame(eachExtcount,rep(median(masterdf$interview_age),(bemax+1)))
	predictDFsom=data.frame(eachSomcount,rep(median(masterdf$interview_age),(bsommax+1)))
	predictDFanx=data.frame(eachAnxcount,rep(median(masterdf$interview_age),(banxmax+1)))
	predictDFtho=data.frame(eachThocount,rep(median(masterdf$interview_age),(bthomax+1)))
	predictDFwit=data.frame(eachWitcount,rep(median(masterdf$interview_age),(bwitmax+1)))
	predictDFsoc=data.frame(eachSoccount,rep(median(masterdf$interview_age),(bsocmax+1)))
	predictDFatt=data.frame(eachAttcount,rep(median(masterdf$interview_age),(battmax+1)))
	predictDFrul=data.frame(eachRulcount,rep(median(masterdf$interview_age),(brulmax+1)))
	predictDFagg=data.frame(eachAggcount,rep(median(masterdf$interview_age),(baggmax+1)))
	predictDFasrp=data.frame(eachasrPcount,rep(median(masterdf$interview_age),(basrpmax+1)))
	predictDFasrint=data.frame(eachasrIntcount,rep(median(masterdf$interview_age),(basrintmax+1)))
	predictDFasrext=data.frame(eachasrExtcount,rep(median(masterdf$interview_age),(basrextmax+1)))
	predictDFasrsom=data.frame(eachasrSomcount,rep(median(masterdf$interview_age),(basrsommax+1)))
	predictDFasranx=data.frame(eachasrAnxcount,rep(median(masterdf$interview_age),(basranxmax+1)))
	predictDFasrtho=data.frame(eachasrThocount,rep(median(masterdf$interview_age),(basrthomax+1)))
	predictDFasrwit=data.frame(eachasrWitcount,rep(median(masterdf$interview_age),(basrwitmax+1)))
	predictDFasratt=data.frame(eachasrAttcount,rep(median(masterdf$interview_age),(basrattmax+1)))
	predictDFasrrul=data.frame(eachasrRulcount,rep(median(masterdf$interview_age),(basrrulmax+1)))
	predictDFasragg=data.frame(eachasrAggcount,rep(median(masterdf$interview_age),(basraggmax+1)))
	# set colnames so predict can work
	colnames(predictDFp)=c('cbcl_scr_syn_totprob_r','interview_age')
	colnames(predictDFint)=c('cbcl_scr_syn_internal_r','interview_age')
	colnames(predictDFext)=c('cbcl_scr_syn_external_r','interview_age')
	colnames(predictDFsom)=c('cbcl_scr_syn_somatic_r','interview_age')
	colnames(predictDFanx)=c('cbcl_scr_syn_anxdep_r','interview_age')
	colnames(predictDFtho)=c('cbcl_scr_syn_thought_r','interview_age')
	colnames(predictDFwit)=c('cbcl_scr_syn_withdep_r','interview_age')
	colnames(predictDFsoc)=c('cbcl_scr_syn_social_r','interview_age')
	colnames(predictDFatt)=c('cbcl_scr_syn_attention_r','interview_age')
	colnames(predictDFrul)=c('cbcl_scr_syn_rulebreak_r','interview_age')
	colnames(predictDFagg)=c('cbcl_scr_syn_aggressive_r','interview_age')
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
	forFitSoc=predict(socgAge,predictDFsoc)
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
	pFit[b,1:(bpmax+1)]=forFitP
	intFit[b,1:(bimax+1)]=forFitInt
	extFit[b,1:(bemax+1)]=forFitExt
	somFit[b,1:(bsommax+1)]=forFitSom
	anxFit[b,1:(banxmax+1)]=forFitAnx
	thoFit[b,1:(bthomax+1)]=forFitTho
	witFit[b,1:(bwitmax+1)]=forFitWit
	socFit[b,1:(bsocmax+1)]=forFitSoc
	attFit[b,1:(battmax+1)]=forFitAtt
	rulFit[b,1:(brulmax+1)]=forFitRul
	aggFit[b,1:(baggmax+1)]=forFitAgg
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
	####### II PREDICT POVERTY INTERACTIONS #######
	# fit models with standalone poverty term
	# cbcl
	pgAge_pov=bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	intgAge_pov=bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	extgAge_pov=bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	somgAge_pov=bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	anxgAge_pov=bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	thogAge_pov=bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	witgAge_pov=bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	socgAge_pov=bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	attgAge_pov=bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	rulgAge_pov=bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	agggAge_pov=bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	# fit versions with poverty interactions on symptoms
	pgAge_povint=bam(g~s(cbcl_scr_syn_totprob_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	intgAge_povint=bam(g~s(cbcl_scr_syn_internal_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	extgAge_povint=bam(g~s(cbcl_scr_syn_external_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	somgAge_povint=bam(g~s(cbcl_scr_syn_somatic_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	anxgAge_povint=bam(g~s(cbcl_scr_syn_anxdep_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	thogAge_povint=bam(g~s(cbcl_scr_syn_thought_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	witgAge_povint=bam(g~s(cbcl_scr_syn_withdep_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	socgAge_povint=bam(g~s(cbcl_scr_syn_social_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	attgAge_povint=bam(g~s(cbcl_scr_syn_attention_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	rulgAge_povint=bam(g~s(cbcl_scr_syn_rulebreak_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	agggAge_povint=bam(g~s(cbcl_scr_syn_aggressive_r,by=poverty,k=4)+s(interview_age,k=4)+poverty,data=bootSamp)
	# fit version on pseudopoverty variable
	pgAge_povpseudo=bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	intgAge_povpseudo=bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	extgAge_povpseudo=bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	somgAge_povpseudo=bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	anxgAge_povpseudo=bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	thogAge_povpseudo=bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	witgAge_povpseudo=bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	socgAge_povpseudo=bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	attgAge_povpseudo=bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	rulgAge_povpseudo=bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	agggAge_povpseudo=bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	# fit version on pseudopoverty variable with interactions
	pgAge_povpseudoint=bam(g~s(cbcl_scr_syn_totprob_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	intgAge_povpseudoint=bam(g~s(cbcl_scr_syn_internal_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	extgAge_povpseudoint=bam(g~s(cbcl_scr_syn_external_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	somgAge_povpseudoint=bam(g~s(cbcl_scr_syn_somatic_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	anxgAge_povpseudoint=bam(g~s(cbcl_scr_syn_anxdep_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	thogAge_povpseudoint=bam(g~s(cbcl_scr_syn_thought_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	witgAge_povpseudoint=bam(g~s(cbcl_scr_syn_withdep_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	socgAge_povpseudoint=bam(g~s(cbcl_scr_syn_social_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	attgAge_povpseudoint=bam(g~s(cbcl_scr_syn_attention_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	rulgAge_povpseudoint=bam(g~s(cbcl_scr_syn_rulebreak_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	agggAge_povpseudoint=bam(g~s(cbcl_scr_syn_aggressive_r,by=pseudopoverty,k=4)+s(interview_age,k=4)+pseudopoverty,data=bootSamp)
	# asr
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
	# cbcl
	pDiff[b]=AIC(pgAge_pov)-AIC(pgAge_povint)
	intDiff[b]=AIC(intgAge_pov)-AIC(intgAge_povint)
	extDiff[b]=AIC(extgAge_pov)-AIC(extgAge_povint)
	somDiff[b]=AIC(somgAge_pov)-AIC(somgAge_povint)
	anxDiff[b]=AIC(anxgAge_pov)-AIC(anxgAge_povint)
	thoDiff[b]=AIC(thogAge_pov)-AIC(thogAge_povint)
	witDiff[b]=AIC(witgAge_pov)-AIC(witgAge_povint)
	socDiff[b]=AIC(socgAge_pov)-AIC(socgAge_povint)
	attDiff[b]=AIC(attgAge_pov)-AIC(attgAge_povint)
	rulDiff[b]=AIC(rulgAge_pov)-AIC(rulgAge_povint)
	aggDiff[b]=AIC(agggAge_pov)-AIC(agggAge_povint)
	pDiffPseudo[b]=AIC(pgAge_povpseudo)-AIC(pgAge_povpseudoint)
	intDiffPseudo[b]=AIC(intgAge_povpseudo)-AIC(intgAge_povpseudoint)
	extDiffPseudo[b]=AIC(extgAge_povpseudo)-AIC(extgAge_povpseudoint)
	somDiffPseudo[b]=AIC(somgAge_povpseudo)-AIC(somgAge_povpseudoint)
	anxDiffPseudo[b]=AIC(anxgAge_povpseudo)-AIC(anxgAge_povpseudoint)
	thoDiffPseudo[b]=AIC(thogAge_povpseudo)-AIC(thogAge_povpseudoint)
	witDiffPseudo[b]=AIC(witgAge_povpseudo)-AIC(witgAge_povpseudoint)
	socDiffPseudo[b]=AIC(socgAge_povpseudo)-AIC(socgAge_povpseudoint)
	attDiffPseudo[b]=AIC(attgAge_povpseudo)-AIC(attgAge_povpseudoint)
	rulDiffPseudo[b]=AIC(rulgAge_povpseudo)-AIC(rulgAge_povpseudoint)
	aggDiffPseudo[b]=AIC(agggAge_povpseudo)-AIC(agggAge_povpseudoint)
	# asr
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
	asrPDiffPseudo[b]=AIC(asrpgAge_povpseudo)-AIC(asrpgAge_povpseudoint)
	asrintDiffPseudo[b]=AIC(asrintgAge_povpseudo)-AIC(asrintgAge_povpseudoint)
	asrextDiffPseudo[b]=AIC(asrextgAge_povpseudo)-AIC(asrextgAge_povpseudoint)
	asrsomDiffPseudo[b]=AIC(asrsomgAge_povpseudo)-AIC(asrsomgAge_povpseudoint)
	asranxDiffPseudo[b]=AIC(asranxgAge_povpseudo)-AIC(asranxgAge_povpseudoint)
	asrthoDiffPseudo[b]=AIC(asrthogAge_povpseudo)-AIC(asrthogAge_povpseudoint)
	asrwitDiffPseudo[b]=AIC(asrwitgAge_povpseudo)-AIC(asrwitgAge_povpseudoint)
	asrattDiffPseudo[b]=AIC(asrattgAge_povpseudo)-AIC(asrattgAge_povpseudoint)
	asrrulDiffPseudo[b]=AIC(asrrulgAge_povpseudo)-AIC(asrrulgAge_povpseudoint)
	asraggDiffPseudo[b]=AIC(asragggAge_povpseudo)-AIC(asragggAge_povpseudoint)
	# make new predict dataframes with poverty variable
	# cbcl
	predictDFp$poverty=1
	predictDFint$poverty=1
	predictDFext$poverty=1
	predictDFsom$poverty=1
	predictDFanx$poverty=1
	predictDFtho$poverty=1	
	predictDFwit$poverty=1
	predictDFsoc$poverty=1
	predictDFatt$poverty=1
	predictDFrul$poverty=1
	predictDFagg$poverty=1
	# asr
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
	# cbcl
	forFitp=predict(pgAge_povint,predictDFp)
	forFitint=predict(intgAge_povint,predictDFint)
	forFitext=predict(extgAge_povint,predictDFext)
	forFitsom=predict(somgAge_povint,predictDFsom)
	forFitanx=predict(anxgAge_povint,predictDFanx)
	forFittho=predict(thogAge_povint,predictDFtho)
	forFitwit=predict(witgAge_povint,predictDFwit)
	forFitsoc=predict(socgAge_povint,predictDFsoc)
	forFitatt=predict(attgAge_povint,predictDFatt)
	forFitrul=predict(rulgAge_povint,predictDFrul)
	forFitagg=predict(agggAge_povint,predictDFagg)
	# asr
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
	# cbcl
	pFitPov[b,1:(bpmax+1)]=forFitp
	intFitPov[b,1:(bimax+1)]=forFitint
	extFitPov[b,1:(bemax+1)]=forFitext
	somFitPov[b,1:(bsommax+1)]=forFitsom
	anxFitPov[b,1:(banxmax+1)]=forFitanx
	thoFitPov[b,1:(bthomax+1)]=forFittho
	witFitPov[b,1:(bwitmax+1)]=forFitwit
	socFitPov[b,1:(bsocmax+1)]=forFitsoc
	attFitPov[b,1:(battmax+1)]=forFitatt
	rulFitPov[b,1:(brulmax+1)]=forFitrul
	aggFitPov[b,1:(baggmax+1)]=forFitagg
	# asr
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
	# cbcl
	predictDFp$poverty=0
	predictDFint$poverty=0
	predictDFext$poverty=0
	predictDFsom$poverty=0
	predictDFanx$poverty=0
	predictDFtho$poverty=0
	predictDFwit$poverty=0
	predictDFsoc$poverty=0
	predictDFatt$poverty=0
	predictDFrul$poverty=0
	predictDFagg$poverty=0
	# asr
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
	# cbcl
	forFitp=predict(pgAge_povint,predictDFp)
	forFitint=predict(intgAge_povint,predictDFint)
	forFitext=predict(extgAge_povint,predictDFext)
	forFitsom=predict(somgAge_povint,predictDFsom)
	forFitanx=predict(anxgAge_povint,predictDFanx)
	forFittho=predict(thogAge_povint,predictDFtho)
	forFitwit=predict(witgAge_povint,predictDFwit)
	forFitsoc=predict(socgAge_povint,predictDFsoc)
	forFitatt=predict(attgAge_povint,predictDFatt)
	forFitrul=predict(rulgAge_povint,predictDFrul)
	forFitagg=predict(agggAge_povint,predictDFagg)
	# asr
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
	# cbcl
	pFitNonPov[b,1:(bpmax+1)]=forFitp
	intFitNonPov[b,1:(bimax+1)]=forFitint
	extFitNonPov[b,1:(bemax+1)]=forFitext
	somFitNonPov[b,1:(bsommax+1)]=forFitsom
	anxFitNonPov[b,1:(banxmax+1)]=forFitanx
	thoFitNonPov[b,1:(bthomax+1)]=forFittho
	witFitNonPov[b,1:(bwitmax+1)]=forFitwit
	socFitNonPov[b,1:(bsocmax+1)]=forFitsoc
	attFitNonPov[b,1:(battmax+1)]=forFitatt
	rulFitNonPov[b,1:(brulmax+1)]=forFitrul
	aggFitNonPov[b,1:(baggmax+1)]=forFitagg
	# asr
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
# cbcl
outdf=data.frame(pDiff,intDiff,extDiff,somDiff,anxDiff,thoDiff,witDiff,socDiff,attDiff,rulDiff,aggDiff)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_cbcl5.rds')
# asr
outdf=data.frame(asrPDiff,asrintDiff,asrextDiff,asrsomDiff,asranxDiff,asrthoDiff,asrwitDiff,asrattDiff,asrrulDiff,asraggDiff)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_asr5.rds')

# save out pseudo versions for comparison if needed
# cbcl
outdf=data.frame(pDiffPseudo,intDiffPseudo,extDiffPseudo,somDiffPseudo,anxDiffPseudo,thoDiffPseudo,witDiffPseudo,socDiffPseudo,attDiffPseudo,rulDiffPseudo,aggDiffPseudo)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_cbclPseudo5.rds')
# asr
outdf=data.frame(asrPDiffPseudo,asrintDiffPseudo,asrextDiffPseudo,asrsomDiffPseudo,asranxDiffPseudo,asrthoDiffPseudo,asrwitDiffPseudo,asrattDiffPseudo,asrrulDiffPseudo,asraggDiffPseudo)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDiffBoots_asrPseudo5.rds')

# save out poverty and nonpoverty fits - cbcl
outdf=data.frame(pFitPov,pFitNonPov,intFitPov,intFitNonPov,extFitPov,extFitNonPov,somFitPov,somFitNonPov,anxFitPov,anxFitNonPov,thoFitPov,thoFitNonPov,witFitPov,witFitNonPov,socFitPov,socFitNonPov,attFitPov,attFitNonPov,rulFitPov,rulFitNonPov,aggFitPov,aggFitNonPov)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_cbcl_pNp5.rds')
# save out poverty and nonpoverty fits - asr
outdf=data.frame(asrPFitPov,asrPFitNonPov,asrintFitPov,asrintFitNonPov,asrextFitPov,asrextFitNonPov,asrsomFitPov,asrsomFitNonPov,asranxFitPov,asranxFitNonPov,asrthoFitPov,asrthoFitNonPov,asrwitFitPov,asrwitFitNonPov,asrattFitPov,asrattFitNonPov,asrrulFitPov,asrrulFitNonPov,asraggFitPov,asraggFitNonPov)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_asr_pNp5.rds')

# save out version with all cbcl and asr fits
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit,asrPFit,asrintFit,asrextFit,asrsomFit,asranxFit,asrthoFit,asrwitFit,asrattFit,asrrulFit,asraggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_cbclasr5.rds')
print('done with g~p fit bootstrapping!')
