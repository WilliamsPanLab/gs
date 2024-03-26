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
        bootSamp$pseudopoverty2=0
	Povcount=length(bootSamp$subjectkey[bootSamp$poverty==1])
	print(Povcount)
        bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
	# make bootsamp just those in pseudopov2 group
	bootSamp=bootSamp[bootSamp$pseudopoverty2==1,]
	print(dim(bootSamp))
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
}
# SAVEOUT
# save out version with all cbcl and asr fits
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit,asrPFit,asrintFit,asrextFit,asrsomFit,asranxFit,asrthoFit,asrwitFit,asrattFit,asrrulFit,asraggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_nMatched_psuedoPov_cbclasr1.rds')
print('done with g~p fit bootstrapping!')
