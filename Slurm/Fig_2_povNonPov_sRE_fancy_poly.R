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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','site','g','subjectkey','interview_age','INR')]
# get length of df for later
lenDF=dim(masterdf)[1]
# set site to facotr
masterdf$site=as.factor(masterdf$site)
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
# initialize version of each for poverty and nonpoverty groups. Note we now need 5 different levels of poverty (pov0,pov1,pov2,pov3,pov4)
pFitPov0=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov0=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov0=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov0=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov0=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov0=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov0=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov0=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov0=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov0=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov0=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
pFitPov1=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov1=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov1=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov1=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov1=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov1=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov1=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov1=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov1=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov1=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov1=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
pFitPov2=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov2=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov2=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov2=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov2=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov2=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov2=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov2=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov2=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov2=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov2=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
pFitPov3=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov3=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov3=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov3=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov3=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov3=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov3=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov3=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov3=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov3=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov3=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
pFitPov4=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov4=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov4=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov4=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov4=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov4=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov4=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov4=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov4=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov4=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov4=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
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
# finally, add pov classification (v5, based on Kim HH et al 2022)
masterdf$poverty=0
masterdf$poverty[masterdf$INR>.5]=1
masterdf$poverty[masterdf$INR>1]=2
masterdf$poverty[masterdf$INR>2]=3
masterdf$poverty[masterdf$INR>4]=4
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
	# get unique instances of subjects, note duplcated and !duplicated are same because 2x instances of each subj
	unqSubjs=bootSamp[duplicated(bootSamp$subjectkey),]
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
	# set age to to median for predict df, also set child symptom score to median for predict df
	#####################################
        # set age to to median for predict df, +1 to max values because 0-max sequence is 1 > max in length.
        predictDFp=data.frame(eachPcount,rep(median(masterdf$interview_age),bpmax+1),rep('site13',bpmax+1))
        predictDFint=data.frame(eachIntcount,rep(median(masterdf$interview_age),bimax+1),rep('site13',bimax+1))
        predictDFext=data.frame(eachExtcount,rep(median(masterdf$interview_age),bemax+1),rep('site13',bemax+1))
        predictDFsom=data.frame(eachSomcount,rep(median(masterdf$interview_age),bsommax+1),rep('site13',bsommax+1))
        predictDFanx=data.frame(eachAnxcount,rep(median(masterdf$interview_age),banxmax+1),rep('site13',banxmax+1))
        predictDFtho=data.frame(eachThocount,rep(median(masterdf$interview_age),bthomax+1),rep('site13',bthomax+1))
        predictDFwit=data.frame(eachWitcount,rep(median(masterdf$interview_age),bwitmax+1),rep('site13',bwitmax+1))
        predictDFsoc=data.frame(eachSoccount,rep(median(masterdf$interview_age),bsocmax+1),rep('site13',bsocmax+1))
        predictDFatt=data.frame(eachAttcount,rep(median(masterdf$interview_age),battmax+1),rep('site13',battmax+1))
        predictDFrul=data.frame(eachRulcount,rep(median(masterdf$interview_age),brulmax+1),rep('site13',brulmax+1))
        predictDFagg=data.frame(eachAggcount,rep(median(masterdf$interview_age),baggmax+1),rep('site13',baggmax+1))
        # set colnames so predict can work
        colnames(predictDFp)=c('cbcl_scr_syn_totprob_r','interview_age','site')
        colnames(predictDFint)=c('cbcl_scr_syn_internal_r','interview_age','site')
        colnames(predictDFext)=c('cbcl_scr_syn_external_r','interview_age','site')
        colnames(predictDFsom)=c('cbcl_scr_syn_somatic_r','interview_age','site')
        colnames(predictDFanx)=c('cbcl_scr_syn_anxdep_r','interview_age','site')
        colnames(predictDFtho)=c('cbcl_scr_syn_thought_r','interview_age','site')
        colnames(predictDFwit)=c('cbcl_scr_syn_withdep_r','interview_age','site')
        colnames(predictDFsoc)=c('cbcl_scr_syn_social_r','interview_age','site')
        colnames(predictDFatt)=c('cbcl_scr_syn_attention_r','interview_age','site')
        colnames(predictDFrul)=c('cbcl_scr_syn_rulebreak_r','interview_age','site')
        colnames(predictDFagg)=c('cbcl_scr_syn_aggressive_r','interview_age','site')	
	####### II PREDICT POVERTY INTERACTIONS #######
	# fit models with standalone poverty term
	# cbcl
	pgAge_pov=bam(g~poly(cbcl_scr_syn_totprob_r,2)*poverty+s(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	intgAge_pov=bam(g~poly(cbcl_scr_syn_internal_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	extgAge_pov=bam(g~poly(cbcl_scr_syn_external_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	somgAge_pov=bam(g~poly(cbcl_scr_syn_somatic_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	anxgAge_pov=bam(g~poly(cbcl_scr_syn_anxdep_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	thogAge_pov=bam(g~poly(cbcl_scr_syn_thought_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	witgAge_pov=bam(g~poly(cbcl_scr_syn_withdep_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	socgAge_pov=bam(g~poly(cbcl_scr_syn_social_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	attgAge_pov=bam(g~poly(cbcl_scr_syn_attention_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	rulgAge_pov=bam(g~poly(cbcl_scr_syn_rulebreak_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	agggAge_pov=bam(g~poly(cbcl_scr_syn_aggressive_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	# make new predict dataframes with poverty variable
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
	# get poverty fits
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	pFitPov0[b,1:(bpmax+1)]=forFitp
	intFitPov0[b,1:(bimax+1)]=forFitint
	extFitPov0[b,1:(bemax+1)]=forFitext
	somFitPov0[b,1:(bsommax+1)]=forFitsom
	anxFitPov0[b,1:(banxmax+1)]=forFitanx
	thoFitPov0[b,1:(bthomax+1)]=forFittho
	witFitPov0[b,1:(bwitmax+1)]=forFitwit
	socFitPov0[b,1:(bsocmax+1)]=forFitsoc
	attFitPov0[b,1:(battmax+1)]=forFitatt
	rulFitPov0[b,1:(brulmax+1)]=forFitrul
	aggFitPov0[b,1:(baggmax+1)]=forFitagg
	# for poverty at level 1
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
	# get poverty fits
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	pFitPov1[b,1:(bpmax+1)]=forFitp
	intFitPov1[b,1:(bimax+1)]=forFitint
	extFitPov1[b,1:(bemax+1)]=forFitext
	somFitPov1[b,1:(bsommax+1)]=forFitsom
	anxFitPov1[b,1:(banxmax+1)]=forFitanx
	thoFitPov1[b,1:(bthomax+1)]=forFittho
	witFitPov1[b,1:(bwitmax+1)]=forFitwit
	socFitPov1[b,1:(bsocmax+1)]=forFitsoc
	attFitPov1[b,1:(battmax+1)]=forFitatt
	rulFitPov1[b,1:(brulmax+1)]=forFitrul
	aggFitPov1[b,1:(baggmax+1)]=forFitagg
	# for poverty at level 2
	predictDFp$poverty=2
	predictDFint$poverty=2
	predictDFext$poverty=2
	predictDFsom$poverty=2
	predictDFanx$poverty=2
	predictDFtho$poverty=2
	predictDFwit$poverty=2
	predictDFsoc$poverty=2
	predictDFatt$poverty=2
	predictDFrul$poverty=2
	predictDFagg$poverty=2
	# get poverty fits
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	pFitPov2[b,1:(bpmax+1)]=forFitp
	intFitPov2[b,1:(bimax+1)]=forFitint
	extFitPov2[b,1:(bemax+1)]=forFitext
	somFitPov2[b,1:(bsommax+1)]=forFitsom
	anxFitPov2[b,1:(banxmax+1)]=forFitanx
	thoFitPov2[b,1:(bthomax+1)]=forFittho
	witFitPov2[b,1:(bwitmax+1)]=forFitwit
	socFitPov2[b,1:(bsocmax+1)]=forFitsoc
	attFitPov2[b,1:(battmax+1)]=forFitatt
	rulFitPov2[b,1:(brulmax+1)]=forFitrul
	aggFitPov2[b,1:(baggmax+1)]=forFitagg
	# for poverty at level 3
	predictDFp$poverty=3
	predictDFint$poverty=3
	predictDFext$poverty=3
	predictDFsom$poverty=3
	predictDFanx$poverty=3
	predictDFtho$poverty=3
	predictDFwit$poverty=3
	predictDFsoc$poverty=3
	predictDFatt$poverty=3
	predictDFrul$poverty=3
	predictDFagg$poverty=3
	# get poverty fits
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	pFitPov3[b,1:(bpmax+1)]=forFitp
	intFitPov3[b,1:(bimax+1)]=forFitint
	extFitPov3[b,1:(bemax+1)]=forFitext
	somFitPov3[b,1:(bsommax+1)]=forFitsom
	anxFitPov3[b,1:(banxmax+1)]=forFitanx
	thoFitPov3[b,1:(bthomax+1)]=forFittho
	witFitPov3[b,1:(bwitmax+1)]=forFitwit
	socFitPov3[b,1:(bsocmax+1)]=forFitsoc
	attFitPov3[b,1:(battmax+1)]=forFitatt
	rulFitPov3[b,1:(brulmax+1)]=forFitrul
	aggFitPov3[b,1:(baggmax+1)]=forFitagg
	# for poverty at level 4
	predictDFp$poverty=4
	predictDFint$poverty=4
	predictDFext$poverty=4
	predictDFsom$poverty=4
	predictDFanx$poverty=4
	predictDFtho$poverty=4
	predictDFwit$poverty=4
	predictDFsoc$poverty=4
	predictDFatt$poverty=4
	predictDFrul$poverty=4
	predictDFagg$poverty=4
	# get poverty fits
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	pFitPov4[b,1:(bpmax+1)]=forFitp
	intFitPov4[b,1:(bimax+1)]=forFitint
	extFitPov4[b,1:(bemax+1)]=forFitext
	somFitPov4[b,1:(bsommax+1)]=forFitsom
	anxFitPov4[b,1:(banxmax+1)]=forFitanx
	thoFitPov4[b,1:(bthomax+1)]=forFittho
	witFitPov4[b,1:(bwitmax+1)]=forFitwit
	socFitPov4[b,1:(bsocmax+1)]=forFitsoc
	attFitPov4[b,1:(battmax+1)]=forFitatt
	rulFitPov4[b,1:(brulmax+1)]=forFitrul
	aggFitPov4[b,1:(baggmax+1)]=forFitagg
}
# SAVEOUT fits for pov 0
saveRDS(pFitPov0,intFitPov0,extFitPov0,somFitPov0,anxFitPov0,thoFitPov0,witFitPov0,socFitPov0,attFitPov0,rulFitPov0,aggFitPov0,file='/oak/stanford/groups/leanew1/users/apines/data/gp/FitsPov0_fancy_poly.rds')
# SAVEOUT fits for pov 1
saveRDS(pFitPov1,intFitPov1,extFitPov1,somFitPov1,anxFitPov1,thoFitPov1,witFitPov1,socFitPov1,attFitPov1,rulFitPov1,aggFitPov1,file='/oak/stanford/groups/leanew1/users/apines/data/gp/FitsPov1_fancy_poly.rds')
# SAVEOUT fits for pov 2
saveRDS(pFitPov2,intFitPov2,extFitPov2,somFitPov2,anxFitPov2,thoFitPov2,witFitPov2,socFitPov2,attFitPov2,rulFitPov2,aggFitPov2,file='/oak/stanford/groups/leanew1/users/apines/data/gp/FitsPov2_fancy_poly.rds')
# SAVEOUT fits for pov 3
saveRDS(pFitPov3,intFitPov3,extFitPov3,somFitPov3,anxFitPov3,thoFitPov3,witFitPov3,socFitPov3,attFitPov3,rulFitPov3,aggFitPov3,file='/oak/stanford/groups/leanew1/users/apines/data/gp/FitsPov3_fancy_poly.rds')
# SAVEOUT fits for pov 4
saveRDS(pFitPov4,intFitPov4,extFitPov4,somFitPov4,anxFitPov4,thoFitPov4,witFitPov4,socFitPov4,attFitPov4,rulFitPov4,aggFitPov4,file='/oak/stanford/groups/leanew1/users/apines/data/gp/FitsPov4_fancy_poly.rds')

print('done with g~p fit bootstrapping!')
