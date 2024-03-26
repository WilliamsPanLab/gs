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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age','site','income','parental_education')]
# set site to facotr
masterdf$site=as.factor(masterdf$site)
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
### initialize cross-boot vectors
# predicted derivatives: set to maximum value for ncol +1 because 0-max is 1 more than max in length
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
pDeriv=matrix(0,nrow=10000,ncol=pMaxVal+1)
intDeriv=matrix(0,nrow=10000,ncol=iMaxVal+1)
extDeriv=matrix(0,nrow=10000,ncol=emaxVal+1)
somDeriv=matrix(0,nrow=10000,ncol=somMaxVal+1)
anxDeriv=matrix(0,nrow=10000,ncol=anxMaxVal+1)
thoDeriv=matrix(0,nrow=10000,ncol=thoMaxVal+1)
witDeriv=matrix(0,nrow=10000,ncol=witMaxVal+1)
socDeriv=matrix(0,nrow=10000,ncol=socMaxVal+1)
attDeriv=matrix(0,nrow=10000,ncol=attMaxVal+1)
rulDeriv=matrix(0,nrow=10000,ncol=rulMaxVal+1)
aggDeriv=matrix(0,nrow=10000,ncol=aggMaxVal+1)
# predicted values: set to maximum value for ncol +1 because 0-max is 1 more than max in length
pFit=matrix(0,nrow=10000,ncol=pMaxVal+1)
intFit=matrix(0,nrow=10000,ncol=iMaxVal+1)
extFit=matrix(0,nrow=10000,ncol=emaxVal+1)
somFit=matrix(0,nrow=10000,ncol=somMaxVal+1)
anxFit=matrix(0,nrow=10000,ncol=anxMaxVal+1)
thoFit=matrix(0,nrow=10000,ncol=thoMaxVal+1)
witFit=matrix(0,nrow=10000,ncol=witMaxVal+1)
socFit=matrix(0,nrow=10000,ncol=socMaxVal+1)
attFit=matrix(0,nrow=10000,ncol=attMaxVal+1)
rulFit=matrix(0,nrow=10000,ncol=rulMaxVal+1)
aggFit=matrix(0,nrow=10000,ncol=aggMaxVal+1)

# for parental education instead of income in-model
edu_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal+1)
edu_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal+1)
edu_extDeriv=matrix(0,nrow=10000,ncol=emaxVal+1)
edu_somDeriv=matrix(0,nrow=10000,ncol=somMaxVal+1)
edu_anxDeriv=matrix(0,nrow=10000,ncol=anxMaxVal+1)
edu_thoDeriv=matrix(0,nrow=10000,ncol=thoMaxVal+1)
edu_witDeriv=matrix(0,nrow=10000,ncol=witMaxVal+1)
edu_socDeriv=matrix(0,nrow=10000,ncol=socMaxVal+1)
edu_attDeriv=matrix(0,nrow=10000,ncol=attMaxVal+1)
edu_rulDeriv=matrix(0,nrow=10000,ncol=rulMaxVal+1)
edu_aggDeriv=matrix(0,nrow=10000,ncol=aggMaxVal+1)
# predicted values: set to maximum value for ncol +1 because 0-max is 1 more than max in length
edu_pFit=matrix(0,nrow=10000,ncol=pMaxVal+1)
edu_intFit=matrix(0,nrow=10000,ncol=iMaxVal+1)
edu_extFit=matrix(0,nrow=10000,ncol=emaxVal+1)
edu_somFit=matrix(0,nrow=10000,ncol=somMaxVal+1)
edu_anxFit=matrix(0,nrow=10000,ncol=anxMaxVal+1)
edu_thoFit=matrix(0,nrow=10000,ncol=thoMaxVal+1)
edu_witFit=matrix(0,nrow=10000,ncol=witMaxVal+1)
edu_socFit=matrix(0,nrow=10000,ncol=socMaxVal+1)
edu_attFit=matrix(0,nrow=10000,ncol=attMaxVal+1)
edu_rulFit=matrix(0,nrow=10000,ncol=rulMaxVal+1)
edu_aggFit=matrix(0,nrow=10000,ncol=aggMaxVal+1)

#### for both in-model
both_pDeriv=matrix(0,nrow=10000,ncol=pMaxVal+1)
both_intDeriv=matrix(0,nrow=10000,ncol=iMaxVal+1)
both_extDeriv=matrix(0,nrow=10000,ncol=emaxVal+1)
both_somDeriv=matrix(0,nrow=10000,ncol=somMaxVal+1)
both_anxDeriv=matrix(0,nrow=10000,ncol=anxMaxVal+1)
both_thoDeriv=matrix(0,nrow=10000,ncol=thoMaxVal+1)
both_witDeriv=matrix(0,nrow=10000,ncol=witMaxVal+1)
both_socDeriv=matrix(0,nrow=10000,ncol=socMaxVal+1)
both_attDeriv=matrix(0,nrow=10000,ncol=attMaxVal+1)
both_rulDeriv=matrix(0,nrow=10000,ncol=rulMaxVal+1)
both_aggDeriv=matrix(0,nrow=10000,ncol=aggMaxVal+1)
# predicted values: set to maximum value for ncol +1 because 0-max is 1 more than max in length
both_pFit=matrix(0,nrow=10000,ncol=pMaxVal+1)
both_intFit=matrix(0,nrow=10000,ncol=iMaxVal+1)
both_extFit=matrix(0,nrow=10000,ncol=emaxVal+1)
both_somFit=matrix(0,nrow=10000,ncol=somMaxVal+1)
both_anxFit=matrix(0,nrow=10000,ncol=anxMaxVal+1)
both_thoFit=matrix(0,nrow=10000,ncol=thoMaxVal+1)
both_witFit=matrix(0,nrow=10000,ncol=witMaxVal+1)
both_socFit=matrix(0,nrow=10000,ncol=socMaxVal+1)
both_attFit=matrix(0,nrow=10000,ncol=attMaxVal+1)
both_rulFit=matrix(0,nrow=10000,ncol=rulMaxVal+1)
both_aggFit=matrix(0,nrow=10000,ncol=aggMaxVal+1)

# make an adjusted R^2 vector for each subscale, for income, education, and both
pAdjR2=rep(0,10000)
intAdjR2=rep(0,10000)
extAdjR2=rep(0,10000)
somAdjR2=rep(0,10000)
anxAdjR2=rep(0,10000)
thoAdjR2=rep(0,10000)
witAdjR2=rep(0,10000)
socAdjR2=rep(0,10000)
attAdjR2=rep(0,10000)
rulAdjR2=rep(0,10000)
aggAdjR2=rep(0,10000)
edu_pAdjR2=rep(0,10000)
edu_intAdjR2=rep(0,10000)
edu_extAdjR2=rep(0,10000)
edu_somAdjR2=rep(0,10000)
edu_anxAdjR2=rep(0,10000)
edu_thoAdjR2=rep(0,10000)
edu_witAdjR2=rep(0,10000)
edu_socAdjR2=rep(0,10000)
edu_attAdjR2=rep(0,10000)
edu_rulAdjR2=rep(0,10000)
edu_aggAdjR2=rep(0,10000)
both_pAdjR2=rep(0,10000)
both_intAdjR2=rep(0,10000)
both_extAdjR2=rep(0,10000)
both_somAdjR2=rep(0,10000)
both_anxAdjR2=rep(0,10000)
both_thoAdjR2=rep(0,10000)
both_witAdjR2=rep(0,10000)
both_socAdjR2=rep(0,10000)
both_attAdjR2=rep(0,10000)
both_rulAdjR2=rep(0,10000)
both_aggAdjR2=rep(0,10000)

# maximum value within-bootstrap
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

### get medians for predict df
# get median income
medIncome=median(masterdf$income)
# get median education 
medEdu=median(as.numeric(masterdf$parental_education))
# set to integer for equivalent fitting
masterdf$parental_education=as.numeric(masterdf$parental_education)

# loop over manual bootstrap
for (b in 1:2){
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
	bsocmax=max(bootSamp$cbcl_scr_syn_social_r)
	battmax=max(bootSamp$cbcl_scr_syn_attention_r)
	brulmax=max(bootSamp$cbcl_scr_syn_rulebreak_r)
	baggmax=max(bootSamp$cbcl_scr_syn_aggressive_r)
	######## I PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable, include income
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	somgAge<-bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	thogAge<-bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	witgAge<-bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	socgAge<-bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	attgAge<-bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4),data=bootSamp)
	## include edu ####
	edu_pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_intgAge<-bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_extgAge<-bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_somgAge<-bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_thogAge<-bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_witgAge<-bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_socgAge<-bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_attgAge<-bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	edu_agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(parental_education,k=4),data=bootSamp)
	# include both income and edu
	both_pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_intgAge<-bam(g~s(cbcl_scr_syn_internal_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_extgAge<-bam(g~s(cbcl_scr_syn_external_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_somgAge<-bam(g~s(cbcl_scr_syn_somatic_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_anxgAge<-bam(g~s(cbcl_scr_syn_anxdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_thogAge<-bam(g~s(cbcl_scr_syn_thought_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_witgAge<-bam(g~s(cbcl_scr_syn_withdep_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_socgAge<-bam(g~s(cbcl_scr_syn_social_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_attgAge<-bam(g~s(cbcl_scr_syn_attention_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_rulgAge<-bam(g~s(cbcl_scr_syn_rulebreak_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
	both_agggAge<-bam(g~s(cbcl_scr_syn_aggressive_r,k=4)+s(interview_age,k=4)+s(site,bs="re")+s(income,k=4)+s(parental_education,k=4),data=bootSamp)
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
	### add median income and edu to predict df ####
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
	# add median income, predict for each symptom count
	predictDFp$income=rep(medIncome,bpmax+1)
	predictDFint$income=rep(medIncome,bimax+1)
	predictDFext$income=rep(medIncome,bemax+1)
	predictDFsom$income=rep(medIncome,bsommax+1)
	predictDFanx$income=rep(medIncome,banxmax+1)
	predictDFtho$income=rep(medIncome,bthomax+1)
	predictDFwit$income=rep(medIncome,bwitmax+1)
	predictDFsoc$income=rep(medIncome,bsocmax+1)
	predictDFatt$income=rep(medIncome,battmax+1)
	predictDFrul$income=rep(medIncome,brulmax+1)
	predictDFagg$income=rep(medIncome,baggmax+1)
	# add median edu, predict for each symptom count
	predictDFp$parental_education=rep(medEdu,bpmax+1)
	predictDFint$parental_education=rep(medEdu,bimax+1)
	predictDFext$parental_education=rep(medEdu,bemax+1)
	predictDFsom$parental_education=rep(medEdu,bsommax+1)
	predictDFanx$parental_education=rep(medEdu,banxmax+1)
	predictDFtho$parental_education=rep(medEdu,bthomax+1)
	predictDFwit$parental_education=rep(medEdu,bwitmax+1)
	predictDFsoc$parental_education=rep(medEdu,bsocmax+1)
	predictDFatt$parental_education=rep(medEdu,battmax+1)
	predictDFrul$parental_education=rep(medEdu,brulmax+1)
	predictDFagg$parental_education=rep(medEdu,baggmax+1)
	# predict: income v4ersion
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
	# predict: edu version
	forFitP=predict(edu_pgAge,predictDFp)
	forFitInt=predict(edu_intgAge,predictDFint)
	forFitExt=predict(edu_extgAge,predictDFext)
	forFitSom=predict(edu_somgAge,predictDFsom)
	forFitAnx=predict(edu_anxgAge,predictDFanx)
	forFitTho=predict(edu_thogAge,predictDFtho)
	forFitWit=predict(edu_witgAge,predictDFwit)
	forFitSoc=predict(edu_socgAge,predictDFsoc)
	forFitAtt=predict(edu_attgAge,predictDFatt)
	forFitRul=predict(edu_rulgAge,predictDFrul)
	forFitAgg=predict(edu_agggAge,predictDFagg)
	# print out fit
	edu_pFit[b,1:(bpmax+1)]=forFitP
	edu_intFit[b,1:(bimax+1)]=forFitInt
	edu_extFit[b,1:(bemax+1)]=forFitExt
	edu_somFit[b,1:(bsommax+1)]=forFitSom
	edu_anxFit[b,1:(banxmax+1)]=forFitAnx
	edu_thoFit[b,1:(bthomax+1)]=forFitTho
	edu_witFit[b,1:(bwitmax+1)]=forFitWit
	edu_socFit[b,1:(bsocmax+1)]=forFitSoc
	edu_attFit[b,1:(battmax+1)]=forFitAtt
	edu_rulFit[b,1:(brulmax+1)]=forFitRul
	edu_aggFit[b,1:(baggmax+1)]=forFitAgg
	# predict: both version
	forFitP=predict(both_pgAge,predictDFp)
	forFitInt=predict(both_intgAge,predictDFint)
	forFitExt=predict(both_extgAge,predictDFext)
	forFitSom=predict(both_somgAge,predictDFsom)
	forFitAnx=predict(both_anxgAge,predictDFanx)
	forFitTho=predict(both_thogAge,predictDFtho)
	forFitWit=predict(both_witgAge,predictDFwit)
	forFitSoc=predict(both_socgAge,predictDFsoc)
	forFitAtt=predict(both_attgAge,predictDFatt)
	forFitRul=predict(both_rulgAge,predictDFrul)
	forFitAgg=predict(both_agggAge,predictDFagg)
	# print out fit
	both_pFit[b,1:(bpmax+1)]=forFitP
	both_intFit[b,1:(bimax+1)]=forFitInt
	both_extFit[b,1:(bemax+1)]=forFitExt
	both_somFit[b,1:(bsommax+1)]=forFitSom
	both_anxFit[b,1:(banxmax+1)]=forFitAnx
	both_thoFit[b,1:(bthomax+1)]=forFitTho
	both_witFit[b,1:(bwitmax+1)]=forFitWit
	both_socFit[b,1:(bsocmax+1)]=forFitSoc
	both_attFit[b,1:(battmax+1)]=forFitAtt
	both_rulFit[b,1:(brulmax+1)]=forFitRul
	both_aggFit[b,1:(baggmax+1)]=forFitAgg
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE,n=bpmax+1)
	forSplineint=derivatives(intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE,n=bimax+1)
	forSplineext=derivatives(extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE,n=bemax+1)
	forSplinesom=derivatives(somgAge,term='s(cbcl_scr_syn_somatic_r)',partial_match = TRUE,n=bsommax+1)
	forSplineanx=derivatives(anxgAge,term='s(cbcl_scr_syn_anxdep_r)',partial_match = TRUE,n=banxmax+1)
	forSplinetho=derivatives(thogAge,term='s(cbcl_scr_syn_thought_r)',partial_match = TRUE,n=bthomax+1)
	forSplinewit=derivatives(witgAge,term='s(cbcl_scr_syn_withdep_r)',partial_match = TRUE,n=bwitmax+1)
	forSplinesoc=derivatives(socgAge,term='s(cbcl_scr_syn_social_r)',partial_match = TRUE,n=bsocmax+1)
	forSplineatt=derivatives(attgAge,term='s(cbcl_scr_syn_attention_r)',partial_match = TRUE,n=battmax+1)
	forSplinerul=derivatives(rulgAge,term='s(cbcl_scr_syn_rulebreak_r)',partial_match = TRUE,n=brulmax+1)
	forSplineagg=derivatives(agggAge,term='s(cbcl_scr_syn_aggressive_r)',partial_match = TRUE,n=baggmax+1)
	# print out fit derivatives
	pDeriv[b,1:(bpmax+1)]=forSplinep$derivative
	intDeriv[b,1:(bimax+1)]=forSplineint$derivative
	extDeriv[b,1:(bemax+1)]=forSplineext$derivative
	somDeriv[b,1:(bsommax+1)]=forSplinesom$derivative
	anxDeriv[b,1:(banxmax+1)]=forSplineanx$derivative
	thoDeriv[b,1:(bthomax+1)]=forSplinetho$derivative
	witDeriv[b,1:(bwitmax+1)]=forSplinewit$derivative
	socDeriv[b,1:(bsocmax+1)]=forSplinesoc$derivative
	attDeriv[b,1:(battmax+1)]=forSplineatt$derivative
	rulDeriv[b,1:(brulmax+1)]=forSplinerul$derivative
	aggDeriv[b,1:(baggmax+1)]=forSplineagg$derivative
	# education version of derivatives
	forSplinep=derivatives(edu_pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE,n=bpmax+1)
	forSplineint=derivatives(edu_intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE,n=bimax+1)
	forSplineext=derivatives(edu_extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE,n=bemax+1)
	forSplinesom=derivatives(edu_somgAge,term='s(cbcl_scr_syn_somatic_r)',partial_match = TRUE,n=bsommax+1)
	forSplineanx=derivatives(edu_anxgAge,term='s(cbcl_scr_syn_anxdep_r)',partial_match = TRUE,n=banxmax+1)
	forSplinetho=derivatives(edu_thogAge,term='s(cbcl_scr_syn_thought_r)',partial_match = TRUE,n=bthomax+1)
	forSplinewit=derivatives(edu_witgAge,term='s(cbcl_scr_syn_withdep_r)',partial_match = TRUE,n=bwitmax+1)
	forSplinesoc=derivatives(edu_socgAge,term='s(cbcl_scr_syn_social_r)',partial_match = TRUE,n=bsocmax+1)
	forSplineatt=derivatives(edu_attgAge,term='s(cbcl_scr_syn_attention_r)',partial_match = TRUE,n=battmax+1)
	forSplinerul=derivatives(edu_rulgAge,term='s(cbcl_scr_syn_rulebreak_r)',partial_match = TRUE,n=brulmax+1)
	forSplineagg=derivatives(edu_agggAge,term='s(cbcl_scr_syn_aggressive_r)',partial_match = TRUE,n=baggmax+1)
	# print out fit derivatives
	edu_pDeriv[b,1:(bpmax+1)]=forSplinep$derivative
	edu_intDeriv[b,1:(bimax+1)]=forSplineint$derivative
	edu_extDeriv[b,1:(bemax+1)]=forSplineext$derivative
	edu_somDeriv[b,1:(bsommax+1)]=forSplinesom$derivative
	edu_anxDeriv[b,1:(banxmax+1)]=forSplineanx$derivative
	edu_thoDeriv[b,1:(bthomax+1)]=forSplinetho$derivative
	edu_witDeriv[b,1:(bwitmax+1)]=forSplinewit$derivative
	edu_socDeriv[b,1:(bsocmax+1)]=forSplinesoc$derivative
	edu_attDeriv[b,1:(battmax+1)]=forSplineatt$derivative
	edu_rulDeriv[b,1:(brulmax+1)]=forSplinerul$derivative
	edu_aggDeriv[b,1:(baggmax+1)]=forSplineagg$derivative
	# both version of derivatives
	forSplinep=derivatives(both_pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE,n=bpmax+1)
	forSplineint=derivatives(both_intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE,n=bimax+1)
	forSplineext=derivatives(both_extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE,n=bemax+1)
	forSplinesom=derivatives(both_somgAge,term='s(cbcl_scr_syn_somatic_r)',partial_match = TRUE,n=bsommax+1)
	forSplineanx=derivatives(both_anxgAge,term='s(cbcl_scr_syn_anxdep_r)',partial_match = TRUE,n=banxmax+1)
	forSplinetho=derivatives(both_thogAge,term='s(cbcl_scr_syn_thought_r)',partial_match = TRUE,n=bthomax+1)
	forSplinewit=derivatives(both_witgAge,term='s(cbcl_scr_syn_withdep_r)',partial_match = TRUE,n=bwitmax+1)
	forSplinesoc=derivatives(both_socgAge,term='s(cbcl_scr_syn_social_r)',partial_match = TRUE,n=bsocmax+1)
	forSplineatt=derivatives(both_attgAge,term='s(cbcl_scr_syn_attention_r)',partial_match = TRUE,n=battmax+1)
	forSplinerul=derivatives(both_rulgAge,term='s(cbcl_scr_syn_rulebreak_r)',partial_match = TRUE,n=brulmax+1)
	forSplineagg=derivatives(both_agggAge,term='s(cbcl_scr_syn_aggressive_r)',partial_match = TRUE,n=baggmax+1)
	# print out fit derivatives
	both_pDeriv[b,1:(bpmax+1)]=forSplinep$derivative
	both_intDeriv[b,1:(bimax+1)]=forSplineint$derivative
	both_extDeriv[b,1:(bemax+1)]=forSplineext$derivative
	both_somDeriv[b,1:(bsommax+1)]=forSplinesom$derivative
	both_anxDeriv[b,1:(banxmax+1)]=forSplineanx$derivative
	both_thoDeriv[b,1:(bthomax+1)]=forSplinetho$derivative
	both_witDeriv[b,1:(bwitmax+1)]=forSplinewit$derivative
	both_socDeriv[b,1:(bsocmax+1)]=forSplinesoc$derivative
	both_attDeriv[b,1:(battmax+1)]=forSplineatt$derivative
	both_rulDeriv[b,1:(brulmax+1)]=forSplinerul$derivative
	both_aggDeriv[b,1:(baggmax+1)]=forSplineagg$derivative
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
	# record adjusted r^2s from models
	pAdjR2[b]=summary(pgAge)$r.sq
	intAdjR2[b]=summary(intgAge)$r.sq
	extAdjR2[b]=summary(extgAge)$r.sq
	somAdjR2[b]=summary(somgAge)$r.sq
	anxAdjR2[b]=summary(anxgAge)$r.sq
	thoAdjR2[b]=summary(thogAge)$r.sq
	witAdjR2[b]=summary(witgAge)$r.sq
	socAdjR2[b]=summary(socgAge)$r.sq
	attAdjR2[b]=summary(attgAge)$r.sq
	rulAdjR2[b]=summary(rulgAge)$r.sq
	aggAdjR2[b]=summary(agggAge)$r.sq
	# record edu adjusted r^2s from models: edu
	edu_pAdjR2[b]=summary(edu_pgAge)$r.sq
	edu_intAdjR2[b]=summary(edu_intgAge)$r.sq
	edu_extAdjR2[b]=summary(edu_extgAge)$r.sq
	edu_somAdjR2[b]=summary(edu_somgAge)$r.sq
	edu_anxAdjR2[b]=summary(edu_anxgAge)$r.sq
	edu_thoAdjR2[b]=summary(edu_thogAge)$r.sq
	edu_witAdjR2[b]=summary(edu_witgAge)$r.sq
	edu_socAdjR2[b]=summary(edu_socgAge)$r.sq
	edu_attAdjR2[b]=summary(edu_attgAge)$r.sq
	edu_rulAdjR2[b]=summary(edu_rulgAge)$r.sq
	edu_aggAdjR2[b]=summary(edu_agggAge)$r.sq
	# record both adjusted r^2s from models: both
	both_pAdjR2[b]=summary(both_pgAge)$r.sq
	both_intAdjR2[b]=summary(both_intgAge)$r.sq
	both_extAdjR2[b]=summary(both_extgAge)$r.sq
	both_somAdjR2[b]=summary(both_somgAge)$r.sq
	both_anxAdjR2[b]=summary(both_anxgAge)$r.sq
	both_thoAdjR2[b]=summary(both_thogAge)$r.sq
	both_witAdjR2[b]=summary(both_witgAge)$r.sq
	both_socAdjR2[b]=summary(both_socgAge)$r.sq
	both_attAdjR2[b]=summary(both_attgAge)$r.sq
	both_rulAdjR2[b]=summary(both_rulgAge)$r.sq
	both_aggAdjR2[b]=summary(both_agggAge)$r.sq
}
# SAVEOUT
# save out INCOME versions
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,somLinBoots,anxLinBoots,thoLinBoots,witLinBoots,socLinBoots,attLinBoots,rulLinBoots,aggLinBoots,pMax,intMax,extMax,somMax,anxMax,thoMax,witMax,socMax,attMax,rulMax,aggMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots_sRE_inc.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,socDeriv,attDeriv,rulDeriv,aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_sRE_inc.rds')
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_sRE_inc.rds')
# save out EDUCATION versions: derivatives
outdf=data.frame(edu_pDeriv,edu_intDeriv,edu_extDeriv,edu_somDeriv,edu_anxDeriv,edu_thoDeriv,edu_witDeriv,edu_socDeriv,edu_attDeriv,edu_rulDeriv,edu_aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_sRE_edu.rds')
# education versions: fit
outdf=data.frame(edu_pFit,edu_intFit,edu_extFit,edu_somFit,edu_anxFit,edu_thoFit,edu_witFit,edu_socFit,edu_attFit,edu_rulFit,edu_aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_sRE_edu.rds')
# save out BOTH versions: derivatives
outdf=data.frame(both_pDeriv,both_intDeriv,both_extDeriv,both_somDeriv,both_anxDeriv,both_thoDeriv,both_witDeriv,both_socDeriv,both_attDeriv,both_rulDeriv,both_aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_sRE_incEdu.rds')
# both versions: fit
outdf=data.frame(both_pFit,both_intFit,both_extFit,both_somFit,both_anxFit,both_thoFit,both_witFit,both_socFit,both_attFit,both_rulFit,both_aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_sRE_incEdu.rds')
# save out adjusted R^2s: income, education, and both
outdf=data.frame(pAdjR2,intAdjR2,extAdjR2,somAdjR2,anxAdjR2,thoAdjR2,witAdjR2,socAdjR2,attAdjR2,rulAdjR2,aggAdjR2,edu_pAdjR2,edu_intAdjR2,edu_extAdjR2,edu_somAdjR2,edu_anxAdjR2,edu_thoAdjR2,edu_witAdjR2,edu_socAdjR2,edu_attAdjR2,edu_rulAdjR2,edu_aggAdjR2,both_pAdjR2,both_intAdjR2,both_extAdjR2,both_somAdjR2,both_anxAdjR2,both_thoAdjR2,both_witAdjR2,both_socAdjR2,both_attAdjR2,both_rulAdjR2,both_aggAdjR2)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpAdjR2Boots_sRE_incEdu.rds')
print('done with g~p fit bootstrapping!')
