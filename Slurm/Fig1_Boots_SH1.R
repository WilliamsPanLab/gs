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
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_SH1_masterdf.rds')
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
# linear? 1xbootlength
plinBoots=rep(0,10000)
intlinBoots=rep(0,10000)
extlinBoots=rep(0,10000)
somLinBoots=rep(0,10000)
anxLinBoots=rep(0,10000)
thoLinBoots=rep(0,10000)
witLinBoots=rep(0,10000)
socLinBoots=rep(0,10000)
attLinBoots=rep(0,10000)
rulLinBoots=rep(0,10000)
aggLinBoots=rep(0,10000)
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
	bsommax=max(bootSamp$cbcl_scr_syn_somatic_r)
	banxmax=max(bootSamp$cbcl_scr_syn_anxdep_r)
	bthomax=max(bootSamp$cbcl_scr_syn_thought_r)
	bwitmax=max(bootSamp$cbcl_scr_syn_withdep_r)
	bsocmax=max(bootSamp$cbcl_scr_syn_social_r)
	battmax=max(bootSamp$cbcl_scr_syn_attention_r)
	brulmax=max(bootSamp$cbcl_scr_syn_rulebreak_r)
	baggmax=max(bootSamp$cbcl_scr_syn_aggressive_r)
	######## I FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(g~cbcl_scr_syn_totprob_r+s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	plinBoots[b]=summary(pgAge)$s.pv[1]
	intgAge<-bam(g~cbcl_scr_syn_internal_r+s(cbcl_scr_syn_internal_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	intlinBoots[b]=summary(intgAge)$s.pv[1]
	extgAge<-bam(g~cbcl_scr_syn_external_r+s(cbcl_scr_syn_external_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	extlinBoots[b]=summary(extgAge)$s.pv[1]
	somgAge<-bam(g~cbcl_scr_syn_somatic_r+s(cbcl_scr_syn_somatic_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	somLinBoots[b]=summary(somgAge)$s.pv[1]
	anxgAge<-bam(g~cbcl_scr_syn_anxdep_r+s(cbcl_scr_syn_anxdep_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	anxLinBoots[b]=summary(anxgAge)$s.pv[1]
	thogAge<-bam(g~cbcl_scr_syn_thought_r+s(cbcl_scr_syn_thought_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	thoLinBoots[b]=summary(thogAge)$s.pv[1]
	witgAge<-bam(g~cbcl_scr_syn_withdep_r+s(cbcl_scr_syn_withdep_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	witLinBoots[b]=summary(witgAge)$s.pv[1]
	socgAge<-bam(g~cbcl_scr_syn_social_r+s(cbcl_scr_syn_social_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	socLinBoots[b]=summary(socgAge)$s.pv[1]
	attgAge<-bam(g~cbcl_scr_syn_attention_r+s(cbcl_scr_syn_attention_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	attLinBoots[b]=summary(attgAge)$s.pv[1]
	rulgAge<-bam(g~cbcl_scr_syn_rulebreak_r+s(cbcl_scr_syn_rulebreak_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	rulLinBoots[b]=summary(rulgAge)$s.pv[1]
	agggAge<-bam(g~cbcl_scr_syn_aggressive_r+s(cbcl_scr_syn_aggressive_r,m=c(2,0))+s(interview_age,k=4),data=bootSamp)
	aggLinBoots[b]=summary(agggAge)$s.pv[1]
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
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
	# set age to to median for predict df, +1 to max values because 0-max sequence is 1 > max in length.
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax+1))
	predictDFint=data.frame(eachIntcount,rep(median(bootSamp$interview_age),bimax+1))
	predictDFext=data.frame(eachExtcount,rep(median(bootSamp$interview_age),bemax+1))
	predictDFsom=data.frame(eachSomcount,rep(median(bootSamp$interview_age),bsommax+1))
	predictDFanx=data.frame(eachAnxcount,rep(median(bootSamp$interview_age),banxmax+1))
	predictDFtho=data.frame(eachThocount,rep(median(bootSamp$interview_age),bthomax+1))
	predictDFwit=data.frame(eachWitcount,rep(median(bootSamp$interview_age),bwitmax+1))
	predictDFsoc=data.frame(eachSoccount,rep(median(bootSamp$interview_age),bsocmax+1))
	predictDFatt=data.frame(eachAttcount,rep(median(bootSamp$interview_age),battmax+1))
	predictDFrul=data.frame(eachRulcount,rep(median(bootSamp$interview_age),brulmax+1))
	predictDFagg=data.frame(eachAggcount,rep(median(bootSamp$interview_age),baggmax+1))
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
}
# SAVEOUT
# save out version with all cbcl factors
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,somLinBoots,anxLinBoots,thoLinBoots,witLinBoots,socLinBoots,attLinBoots,rulLinBoots,aggLinBoots,pMax,intMax,extMax,somMax,anxMax,thoMax,witMax,socMax,attMax,rulMax,aggMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gp_SH1_Boots.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,socDeriv,attDeriv,rulDeriv,aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gp_SH1_DerivBoots.rds')
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gp_SH1_FitBoots.rds')

print('done with g~p fit bootstrapping!')
