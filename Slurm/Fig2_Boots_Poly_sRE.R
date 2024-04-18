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

get_derivs <- function(fit) {
  derivative_arr=rep(0,length(fit)-1)
  # Calculate the derivative for each column
  for (i in 1:(length(fit) - 1)) {
    # Calculate the differences in x (assuming a constant difference)
    dx <- 1
    # Calculate the differences in y (predicted values)
    dy <- fit[i + 1] - fit[i]
    # Calculate the derivatives (slopes)
    derivatives <- dy / dx
    # Store the derivatives in the derivative matrix
    derivative_arr[i] <- derivatives
  }
  return(derivative_arr)
}

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)
# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age','site')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert site to factor
masterdf$site<-as.factor(masterdf$site)
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
# note deriv is set to one shorter than gratia method: doesn't resample across caclulated number-to-number difference
pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
somDeriv=matrix(0,nrow=10000,ncol=somMaxVal)
anxDeriv=matrix(0,nrow=10000,ncol=anxMaxVal)
thoDeriv=matrix(0,nrow=10000,ncol=thoMaxVal)
witDeriv=matrix(0,nrow=10000,ncol=witMaxVal)
socDeriv=matrix(0,nrow=10000,ncol=socMaxVal)
attDeriv=matrix(0,nrow=10000,ncol=attMaxVal)
rulDeriv=matrix(0,nrow=10000,ncol=rulMaxVal)
aggDeriv=matrix(0,nrow=10000,ncol=aggMaxVal)
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
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
	pgAge<-bam(g~poly(cbcl_scr_syn_totprob_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	intgAge<-bam(g~poly(cbcl_scr_syn_internal_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	extgAge<-bam(g~poly(cbcl_scr_syn_external_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	somgAge<-bam(g~poly(cbcl_scr_syn_somatic_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	anxgAge<-bam(g~poly(cbcl_scr_syn_anxdep_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	thogAge<-bam(g~poly(cbcl_scr_syn_thought_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	witgAge<-bam(g~poly(cbcl_scr_syn_withdep_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	socgAge<-bam(g~poly(cbcl_scr_syn_social_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	attgAge<-bam(g~poly(cbcl_scr_syn_attention_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	rulgAge<-bam(g~poly(cbcl_scr_syn_rulebreak_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
	agggAge<-bam(g~poly(cbcl_scr_syn_aggressive_r,2)+poly(interview_age,2)+s(site,bs="re"),data=bootSamp)
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
        predictDFp=data.frame(eachPcount,rep(median(masterdf$interview_age),bpmax+1),rep('site16',bpmax+1))
        predictDFint=data.frame(eachIntcount,rep(median(masterdf$interview_age),bimax+1),rep('site16',bimax+1))
        predictDFext=data.frame(eachExtcount,rep(median(masterdf$interview_age),bemax+1),rep('site16',bemax+1))
        predictDFsom=data.frame(eachSomcount,rep(median(masterdf$interview_age),bsommax+1),rep('site16',bsommax+1))
        predictDFanx=data.frame(eachAnxcount,rep(median(masterdf$interview_age),banxmax+1),rep('site16',banxmax+1))
        predictDFtho=data.frame(eachThocount,rep(median(masterdf$interview_age),bthomax+1),rep('site16',bthomax+1))
        predictDFwit=data.frame(eachWitcount,rep(median(masterdf$interview_age),bwitmax+1),rep('site16',bwitmax+1))
        predictDFsoc=data.frame(eachSoccount,rep(median(masterdf$interview_age),bsocmax+1),rep('site16',bsocmax+1))
        predictDFatt=data.frame(eachAttcount,rep(median(masterdf$interview_age),battmax+1),rep('site16',battmax+1))
        predictDFrul=data.frame(eachRulcount,rep(median(masterdf$interview_age),brulmax+1),rep('site16',brulmax+1))
        predictDFagg=data.frame(eachAggcount,rep(median(masterdf$interview_age),baggmax+1),rep('site16',baggmax+1))
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
	forSplinep=get_derivs(forFitP)
	forSplineint=get_derivs(forFitInt)
	forSplineext=get_derivs(forFitExt)
	forSplinesom=get_derivs(forFitSom)
	forSplineanx=get_derivs(forFitAnx)
	forSplinetho=get_derivs(forFitTho)
	forSplinewit=get_derivs(forFitWit)
	forSplinesoc=get_derivs(forFitSoc)
	forSplineatt=get_derivs(forFitAtt)
	forSplinerul=get_derivs(forFitRul)
	forSplineagg=get_derivs(forFitAgg)
	# print out fit derivatives
	pDeriv[b,1:(bpmax)]=forSplinep
	intDeriv[b,1:(bimax)]=forSplineint
	extDeriv[b,1:(bemax)]=forSplineext
	somDeriv[b,1:(bsommax)]=forSplinesom
	anxDeriv[b,1:(banxmax)]=forSplineanx
	thoDeriv[b,1:(bthomax)]=forSplinetho
	witDeriv[b,1:(bwitmax)]=forSplinewit
	socDeriv[b,1:(bsocmax)]=forSplinesoc
	attDeriv[b,1:(battmax)]=forSplineatt
	rulDeriv[b,1:(brulmax)]=forSplinerul
	aggDeriv[b,1:(baggmax)]=forSplineagg
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
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots_poly_sRE.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv,somDeriv,anxDeriv,thoDeriv,witDeriv,socDeriv,attDeriv,rulDeriv,aggDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots_poly_sRE.rds')
outdf=data.frame(pFit,intFit,extFit,somFit,anxFit,thoFit,witFit,socFit,attFit,rulFit,aggFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_poly_sRE.rds')

print('done with g~p fit bootstrapping!')
