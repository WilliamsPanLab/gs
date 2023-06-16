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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','parentPcount','g','subjectkey','interview_age')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r=as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r=as.numeric(masterdf$cbcl_scr_syn_external_r)
masterdf$parentPcount=as.numeric(masterdf$parentPcount)
### initialize cross-boot vectors
# linear? 1xbootlength
plinBoots=rep(0,10000)
intlinBoots=rep(0,10000)
extlinBoots=rep(0,10000)
### NOTE THAT PARENT P IS ACTUALLY G AS OUTCOME VARIABLE
parentPlinBoots=rep(0,10000)
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$parentPcount)
# and initialize derivative vectors
pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
intDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
extDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
# again, parent p nomenclature represents g as outcome variable
parentPderiv=matrix(0,nrow=10000,ncol=pMaxVal)
# predicted values: set to maximum value for ncol
pFit=matrix(0,nrow=10000,ncol=pMaxVal)
intFit=matrix(0,nrow=10000,ncol=pMaxVal)
extFit=matrix(0,nrow=10000,ncol=pMaxVal)
parentPFit=matrix(0,nrow=10000,ncol=pMaxVal)
# deviance explained in g from parent P vs child measures
pDevExpl_g=rep(0,10000)
intDevExpl_g=rep(0,10000)
extDevExpl_g=rep(0,10000)
parentPDevExpl_g=rep(0,10000)
# note only 1 max because it is parent P as predictor
pMax=rep(0,10000)
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
	bpmax=max(bootSamp$parentPcount)
	######## I FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(cbcl_scr_syn_totprob_r~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age),data=bootSamp,family=nb())
	plinBoots[b]=summary(pgAge)$s.pv[1]
	intgAge<-bam(cbcl_scr_syn_internal_r~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age),data=bootSamp,family=nb())
	intlinBoots[b]=summary(intgAge)$s.pv[1]
	extgAge<-bam(cbcl_scr_syn_external_r~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age),data=bootSamp,family=nb())
	extlinBoots[b]=summary(extgAge)$s.pv[1]
	parentPgAge<-bam(g~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age),data=bootSamp)
	parentPlinBoots[b]=summary(parentPgAge)$s.pv[1]
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
	pgAge<-bam(cbcl_scr_syn_totprob_r~s(parentPcount)+s(interview_age),data=bootSamp,family=nb())
	intgAge<-bam(cbcl_scr_syn_internal_r~s(parentPcount)+s(interview_age),data=bootSamp,family=nb())
	extgAge<-bam(cbcl_scr_syn_external_r~s(parentPcount)+s(interview_age),data=bootSamp,family=nb())
	gParentp<-bam(g~s(parentPcount)+s(interview_age),data=bootSamp)
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(1:bpmax)
	# set age to to median for predict df
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),bpmax))
	# set colnames so predict can work
	colnames(predictDFp)=c('parentPcount','interview_age')
	# predict
	forFitP=predict(pgAge,predictDFp)
	forFitint=predict(intgAge,predictDFp)
	forFitext=predict(extgAge,predictDFp)
	forFitParentP=predict(gParentp,predictDFp)
	# print out fit
	pFit[b,1:bpmax]=forFitP
	intFit[b,1:bpmax]=forFitint
	extFit[b,1:bpmax]=forFitext
	parentPFit[b,1:bpmax]=forFitParentP
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(parentPcount)',partial_match = TRUE,n=bpmax)
	forSplineint=derivatives(intgAge,term='s(parentPcount)',partial_match = TRUE,n=bpmax)
	forSplineext=derivatives(extgAge,term='s(parentPcount)',partial_match = TRUE,n=bpmax)
	forSplineParentP=derivatives(gParentp,term='s(parentPcount)',partial_match = TRUE,n=bpmax)
	# print out fit derivatives
	pDeriv[b,1:bpmax]=forSplinep$derivative
	intDeriv[b,1:bpmax]=forSplineint$derivative
	extDeriv[b,1:bpmax]=forSplineext$derivative
	parentPderiv[b,1:bpmax]=forSplineParentP$derivative
	####### III COMPARE DEVIANCE EXPLAINED IN G FROM PARENT P VS CHILD P, INT EXT
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	pDevExpl_g[b]=summary(pgAge)$dev.expl
	intDevExpl_g[b]=summary(intgAge)$dev.expl
	extDevExpl_g[b]=summary(extgAge)$dev.expl
	parentPDevExpl_g[b]=summary(gParentp)$dev.expl
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
}
# SAVEOUT
# save out version with all cbcl factors
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,parentPlinBoots)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gParentpBoots.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv,parentPderiv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gParentpDerivBoots.rds')
outdf=data.frame(pFit,intFit,extFit,parentPFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gParentpFitBoots.rds')
outdf=data.frame(pDevExpl_g,intDevExpl_g,extDevExpl_g,parentPDevExpl_g)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gParentpDevExplBoots.rds')
print('done with g~p fit bootstrapping!')
