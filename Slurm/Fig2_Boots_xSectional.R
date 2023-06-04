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
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','g','subjectkey','interview_age','Grades')]
# get length of df for later
lenDF=dim(masterdf)[1]
### initialize cross-boot vectors
# linear? 1xbootlength
plinBoots=rep(0,10000)
intlinBoots=rep(0,10000)
extlinBoots=rep(0,10000)
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
iMaxVal=max(masterdf$cbcl_scr_syn_internal_r)
emaxVal=max(masterdf$cbcl_scr_syn_external_r)
pDeriv=matrix(0,nrow=10000,ncol=pMaxVal)
intDeriv=matrix(0,nrow=10000,ncol=iMaxVal)
extDeriv=matrix(0,nrow=10000,ncol=emaxVal)
pDerivRaw=matrix(0,nrow=10000,ncol=pMaxVal)
intDerivRaw=matrix(0,nrow=10000,ncol=iMaxVal)
extDerivRaw=matrix(0,nrow=10000,ncol=emaxVal)
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
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
	######## I FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(g~cbcl_scr_syn_totprob_r+s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age),data=bootSamp)
	plinBoots[b]=summary(pgAge)$s.pv[1]
	intgAge<-bam(g~cbcl_scr_syn_internal_r+s(cbcl_scr_syn_internal_r,m=c(2,0))+s(interview_age),data=bootSamp)
	intlinBoots[b]=summary(intgAge)$s.pv[1]
	extgAge<-bam(g~cbcl_scr_syn_external_r+s(cbcl_scr_syn_external_r,m=c(2,0))+s(interview_age),data=bootSamp)
	extlinBoots[b]=summary(extgAge)$s.pv[1]
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	# use PREDICTED VALUES of model fit for saving
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE,n=bpmax)
	forSplineint=derivatives(intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE,n=bimax)
	forSplineext=derivatives(extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE,n=bemax)
	# print out unconverted version
	pDerivRaw[b,1:bpmax]=forSplinep$derivative
	intDerivRaw[b,1:bimax]=forSplineint$derivative
	extDerivRaw[b,1:bemax]=forSplineext$derivative
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
	intMax[b]=bimax
	extMax[b]=bemax
}
# SAVEOUT
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,pMax,intMax,extMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots.rds')
outdf=data.frame(pDerivRaw,intDerivRaw,extDerivRaw)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpRawDerivBoots.rds')
print('done with g~p fit bootstrapping!')
