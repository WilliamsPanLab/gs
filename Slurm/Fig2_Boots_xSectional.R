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
# get 99.9th percentiles for scaling 1:200 matrix
# calculate 99.9th percentile of P int and Ext
p999<-quantile(masterdf$cbcl_scr_syn_totprob_r,.999)
i999<-quantile(masterdf$cbcl_scr_syn_internal_r,.999)
e999<-quantile(masterdf$cbcl_scr_syn_external_r,.999)
g999<-quantile(masterdf$g,.999)
gLow999<-quantile(masterdf$g,.001)
### initialize cross-boot vectors
# linear? 1xbootlength
plinBoots=rep(0,10000)
intlinBoots=rep(0,10000)
extlinBoots=rep(0,10000)
# predicted derivatives? note 200 will scale to maximum of boot samp.
pDeriv=matrix(0,nrow=10000,ncol=200)
intDeriv=matrix(0,nrow=10000,ncol=200)
extDeriv=matrix(0,nrow=10000,ncol=200)
pDerivRaw=matrix(0,nrow=10000,ncol=200)
intDerivRaw=matrix(0,nrow=10000,ncol=200)
extDerivRaw=matrix(0,nrow=10000,ncol=200)
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
# DevExpl in p, int, and ext by grades and g respectively, across log transform and nb
pDevExpl_g=rep(0,10000)
intDevExpl_g=rep(0,10000)
extDevExpl_g=rep(0,10000)
pDevExpl_grades=rep(0,10000)
intDevExpl_grades=rep(0,10000)
extDevExpl_grades=rep(0,10000)
# with negative binomial
nbpDevExpl_g=rep(0,10000)
nbintDevExpl_g=rep(0,10000)
nbextDevExpl_g=rep(0,10000)
nbpDevExpl_grades=rep(0,10000)
nbintDevExpl_grades=rep(0,10000)
nbextDevExpl_grades=rep(0,10000)
# with log transform +1
logpDevExpl_g=rep(0,10000)
logintDevExpl_g=rep(0,10000)
logextDevExpl_g=rep(0,10000)
logpDevExpl_grades=rep(0,10000)
logintDevExpl_grades=rep(0,10000)
logextDevExpl_grades=rep(0,10000)
# conduct log transforms
masterdf$cbcl_scr_syn_totprob_rL=log(masterdf$cbcl_scr_syn_totprob_r+1)
masterdf$cbcl_scr_syn_internal_rL=log(masterdf$cbcl_scr_syn_internal_r+1)
masterdf$cbcl_scr_syn_external_rL=log(masterdf$cbcl_scr_syn_external_r+1)
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
	# use DERIVATIVES of model fit for saving
	forSplinep=derivatives(pgAge,term='s(cbcl_scr_syn_totprob_r)',partial_match = TRUE)
	forSplineint=derivatives(intgAge,term='s(cbcl_scr_syn_internal_r)',partial_match = TRUE)
	forSplineext=derivatives(extgAge,term='s(cbcl_scr_syn_external_r)',partial_match = TRUE)
	# print out unconverted version
	pDerivRaw[b,]=forSplinep$data
	intDerivRaw[b,]=forSplineint$data
	extDerivRaw[b,]=forSplineext$data
	# print out max of unconverted versions to anchor em later
	pMax[b]=max(bootSamp$cbcl_scr_syn_totprob_r)
	intMax[b]=max(bootSamp$cbcl_scr_syn_internal_r)
	extMax[b]=max(bootSamp$	cbcl_scr_syn_external_r)
	# convert to reflect 0-99.9th percentile of true sample scale
	convSplinep<-approx(x=forSplinep$data,y=forSplinep$derivative,xout=seq(0,p999,length.out=200))
	convSplineint<-approx(x=forSplineint$data,y=forSplineint$derivative,xout=seq(0,i999,length.out=200))
	convSplineext<-approx(x=forSplineext$data,y=forSplineext$derivative,xout=seq(0,e999,length.out=200))
	# save out forspline to a matrix
	pDeriv[b,]=convSplinep$y
	intDeriv[b,]=convSplineint$y
	extDeriv[b,]=convSplineext$y
	######## III GRADES DEV EXPLAINED VS G DEV EXPLAINED : symptoms as outcome
	# ACROSS P INT EXT
	pgAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age),data=bootSamp)
	intgAge<-bam(cbcl_scr_syn_internal_r~s(g)+s(interview_age),data=bootSamp)
	extgAge<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age),data=bootSamp)
	pDevExpl_g[b]=summary(pgAge)$dev.expl
	intDevExpl_g[b]=summary(intgAge)$dev.expl
	extDevExpl_g[b]=summary(extgAge)$dev.expl
	# negative binomial link function version
	pgAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age),data=bootSamp,family=nb())
	intgAge<-bam(cbcl_scr_syn_internal_r~s(g)+s(interview_age),data=bootSamp,family=nb())
	extgAge<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age),data=bootSamp,family=nb())
	nbpDevExpl_g[b]=summary(pgAge)$dev.expl
	nbintDevExpl_g[b]=summary(intgAge)$dev.expl
	nbextDevExpl_g[b]=summary(extgAge)$dev.expl
	# log transformed version
	pgAge<-bam(cbcl_scr_syn_totprob_rL~s(g)+s(interview_age),data=bootSamp)
	intgAge<-bam(cbcl_scr_syn_internal_rL~s(g)+s(interview_age),data=bootSamp)
	extgAge<-bam(cbcl_scr_syn_external_rL~s(g)+s(interview_age),data=bootSamp)
	logpDevExpl_g[b]=summary(pgAge)$dev.expl
	logintDevExpl_g[b]=summary(intgAge)$dev.expl
	logextDevExpl_g[b]=summary(extgAge)$dev.expl
	### now for grades ###
	pgAge<-bam(cbcl_scr_syn_totprob_r~Grades+s(interview_age),data=bootSamp)
        intgAge<-bam(cbcl_scr_syn_internal_r~Grades+s(interview_age),data=bootSamp)
        extgAge<-bam(cbcl_scr_syn_external_r~Grades+s(interview_age),data=bootSamp)
        pDevExpl_grades[b]=summary(pgAge)$dev.expl
        intDevExpl_grades[b]=summary(intgAge)$dev.expl
        extDevExpl_grades[b]=summary(extgAge)$dev.expl
        # negative binomial link function version
        pgAge<-bam(cbcl_scr_syn_totprob_r~Grades+s(interview_age),data=bootSamp,family=nb())
        intgAge<-bam(cbcl_scr_syn_internal_r~Grades+s(interview_age),data=bootSamp,family=nb())
        extgAge<-bam(cbcl_scr_syn_external_r~Grades+s(interview_age),data=bootSamp,family=nb())
        nbpDevExpl_grades[b]=summary(pgAge)$dev.expl
        nbintDevExpl_grades[b]=summary(intgAge)$dev.expl
        nbextDevExpl_grades[b]=summary(extgAge)$dev.expl
        # log transformed version
        pgAge<-bam(cbcl_scr_syn_totprob_rL~Grades+s(interview_age),data=bootSamp)
        intgAge<-bam(cbcl_scr_syn_internal_rL~Grades+s(interview_age),data=bootSamp)
        extgAge<-bam(cbcl_scr_syn_external_rL~Grades+s(interview_age),data=bootSamp)
        logpDevExpl_grades[b]=summary(pgAge)$dev.expl
        logintDevExpl_grades[b]=summary(intgAge)$dev.expl
        logextDevExpl_grades[b]=summary(extgAge)$dev.expl
}
# SAVEOUT
outdf=data.frame(plinBoots,intlinBoots,extlinBoots,pDevExpl_g,intDevExpl_g,extDevExpl_g,nbpDevExpl_g,nbintDevExpl_g,nbextDevExpl_g,logpDevExpl_g,logintDevExpl_g,logextDevExpl_g,pDevExpl_grades,intDevExpl_grades,extDevExpl_grades,nbpDevExpl_grades,nbintDevExpl_grades,nbextDevExpl_grades,logpDevExpl_grades,logintDevExpl_grades,logextDevExpl_grades,pMax,intMax,extMax)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots.rds')
outdf=data.frame(pDeriv,intDeriv,extDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpDerivBoots.rds')
outdf=data.frame(pDerivRaw,intDerivRaw,extDerivRaw)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpRawDerivBoots.rds')
print('done with g~p fit bootstrapping!')
