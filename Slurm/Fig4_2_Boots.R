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

# load in data: additional exclusions with 4_2
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf2.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('totcount_y','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','parentPcount','g','subjectkey','interview_age')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r=as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r=as.numeric(masterdf$cbcl_scr_syn_external_r)
masterdf$parentPcount=as.numeric(masterdf$parentPcount)
masterdf$totcount_y=as.numeric(masterdf$totcount_y)
### initialize cross-boot vectors
# deviance explained in g from parent P vs child measures
pDevExpl_g=rep(0,10000)
intDevExpl_g=rep(0,10000)
extDevExpl_g=rep(0,10000)
parentPDevExpl_g=rep(0,10000)
p_sr_DevExpl_g=rep(0,10000)
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
	####### III COMPARE DEVIANCE EXPLAINED IN G FROM PARENT P VS CHILD P, INT EXT
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	intgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age),data=bootSamp)
	extgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	gParentp=bam(g~s(parentPcount)+s(interview_age),data=bootSamp)
	SRpgAge<-bam(g~s(totcount_y)+s(interview_age),data=bootSamp)
	pDevExpl_g[b]=summary(pgAge)$dev.expl
	intDevExpl_g[b]=summary(intgAge)$dev.expl
	extDevExpl_g[b]=summary(extgAge)$dev.expl
	parentPDevExpl_g[b]=summary(gParentp)$dev.expl
	p_sr_DevExpl_g[b]=summary(SRpgAge)$dev.expl
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
}
# SAVEOUT
# save out version with all cbcl factors
outdf=data.frame(pDevExpl_g,intDevExpl_g,extDevExpl_g,parentPDevExpl_g,p_sr_DevExpl_g)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gParentpDevExplBoots_df2.rds'
print('done with g~p fit bootstrapping!')
