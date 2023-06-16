# read in masterdf and calculate deviance explained by subsequent additions to model
library(mgcv)
library(patchwork)
library(gratia)
library(scales)
library(dplyr)
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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','parentPcount','g','subjectkey','interview_age','sex','income')]
# get length of df for later
lenDF=dim(masterdf)[1]
# will need to get full and reduced models for each boot, as well as a null distribution
# in addition to derivatives and fits, save F values of interaction + null distribution F values
# interactions to be tested: p*sex, p*poverty, p*sex*poverty
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

# convert parentPcount to numeric
masterdf$parentPcount=as.numeric(masterdf$parentPcount)
# sex (g in for x to get co-pilot around its NSFW filter) and poverty to factors
masterdf$seg<-as.ordered(masterdf$sex)
masterdf$income<-as.numeric(masterdf$income)
masterdf$poverty=0
masterdf$poverty[masterdf$income<5]=1
masterdf$poverty=as.ordered(masterdf$poverty)
# initialize a vector for each of 10k bootstraps that will hold deviance unexplained by omission of terms from full model
devExplBoots_p=rep(0,10000)
devExplBoots_intext=rep(0,10000)
devExplBoots_seg=rep(0,10000)
devExplBoots_pov=rep(0,10000)
devExplBoots_segpov=rep(0,10000)
devExplBoots_parentP=rep(0,10000)

set.seed(1)
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
	# deviance xplained by p only
	pMod<-bam(g~s(cbcl_scr_syn_totprob_r),data=bootSamp)
	devExplBoots_p[b]<-summary(pMod)$dev.expl
	# deviance explained by p, internalizing and externalizing
	devExplBoots_intext[b]<-summary(bam(g~s(cbcl_scr_syn_totprob_r)+s(cbcl_scr_syn_internal_r)+s(cbcl_scr_syn_external_r),data=bootSamp))$dev.expl
	# deviance explained by p, internalizing, externalizing, and seg/seg interaction
	devExplBoots_seg[b]<-summary(bam(g~s(cbcl_scr_syn_totprob_r)+s(cbcl_scr_syn_internal_r)+s(cbcl_scr_syn_external_r)+seg+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_external_r,by=seg),data=bootSamp))$dev.expl
	# deviance explained by all of that and poverty status/poverty status interaction
	devExplBoots_pov[b]<-summary(bam(g~s(cbcl_scr_syn_totprob_r)+s(cbcl_scr_syn_internal_r)+s(cbcl_scr_syn_external_r)+seg+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_external_r,by=seg)+poverty+s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_external_r,by=poverty),data=bootSamp))$dev.expl
	# deviance explained by all of that and triple interactions with seg and poverty for p int ext
	devExplBoots_segpov[b]<-summary(bam(g~s(cbcl_scr_syn_totprob_r)+s(cbcl_scr_syn_internal_r)+s(cbcl_scr_syn_external_r)+seg+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_external_r,by=seg)+poverty+s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=interaction(seg, poverty))+s(cbcl_scr_syn_internal_r,by=interaction(seg, poverty))+s(cbcl_scr_syn_external_r,by=interaction(seg, poverty)),data=bootSamp))$dev.expl
	# deviance explained by all of that and parentPcount
	devExplBoots_parentP[b]<-summary(bam(g~s(cbcl_scr_syn_totprob_r)+s(cbcl_scr_syn_internal_r)+s(cbcl_scr_syn_external_r)+seg+s(cbcl_scr_syn_totprob_r,by=seg)+s(cbcl_scr_syn_internal_r,by=seg)+s(cbcl_scr_syn_external_r,by=seg)+poverty+s(cbcl_scr_syn_totprob_r,by=poverty)+s(cbcl_scr_syn_internal_r,by=poverty)+s(cbcl_scr_syn_external_r,by=poverty)+s(cbcl_scr_syn_totprob_r,by=interaction(seg, poverty))+s(cbcl_scr_syn_internal_r,by=interaction(seg, poverty))+s(cbcl_scr_syn_external_r,by=interaction(seg, poverty))+s(parentPcount),data=bootSamp))$dev.expl
}
# SAVEOUT
# saveout all deviance explained vectors in one dataframe
outdf=data.frame(devExplBoots_p,devExplBoots_intext,devExplBoots_seg,devExplBoots_pov,devExplBoots_segpov,devExplBoots_parentP)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3-5DevExpl.rds')
