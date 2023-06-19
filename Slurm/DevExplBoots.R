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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','parentPcount','g','subjectkey','interview_age','sex','income','Grades','pc2','pc3')]
# get length of df for later
lenDF=dim(masterdf)[1]
# will need to get full and reduced models for each boot, as well as a null distribution
# in addition to derivatives and fits, save F values of interaction + null distribution F values
# interactions to be tested: p*sex, p*poverty, p*sex*poverty
# convert cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
# convert parentPcount to numeric
masterdf$parentPcount=as.numeric(masterdf$parentPcount)
# initialize a vector for each of 10k bootstraps that will hold deviance unexplained by omission of terms from full model
devExplBoots_g=rep(0,10000)
devExplBoots_pcs=rep(0,10000)
devExplBoots_Grades=rep(0,10000)
devExplBoots_gparentP=rep(0,10000)
devExplBoots_pcsparentP=rep(0,10000)
devExplBoots_GradesparentP=rep(0,10000)
set.seed(1)
for (b in 1:1000){
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
	# deviance explained in p
	gmod=bam(cbcl_scr_syn_totprob_r~s(g),data=bootSamp,family=nb())
	devExplBoots_g[b]=summary(gmod)$dev.expl
	# deviance explained by all 3 pcs
	pcsMod=bam(cbcl_scr_syn_totprob_r~s(g)+s(pc2)+s(pc3),data=bootSamp,family=nb())
	devExplBoots_pcs[b]=summary(pcsMod)$dev.expl
	# deviance explained by Grades
	GradesMod=bam(cbcl_scr_syn_totprob_r~s(g)+Grades,data=bootSamp,family=nb())
	devExplBoots_Grades[b]=summary(GradesMod)$dev.expl
	# deviance explained by g and parentPcount
	gparentPMod=bam(cbcl_scr_syn_totprob_r~s(g)+s(parentPcount),data=bootSamp,family=nb())
	devExplBoots_gparentP[b]=summary(gparentPMod)$dev.expl
	# deviance explained by g, parentPcount, and pcs
	pcsparentPMod=bam(cbcl_scr_syn_totprob_r~s(g)+s(parentPcount)+s(pc2)+s(pc3),data=bootSamp,family=nb())
	devExplBoots_pcsparentP[b]=summary(pcsparentPMod)$dev.expl
	# deviance explained by parentPcount and Grades
	GradesparentPMod=bam(cbcl_scr_syn_totprob_r~s(parentPcount)+Grades,data=bootSamp,family=nb())
	devExplBoots_GradesparentP[b]=summary(GradesparentPMod)$dev.expl
}
# SAVEOUT
# saveout all deviance explained vectors in one dataframe
outdf=data.frame(devExplBoots_g,devExplBoots_pcs,devExplBoots_Grades,devExplBoots_gparentP,devExplBoots_pcsparentP,devExplBoots_GradesparentP)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3-5DevExpl.rds')
