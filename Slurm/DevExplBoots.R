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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','parentPcount','g','subjectkey','interview_age','sex','income','Grades')]
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
devExplBoots_Grades=rep(0,10000)
devExplBoots_gparentP=rep(0,10000)
devExplBoots_GradesparentP=rep(0,10000)
# initialize AIC difference vectors for deviance explained full vs. reduced
AICdiff_g_v_grades=rep(0,10000)
AICdiff_g_v_gParentP=rep(0,10000)
AICdiff_g_v_gradesParentP=rep(0,10000)
AICdiff_grades_v_gParentP=rep(0,10000)
AICdiff_grades_v_gradesParentP=rep(0,10000)
AICdiff_gParentP_v_gradesParentP=rep(0,10000)
p_permutedTstats_AvB=rep(0,10000)
p_permutedTstats_AvC=rep(0,10000)
p_permutedTstats_AvD=rep(0,10000)
p_permutedTstats_BvC=rep(0,10000)
p_permutedTstats_BvD=rep(0,10000)
p_permutedTstats_CvD=rep(0,10000)
int_permutedTstats_AvB=rep(0,10000)
int_permutedTstats_AvC=rep(0,10000)
int_permutedTstats_AvD=rep(0,10000)
int_permutedTstats_BvC=rep(0,10000)
int_permutedTstats_BvD=rep(0,10000)
int_permutedTstats_CvD=rep(0,10000)
ext_permutedTstats_AvB=rep(0,10000)
ext_permutedTstats_AvC=rep(0,10000)
ext_permutedTstats_AvD=rep(0,10000)
ext_permutedTstats_BvC=rep(0,10000)
ext_permutedTstats_BvD=rep(0,10000)
ext_permutedTstats_CvD=rep(0,10000)
# collapse all 4s and 5s grades into 4
masterdf$Grades[masterdf$Grades==5]=4

set.seed(1)
for (b in 1:1250){
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
	#### make psuedogrades (subject-wise permutations, not row-wise) to get permuted t-stats between grades on cbcl scores
	# get unique instances of subjects, note duplcated and !duplicated are same because 2x instances of each subj
	unqSubjs=bootSamp[duplicated(bootSamp$subjectkey),]
	# get shuffled version of Grades column as pseudogrades for each unqiue subject
	unqSubjs$Pseudogrades=sample(unqSubjs$Grades)
	# re-initialize pseudogrades in bootSamp
	bootSamp$Pseudogrades = bootSamp$Grades
	# apply each subjects pseudograde back to full samlpe with both observations per subject
	for (i in 1:length(unqSubjs$subjectkey)){
		bootSamp$Pseudogrades[bootSamp$subjectkey==unqSubjs$subjectkey[i]]=unqSubjs$Pseudogrades[i]
	}
	# record permuted-grades t-stats between individual pairs of grades on cbcl scores for null reference distribution
	# get t-stat between A and B
	p_permutedTstats_AvB[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_AvB[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_AvB[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_external_r'])$statistic
	# get t-stat between A and C
	p_permutedTstats_AvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_AvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_AvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_external_r'])$statistic
	# get t-stat between A and D
	p_permutedTstats_AvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_AvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_AvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==1,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_external_r'])$statistic
	# get t-stat between B and C
	p_permutedTstats_BvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_BvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_BvC[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_external_r'])$statistic
	# get t-stat between B and D
	p_permutedTstats_BvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_BvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_BvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==2,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_external_r'])$statistic
	# get t-stat between C and D
	p_permutedTstats_CvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_totprob_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_totprob_r'])$statistic
	int_permutedTstats_CvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_internal_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_internal_r'])$statistic
	ext_permutedTstats_CvD[b]=t.test(bootSamp[bootSamp$Pseudogrades==3,'cbcl_scr_syn_external_r'],bootSamp[bootSamp$Pseudogrades==4,'cbcl_scr_syn_external_r'])$statistic

	# deviance explained in p
	gmod=bam(cbcl_scr_syn_totprob_r~s(g,k=4),data=bootSamp,family=nb())
	devExplBoots_g[b]=summary(gmod)$dev.expl
	# deviance explained by Grades
	GradesMod=bam(cbcl_scr_syn_totprob_r~s(g,k=4)+Grades,data=bootSamp,family=nb())
	devExplBoots_Grades[b]=summary(GradesMod)$dev.expl
	# deviance explained by g and parentPcount
	gparentPMod=bam(cbcl_scr_syn_totprob_r~s(g,k=4)+s(parentPcount,k=4),data=bootSamp,family=nb())
	devExplBoots_gparentP[b]=summary(gparentPMod)$dev.expl
	# deviance explained by parentPcount and Grades
	GradesparentPMod=bam(cbcl_scr_syn_totprob_r~s(parentPcount,k=4)+Grades,data=bootSamp,family=nb())
	devExplBoots_GradesparentP[b]=summary(GradesparentPMod)$dev.expl
	#### get AIC differences between models
	AICdiff_g_v_grades[b]=AIC(gmod)-AIC(GradesMod)
	AICdiff_g_v_gParentP[b]=AIC(gmod)-AIC(gparentPMod)
	AICdiff_g_v_gradesParentP[b]=AIC(gmod)-AIC(GradesparentPMod)
	AICdiff_grades_v_gParentP[b]=AIC(GradesMod)-AIC(gparentPMod)
	AICdiff_grades_v_gradesParentP[b]=AIC(GradesMod)-AIC(GradesparentPMod)
	AICdiff_gParentP_v_gradesParentP[b]=AIC(gparentPMod)-AIC(GradesparentPMod)
}
#### get true difference-of-grades t-stat
p_permutedTstats_AvB[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_totprob_r'],masterdf[masterdf$Grades==2,'cbcl_scr_syn_totprob_r'])$statistic
int_permutedTstats_AvB[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_internal_r'],masterdf[masterdf$Grades==2,'cbcl_scr_syn_internal_r'])$statistic
ext_permutedTstats_AvB[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_external_r'],masterdf[masterdf$Grades==2,'cbcl_scr_syn_external_r'])$statistic
p_permutedTstats_AvC[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_totprob_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_totprob_r'])$statistic
int_permutedTstats_AvC[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_internal_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_internal_r'])$statistic
ext_permutedTstats_AvC[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_external_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_external_r'])$statistic
p_permutedTstats_AvD[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_totprob_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_totprob_r'])$statistic
int_permutedTstats_AvD[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_internal_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_internal_r'])$statistic
ext_permutedTstats_AvD[10001]=t.test(masterdf[masterdf$Grades==1,'cbcl_scr_syn_external_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_external_r'])$statistic
p_permutedTstats_BvC[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_totprob_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_totprob_r'])$statistic
int_permutedTstats_BvC[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_internal_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_internal_r'])$statistic
ext_permutedTstats_BvC[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_external_r'],masterdf[masterdf$Grades==3,'cbcl_scr_syn_external_r'])$statistic
p_permutedTstats_BvD[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_totprob_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_totprob_r'])$statistic
int_permutedTstats_BvD[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_internal_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_internal_r'])$statistic
ext_permutedTstats_BvD[10001]=t.test(masterdf[masterdf$Grades==2,'cbcl_scr_syn_external_r'],masterdf[masterdf$Grades==4,'cbcl_scr_syn_external_r'])$statistic
# SAVEOUT

# saveout all deviance explained vectors in one dataframe
outdf=data.frame(devExplBoots_g,devExplBoots_Grades,devExplBoots_gparentP,devExplBoots_GradesparentP)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3-5DevExpl1.rds')

# saveout all AIC differences in one dataframe
outdf=data.frame(AICdiff_g_v_grades,AICdiff_g_v_gParentP,AICdiff_g_v_gradesParentP,AICdiff_grades_v_gParentP,AICdiff_grades_v_gradesParentP,AICdiff_gParentP_v_gradesParentP)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3-5AICdiff1.rds')

# saveout all t-statistics in one dataframe
outdf=data.frame(p_permutedTstats_AvB,int_permutedTstats_AvB,ext_permutedTstats_AvB,p_permutedTstats_AvC,int_permutedTstats_AvC,ext_permutedTstats_AvC,p_permutedTstats_AvD,int_permutedTstats_AvD,ext_permutedTstats_AvD,p_permutedTstats_BvC,int_permutedTstats_BvC,ext_permutedTstats_BvC,p_permutedTstats_BvD,int_permutedTstats_BvD,ext_permutedTstats_BvD)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3-5Tstats1.rds')
