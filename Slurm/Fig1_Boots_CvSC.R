library(mgcv)
library(pammtools)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)
# load in clinical vs subclinical cutoffs
CvSC=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_ClinCutoffs.rds')
# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert all cbcl scores to numeric
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
# initialize output betas
pBetaDiff=rep(0,10000)
intBetaDiff=rep(0,10000)
extBetaDiff=rep(0,10000)
somBetaDiff=rep(0,10000)
anxBetaDiff=rep(0,10000)
thoBetaDiff=rep(0,10000)
witBetaDiff=rep(0,10000)
socBetaDiff=rep(0,10000)
attBetaDiff=rep(0,10000)
rulBetaDiff=rep(0,10000)
aggBetaDiff=rep(0,10000)
# get number of subjects below and above cutoff for each cbcl factor
pAbove=sum(masterdf$cbcl_scr_syn_totprob_r>CvSC$Pc)
pBelow=sum(masterdf$cbcl_scr_syn_totprob_r<CvSC$Pbc)
intAbove=sum(masterdf$cbcl_scr_syn_internal_r>CvSC$Ic)
intBelow=sum(masterdf$cbcl_scr_syn_internal_r<CvSC$Ibc)
extAbove=sum(masterdf$cbcl_scr_syn_external_r>CvSC$Ec)
extBelow=sum(masterdf$cbcl_scr_syn_external_r<CvSC$Ebc)
somAbove=sum(masterdf$cbcl_scr_syn_somatic_r>CvSC$SomC)
somBelow=sum(masterdf$cbcl_scr_syn_somatic_r<CvSC$SomBc)
anxAbove=sum(masterdf$cbcl_scr_syn_anxdep_r>CvSC$AnxC)
anxBelow=sum(masterdf$cbcl_scr_syn_anxdep_r<CvSC$AnxBc)
thoAbove=sum(masterdf$cbcl_scr_syn_thought_r>CvSC$ThoC)
thoBelow=sum(masterdf$cbcl_scr_syn_thought_r<CvSC$ThoBc)
witAbove=sum(masterdf$cbcl_scr_syn_withdep_r>CvSC$WitC)
witBelow=sum(masterdf$cbcl_scr_syn_withdep_r<CvSC$WitBc)
socAbove=sum(masterdf$cbcl_scr_syn_social_r>CvSC$SocC)
socBelow=sum(masterdf$cbcl_scr_syn_social_r<CvSC$SocBc)
attAbove=sum(masterdf$cbcl_scr_syn_attention_r>CvSC$AttC)
attBelow=sum(masterdf$cbcl_scr_syn_attention_r<CvSC$AttBc)
rulAbove=sum(masterdf$cbcl_scr_syn_rulebreak_r>CvSC$RulC)
rulBelow=sum(masterdf$cbcl_scr_syn_rulebreak_r<CvSC$RulBc)
aggAbove=sum(masterdf$cbcl_scr_syn_aggressive_r>CvSC$AggC)
aggBelow=sum(masterdf$cbcl_scr_syn_aggressive_r<CvSC$AggBc)
# separate DF by subclinical or clinical
pAboveDF=masterdf[masterdf$cbcl_scr_syn_totprob_r>CvSC$Pc,]
pBelowDF=masterdf[masterdf$cbcl_scr_syn_totprob_r<CvSC$Pbc,]
intAboveDF=masterdf[masterdf$cbcl_scr_syn_internal_r>CvSC$Ic,]
intBelowDF=masterdf[masterdf$cbcl_scr_syn_internal_r<CvSC$Ibc,]
extAboveDF=masterdf[masterdf$cbcl_scr_syn_external_r>CvSC$Ec,]
extBelowDF=masterdf[masterdf$cbcl_scr_syn_external_r<CvSC$Ebc,]
somAboveDF=masterdf[masterdf$cbcl_scr_syn_somatic_r>CvSC$SomC,]
somBelowDF=masterdf[masterdf$cbcl_scr_syn_somatic_r<CvSC$SomBc,]
anxAboveDF=masterdf[masterdf$cbcl_scr_syn_anxdep_r>CvSC$AnxC,]
anxBelowDF=masterdf[masterdf$cbcl_scr_syn_anxdep_r<CvSC$AnxBc,]
thoAboveDF=masterdf[masterdf$cbcl_scr_syn_thought_r>CvSC$ThoC,]
thoBelowDF=masterdf[masterdf$cbcl_scr_syn_thought_r<CvSC$ThoBc,]
witAboveDF=masterdf[masterdf$cbcl_scr_syn_withdep_r>CvSC$WitC,]
witBelowDF=masterdf[masterdf$cbcl_scr_syn_withdep_r<CvSC$WitBc,]
socAboveDF=masterdf[masterdf$cbcl_scr_syn_social_r>CvSC$SocC,]
socBelowDF=masterdf[masterdf$cbcl_scr_syn_social_r<CvSC$SocBc,]
attAboveDF=masterdf[masterdf$cbcl_scr_syn_attention_r>CvSC$AttC,]
attBelowDF=masterdf[masterdf$cbcl_scr_syn_attention_r<CvSC$AttBc,]
rulAboveDF=masterdf[masterdf$cbcl_scr_syn_rulebreak_r>CvSC$RulC,]
rulBelowDF=masterdf[masterdf$cbcl_scr_syn_rulebreak_r<CvSC$RulBc,]
aggAboveDF=masterdf[masterdf$cbcl_scr_syn_aggressive_r>CvSC$AggC,]
aggBelowDF=masterdf[masterdf$cbcl_scr_syn_aggressive_r<CvSC$AggBc,]
# get real betas for linear fit within each cbcl subscale
PClinBeta=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pAboveDF)$coefficients[2]
PsubBeta=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pBelowDF)$coefficients[2]
IntClinBeta=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intAboveDF)$coefficients[2]
IntsubBeta=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intBelowDF)$coefficients[2]
ExtClinBeta=lm(g~cbcl_scr_syn_external_r+interview_age,data=extAboveDF)$coefficients[2]
ExtsubBeta=lm(g~cbcl_scr_syn_external_r+interview_age,data=extBelowDF)$coefficients[2]
SomClinBeta=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somAboveDF)$coefficients[2]
SomsubBeta=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somBelowDF)$coefficients[2]
AnxClinBeta=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxAboveDF)$coefficients[2]
AnxsubBeta=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxBelowDF)$coefficients[2]
ThoClinBeta=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoAboveDF)$coefficients[2]
ThosubBeta=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoBelowDF)$coefficients[2]
WitClinBeta=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witAboveDF)$coefficients[2]
WitsubBeta=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witBelowDF)$coefficients[2]
SocClinBeta=lm(g~cbcl_scr_syn_social_r+interview_age,data=socAboveDF)$coefficients[2]
SocsubBeta=lm(g~cbcl_scr_syn_social_r+interview_age,data=socBelowDF)$coefficients[2]
AttClinBeta=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attAboveDF)$coefficients[2]
AttsubBeta=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attBelowDF)$coefficients[2]
RulClinBeta=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulAboveDF)$coefficients[2]
RulsubBeta=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulBelowDF)$coefficients[2]
AggClinBeta=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggAboveDF)$coefficients[2]
AggsubBeta=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggBelowDF)$coefficients[2]
# get observed difference of betas within each subscale
Pdiff=PsubBeta-PClinBeta
Intdiff=IntsubBeta-IntClinBeta
Extdiff=ExtsubBeta-ExtClinBeta
Somdiff=SomsubBeta-SomClinBeta
Anxdiff=AnxsubBeta-AnxClinBeta
Thodiff=ThosubBeta-ThoClinBeta
Witdiff=WitsubBeta-WitClinBeta
Socdiff=SocsubBeta-SocClinBeta
Attdiff=AttsubBeta-AttClinBeta
Ruldiff=RulsubBeta-RulClinBeta
Aggdiff=AggsubBeta-AggClinBeta
# for each permutation
for (b in 1:10000){
	print(b)
	# make a random clin vs. subclin split based on real number of rows beneath and above thresholds
	pAboveDF=masterdf[sample(1:nrow(masterdf),pAbove),]
	pBelowDF=masterdf[sample(1:nrow(masterdf),pBelow),]
	intAboveDF=masterdf[sample(1:nrow(masterdf),intAbove),]
	intBelowDF=masterdf[sample(1:nrow(masterdf),intBelow),]
	extAboveDF=masterdf[sample(1:nrow(masterdf),extAbove),]
	extBelowDF=masterdf[sample(1:nrow(masterdf),extBelow),]
	somAboveDF=masterdf[sample(1:nrow(masterdf),somAbove),]
	somBelowDF=masterdf[sample(1:nrow(masterdf),somBelow),]
	anxAboveDF=masterdf[sample(1:nrow(masterdf),anxAbove),]
	anxBelowDF=masterdf[sample(1:nrow(masterdf),anxBelow),]
	thoAboveDF=masterdf[sample(1:nrow(masterdf),thoAbove),]
	thoBelowDF=masterdf[sample(1:nrow(masterdf),thoBelow),]
	witAboveDF=masterdf[sample(1:nrow(masterdf),witAbove),]
	witBelowDF=masterdf[sample(1:nrow(masterdf),witBelow),]
	socAboveDF=masterdf[sample(1:nrow(masterdf),socAbove),]
	socBelowDF=masterdf[sample(1:nrow(masterdf),socBelow),]
	attAboveDF=masterdf[sample(1:nrow(masterdf),attAbove),]
	attBelowDF=masterdf[sample(1:nrow(masterdf),attBelow),]
	rulAboveDF=masterdf[sample(1:nrow(masterdf),rulAbove),]
	rulBelowDF=masterdf[sample(1:nrow(masterdf),rulBelow),]
	aggAboveDF=masterdf[sample(1:nrow(masterdf),aggAbove),]
	aggBelowDF=masterdf[sample(1:nrow(masterdf),aggBelow),]
	# get betas for each permuted split
	pAboveBeta=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pAboveDF)$coefficients[2]
	pBelowBeta=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pBelowDF)$coefficients[2]
	intAboveBeta=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intAboveDF)$coefficients[2]
	intBelowBeta=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intBelowDF)$coefficients[2]
	extAboveBeta=lm(g~cbcl_scr_syn_external_r+interview_age,data=extAboveDF)$coefficients[2]
	extBelowBeta=lm(g~cbcl_scr_syn_external_r+interview_age,data=extBelowDF)$coefficients[2]
	somAboveBeta=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somAboveDF)$coefficients[2]
	somBelowBeta=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somBelowDF)$coefficients[2]
	anxAboveBeta=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxAboveDF)$coefficients[2]
	anxBelowBeta=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxBelowDF)$coefficients[2]
	thoAboveBeta=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoAboveDF)$coefficients[2]
	thoBelowBeta=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoBelowDF)$coefficients[2]
	witAboveBeta=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witAboveDF)$coefficients[2]
	witBelowBeta=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witBelowDF)$coefficients[2]
	socAboveBeta=lm(g~cbcl_scr_syn_social_r+interview_age,data=socAboveDF)$coefficients[2]
	socBelowBeta=lm(g~cbcl_scr_syn_social_r+interview_age,data=socBelowDF)$coefficients[2]
	attAboveBeta=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attAboveDF)$coefficients[2]
	attBelowBeta=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attBelowDF)$coefficients[2]
	rulAboveBeta=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulAboveDF)$coefficients[2]
	rulBelowBeta=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulBelowDF)$coefficients[2]
	aggAboveBeta=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggAboveDF)$coefficients[2]
	aggBelowBeta=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggBelowDF)$coefficients[2]
	# get permuted differences
	pBetaDiff[b]=pBelowBeta-pAboveBeta
	intBetaDiff[b]=intBelowBeta-intAboveBeta
	extBetaDiff[b]=extBelowBeta-extAboveBeta
	somBetaDiff[b]=somBelowBeta-somAboveBeta
	anxBetaDiff[b]=anxBelowBeta-anxAboveBeta
	thoBetaDiff[b]=thoBelowBeta-thoAboveBeta
	witBetaDiff[b]=witBelowBeta-witAboveBeta
	socBetaDiff[b]=socBelowBeta-socAboveBeta
	attBetaDiff[b]=attBelowBeta-attAboveBeta
	rulBetaDiff[b]=rulBelowBeta-rulAboveBeta
	aggBetaDiff[b]=aggBelowBeta-aggAboveBeta
}
# add real differences as 10,001th value
pBetaDiff[10001]=Pdiff
intBetaDiff[10001]=Intdiff
extBetaDiff[10001]=Extdiff
somBetaDiff[10001]=Somdiff
anxBetaDiff[10001]=Anxdiff
thoBetaDiff[10001]=Thodiff
witBetaDiff[10001]=Witdiff
socBetaDiff[10001]=Socdiff
attBetaDiff[10001]=Attdiff
rulBetaDiff[10001]=Ruldiff
aggBetaDiff[10001]=Aggdiff
# save out all difference vectors in one dataframe
outdf=data.frame(pBetaDiff,intBetaDiff,extBetaDiff,somBetaDiff,anxBetaDiff,thoBetaDiff,witBetaDiff,socBetaDiff,attBetaDiff,rulBetaDiff,aggBetaDiff)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gp_CvSC_diffs.rds')

print('done with g~p fit bootstrapping!')

