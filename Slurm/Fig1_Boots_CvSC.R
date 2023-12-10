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
# and output distributions of clin and subclin betas
pClinBeta=rep(0,10000)
intClinBeta=rep(0,10000)
extClinBeta=rep(0,10000)
somClinBeta=rep(0,10000)
anxClinBeta=rep(0,10000)
thoClinBeta=rep(0,10000)
witClinBeta=rep(0,10000)
socClinBeta=rep(0,10000)
attClinBeta=rep(0,10000)
rulClinBeta=rep(0,10000)
aggClinBeta=rep(0,10000)
pSubclinBeta=rep(0,10000)
intSubclinBeta=rep(0,10000)
extSubclinBeta=rep(0,10000)
somSubclinBeta=rep(0,10000)
anxSubclinBeta=rep(0,10000)
thoSubclinBeta=rep(0,10000)
witSubclinBeta=rep(0,10000)
socSubclinBeta=rep(0,10000)
attSubclinBeta=rep(0,10000)
rulSubclinBeta=rep(0,10000)
aggSubclinBeta=rep(0,10000)
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
	# get subjects to include in this bootstrap
        BootSubjs=sample(subjs,numSubjs,replace=T)
        ### inefficient but interpretable loop
        # Create an empty dataframe to store the resampled observations
        bootSamp <- data.frame()
        for (j in 1:length(BootSubjs)){
                subject_obs <- masterdf[masterdf$subjectkey == BootSubjs[j], ]
                bootSamp <- rbind(bootSamp, subject_obs)
        }
	# ############ ######## sample split based on randomly sampling subjects, not rows
	unqSubjs=bootSamp[duplicated(bootSamp$subjectkey),]
        # get counts of # in clinical threshold for this boot
	subjects_pAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_totprob_r > CvSC$Pc]
	subjects_intAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_internal_r > CvSC$Ic]
	subjects_extAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_external_r > CvSC$Ec]
	subjects_somAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_somatic_r > CvSC$SomC]
	subjects_anxAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_anxdisord_r > CvSC$AnxC]
	subjects_thoAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_thought_r > CvSC$ThoC]
	subjects_witAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_withdep_r > CvSC$WitC]
	subjects_socAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_sociavoid_r > CvSC$SocC]
	subjects_attAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_attent_r > CvSC$AttC]
	subjects_rulAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_rulebreak_r > CvSC$RulC]
	subjects_aggAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_aggressive_r > CvSC$AggC]
	# get counts of # in clinical threshold for this boot
	nsubjects_pAbove = length(subjects_pAbove)
	nsubjects_intAbove = length(subjects_intAbove)
	nsubjects_extAbove = length(subjects_extAbove)
	nsubjects_somAbove = length(subjects_somAbove)
	nsubjects_anxAbove = length(subjects_anxAbove)
	nsubjects_thoAbove = length(subjects_thoAbove)
	nsubjects_witAbove = length(subjects_witAbove)
	nsubjects_socAbove = length(subjects_socAbove)
	nsubjects_attAbove = length(subjects_attAbove)
	nsubjects_rulAbove = length(subjects_rulAbove)
	nsubjects_aggAbove = length(subjects_aggAbove)
	# sample subjects from this boot n_clinical times for pseudo-clinical group
	subjects_ppAbove = sample(unqSubjs$subjectkey, nsubjects_pAbove)
	subjects_pintAbove = sample(unqSubjs$subjectkey, nsubjects_intAbove)
	subjects_pextAbove = sample(unqSubjs$subjectkey, nsubjects_extAbove)
	subjects_psomAbove = sample(unqSubjs$subjectkey, nsubjects_somAbove)
	subjects_panxAbove = sample(unqSubjs$subjectkey, nsubjects_anxAbove)
	subjects_pthoAbove = sample(unqSubjs$subjectkey, nsubjects_thoAbove)
	subjects_pwitAbove = sample(unqSubjs$subjectkey, nsubjects_witAbove)
	subjects_psocAbove = sample(unqSubjs$subjectkey, nsubjects_socAbove)
	subjects_pattAbove = sample(unqSubjs$subjectkey, nsubjects_attAbove)
	subjects_prulAbove = sample(unqSubjs$subjectkey, nsubjects_rulAbove)
	subjects_paggAbove = sample(unqSubjs$subjectkey, nsubjects_aggAbove)
	# create a vector to assign pseudoclinical to entire subjects
	bootSamp$ppAbove = 0
	bootSamp$pintAbove = 0
	bootSamp$pextAbove = 0
	bootSamp$psomAbove = 0
	bootSamp$panxAbove = 0
	bootSamp$pthoAbove = 0
	bootSamp$pwitAbove = 0
	bootSamp$psocAbove = 0
	bootSamp$pattAbove = 0
	bootSamp$prulAbove = 0
	bootSamp$paggAbove = 0
	#randomly assign n_clinical people to pseudoclinical
	bootSamp$ppAbove[bootSamp$subjectkey %in% subjects_ppAbove] = 1
	bootSamp$pintAbove[bootSamp$subjectkey %in% subjects_pintAbove] = 1
	bootSamp$pextAbove[bootSamp$subjectkey %in% subjects_pextAbove] = 1
	bootSamp$psomAbove[bootSamp$subjectkey %in% subjects_psomAbove] = 1
	bootSamp$panxAbove[bootSamp$subjectkey %in% subjects_panxAbove] = 1
	bootSamp$pthoAbove[bootSamp$subjectkey %in% subjects_pthoAbove] = 1
	bootSamp$pwitAbove[bootSamp$subjectkey %in% subjects_pwitAbove] = 1
	bootSamp$psocAbove[bootSamp$subjectkey %in% subjects_psocAbove] = 1
	bootSamp$pattAbove[bootSamp$subjectkey %in% subjects_pattAbove] = 1
	bootSamp$prulAbove[bootSamp$subjectkey %in% subjects_prulAbove] = 1
	bootSamp$paggAbove[bootSamp$subjectkey %in% subjects_paggAbove] = 1
	# turn into pAboveDf and pBelowDf
	pAboveDF = bootSamp[bootSamp$ppAbove == 1,]
	pBelowDF = bootSamp[bootSamp$ppAbove == 0,]
	intAboveDF = bootSamp[bootSamp$pintAbove == 1,]
	intBelowDF = bootSamp[bootSamp$pintAbove == 0,]
	extAboveDF = bootSamp[bootSamp$pextAbove == 1,]
	extBelowDF = bootSamp[bootSamp$pextAbove == 0,]
	somAboveDF = bootSamp[bootSamp$psomAbove == 1,]
	somBelowDF = bootSamp[bootSamp$psomAbove == 0,]
	anxAboveDF = bootSamp[bootSamp$panxAbove == 1,]
	anxBelowDF = bootSamp[bootSamp$panxAbove == 0,]
	thoAboveDF = bootSamp[bootSamp$pthoAbove == 1,]
	thoBelowDF = bootSamp[bootSamp$pthoAbove == 0,]
	witAboveDF = bootSamp[bootSamp$pwitAbove == 1,]
	witBelowDF = bootSamp[bootSamp$pwitAbove == 0,]
	socAboveDF = bootSamp[bootSamp$psocAbove == 1,]
	socBelowDF = bootSamp[bootSamp$psocAbove == 0,]
	attAboveDF = bootSamp[bootSamp$pattAbove == 1,]
	attBelowDF = bootSamp[bootSamp$pattAbove == 0,]
	rulAboveDF = bootSamp[bootSamp$prulAbove == 1,]
	rulBelowDF = bootSamp[bootSamp$prulAbove == 0,]
	aggAboveDF = bootSamp[bootSamp$paggAbove == 1,]
	aggBelowDF = bootSamp[bootSamp$paggAbove == 0,]
	########## ########## ############### ##########

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
	# now use bootstrap sample to record clin beta and subclin betas for this iteration
	# make real clin vs. subclin split based on thresholds
	pClinDF=bootSamp[bootSamp$cbcl_scr_syn_totprob_r>CvSC$Pc,]
	pSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_totprob_r<CvSC$Pbc,]
	intClinDF=bootSamp[bootSamp$cbcl_scr_syn_internal_r>CvSC$Ic,]
	intSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_internal_r<CvSC$Ibc,]
	extClinDF=bootSamp[bootSamp$cbcl_scr_syn_external_r>CvSC$Ec,]
	extSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_external_r<CvSC$Ebc,]
	somClinDF=bootSamp[bootSamp$cbcl_scr_syn_somatic_r>CvSC$SomC,]
	somSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_somatic_r<CvSC$SomBc,]
	anxClinDF=bootSamp[bootSamp$cbcl_scr_syn_anxdep_r>CvSC$AnxC,]
	anxSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_anxdep_r<CvSC$AnxBc,]
	thoClinDF=bootSamp[bootSamp$cbcl_scr_syn_thought_r>CvSC$ThC,]
	thoSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_thought_r<CvSC$ThBc,]
	witClinDF=bootSamp[bootSamp$cbcl_scr_syn_withdep_r>CvSC$WitC,]
	witSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_withdep_r<CvSC$WitBc,]
	socClinDF=bootSamp[bootSamp$cbcl_scr_syn_social_r>CvSC$SocC,]
	socSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_social_r<CvSC$SocBc,]
	attClinDF=bootSamp[bootSamp$cbcl_scr_syn_attention_r>CvSC$AttC,]
	attSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_attention_r<CvSC$AttBc,]
	rulClinDF=bootSamp[bootSamp$cbcl_scr_syn_rulebreak_r>CvSC$RulC,]
	rulSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_rulebreak_r<CvSC$RulBc,]
	aggClinDF=bootSamp[bootSamp$cbcl_scr_syn_aggressive_r>CvSC$AggC,]
	aggSubclinDF=bootSamp[bootSamp$cbcl_scr_syn_aggressive_r<CvSC$AggBc,]
	# get clin and subclin betas
	pClinBeta[b]=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pClinDF)$coefficients[2]
	pSubclinBeta[b]=lm(g~cbcl_scr_syn_totprob_r+interview_age,data=pSubclinDF)$coefficients[2]
	intClinBeta[b]=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intClinDF)$coefficients[2]
	intSubclinBeta[b]=lm(g~cbcl_scr_syn_internal_r+interview_age,data=intSubclinDF)$coefficients[2]
	extClinBeta[b]=lm(g~cbcl_scr_syn_external_r+interview_age,data=extClinDF)$coefficients[2]
	extSubclinBeta[b]=lm(g~cbcl_scr_syn_external_r+interview_age,data=extSubclinDF)$coefficients[2]
	somClinBeta[b]=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somClinDF)$coefficients[2]
	somSubclinBeta[b]=lm(g~cbcl_scr_syn_somatic_r+interview_age,data=somSubclinDF)$coefficients[2]
	anxClinBeta[b]=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxClinDF)$coefficients[2]
	anxSubclinBeta[b]=lm(g~cbcl_scr_syn_anxdep_r+interview_age,data=anxSubclinDF)$coefficients[2]
	thoClinBeta[b]=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoClinDF)$coefficients[2]
	thoSubclinBeta[b]=lm(g~cbcl_scr_syn_thought_r+interview_age,data=thoSubclinDF)$coefficients[2]
	witClinBeta[b]=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witClinDF)$coefficients[2]
	witSubclinBeta[b]=lm(g~cbcl_scr_syn_withdep_r+interview_age,data=witSubclinDF)$coefficients[2]
	socClinBeta[b]=lm(g~cbcl_scr_syn_social_r+interview_age,data=socClinDF)$coefficients[2]
	socSubclinBeta[b]=lm(g~cbcl_scr_syn_social_r+interview_age,data=socSubclinDF)$coefficients[2]
	attClinBeta[b]=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attClinDF)$coefficients[2]
	attSubclinBeta[b]=lm(g~cbcl_scr_syn_attention_r+interview_age,data=attSubclinDF)$coefficients[2]
	rulClinBeta[b]=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulClinDF)$coefficients[2]
	rulSubclinBeta[b]=lm(g~cbcl_scr_syn_rulebreak_r+interview_age,data=rulSubclinDF)$coefficients[2]
	aggClinBeta[b]=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggClinDF)$coefficients[2]
	aggSubclinBeta[b]=lm(g~cbcl_scr_syn_aggressive_r+interview_age,data=aggSubclinDF)$coefficients[2]
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
# save out all bootstrapped betas in one dataframe
outdf=data.frame(pClinBeta,pSubclinBeta,intClinBeta,intSubclinBeta,extClinBeta,extSubclinBeta,somClinBeta,somSubclinBeta,anxClinBeta,anxSubclinBeta,thoClinBeta,thoSubclinBeta,witClinBeta,witSubclinBeta,socClinBeta,socSubclinBeta,attClinBeta,attSubclinBeta,rulClinBeta,rulSubclinBeta,aggClinBeta,aggSubclinBeta)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gp_CvSC_bootBetas.rds')
print('done with g~p fit bootstrapping!')

