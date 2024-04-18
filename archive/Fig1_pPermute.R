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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age','cbcl_scr_syn_totprob_t')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
### initialize cross-boot vectors
gp_cor_permuted=rep(0,10000)
gp_cor_permuted_sc=rep(0,10000)
gp_cor_permuted_c=rep(0,10000)

# derive clinical vs. subclinical cutoff scores
masterdfP_bc<-masterdf[masterdf$cbcl_scr_syn_totprob_t==65,]
masterdfP_c<-masterdf[masterdf$cbcl_scr_syn_totprob_t==69,]
# borderline clinical and clinical cutoffs
Pbc=mean(masterdfP_bc$cbcl_scr_syn_totprob_r)
Pc=mean(masterdfP_c$cbcl_scr_syn_totprob_r)
# get true cor(g,p) values, controlling for age
modelforresids<-gam(g~s(interview_age),data=masterdf)
masterdf$resids_g<-modelforresids$residuals
corrResult_full=cor.test(masterdf$cbcl_scr_syn_totprob_r,masterdf$resids_g)
# now in subclinical and clinical subsets
clindf<-masterdf[as.numeric(masterdf$cbcl_scr_syn_totprob_r)>Pc,]
sclindf<-masterdf[as.numeric(masterdf$cbcl_scr_syn_totprob_r)<Pbc,]
modelforresids<-gam(g~s(interview_age),data=clindf)
clindf$resids<-modelforresids$residuals
modelforresids<-gam(g~s(interview_age),data=sclindf)
sclindf$resids<-modelforresids$residuals
corrResult_sc=cor.test(sclindf$cbcl_scr_syn_totprob_r,sclindf$resids)
corrResult_c=cor.test(clindf$cbcl_scr_syn_totprob_r,clindf$resids)

# get counts of # in clinical threshold
unqSubjs=masterdf[duplicated(masterdf$subjectkey),]
subjects_pAbove = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_totprob_r > Pc]
nsubjects_pAbove = length(subjects_pAbove)
subjects_pBelow = unqSubjs$subjectkey[unqSubjs$cbcl_scr_syn_totprob_r < Pbc]
nsubjects_pBelow = length(subjects_pBelow)

# initialize a permuted p-factor column
masterdf$perm_p=NULL

# loop over manual bootstrap
for (b in 1:10000){
	print(b)
	# permute p factor values for each subject rather than each observation
	# initialize a matrix to store original values
	original_matrix <- matrix(NA, nrow = numSubjs, ncol = 2)
	for (subj_idx in 1:numSubjs) {
		# get this subjects observations
		p_values <- masterdf[masterdf$subjectkey == subjs[subj_idx], "cbcl_scr_syn_totprob_r"]
		original_matrix[subj_idx, ] <- p_values
	}
	# shuffle p-factor values, subjectwise
	permuted_matrix <- original_matrix[sample(numSubjs), ]
	# Insert permuted values into masterdf
    	for (subj_idx in 1:numSubjs) {
    	    # Get the rows corresponding to this subject in masterdf
    	    subj_rows <- which(masterdf$subjectkey == subjs[subj_idx])
    	    # Insert permuted values into the two rows
    	    masterdf[subj_rows, "perm_p"] <- permuted_matrix[subj_idx, ]
    	}
	# get null correlation
	gp_cor_permuted[b]=cor.test(masterdf$perm_p,masterdf$resids_g)$estimate
	# split and obtain null correlations
	clindf_p<-masterdf[as.numeric(masterdf$perm_p)>Pc,]
	sclindf_p<-masterdf[as.numeric(masterdf$perm_p)<Pbc,]
	gp_cor_permuted_sc[b]=cor.test(sclindf_p$perm_p,sclindf_p$resids_g)$estimate
	gp_cor_permuted_c[b]=cor.test(clindf_p$perm_p,clindf_p$resids_g)$estimate	
}
# insert true correlation values to end of null correlation vectors
gp_cor_permuted[10001]=corrResult_full$estimate
gp_cor_permuted_sc[10001]=corrResult_sc$estimate
gp_cor_permuted_c[10001]=corrResult_c$estimate
# SAVEOUT permuted and true values
outdf=data.frame(gp_cor_permuted,gp_cor_permuted_sc,gp_cor_permuted_c)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots_pfactorCorr.rds')

print('done with g~p fit bootstrapping!')
