library(mgcv)
# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
devExplained_full=rep(0,1000)
devExplained_NoCbcl61=rep(0,1000)
devExplained_NoTot=rep(0,1000)
p_devExplained_full=rep(0,1000)
p_devExplained_NoCbcl61=rep(0,1000)
p_devExplained_NoG=rep(0,1000)
# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','g','subjectkey','interview_age','cbcl_q61_p')]
# loop over manual bootstrap
for (b in 1:10000){
	print(b)
	# get subjects to include in this bootstrap
	BootSubjs=sample(subjs,numSubjs,replace=T)
	### inefficient but interpretable loop
	# Create an empty dataframe to store the resampled observations
	resampled_df <- data.frame()
	for (j in 1:length(BootSubjs)){
		subject_obs <- masterdf[masterdf$subjectkey == BootSubjs[j], ]
		resampled_df <- rbind(resampled_df, subject_obs)
	}
	# bootstrap sample
	bootSamp=resampled_df
	#### gpAge_independent splines
	gpAge_full<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age)+cbcl_q61_p,data=bootSamp)
	devExplained_full[b]<-summary(gpAge_full)$dev.expl
	# fit version with cbcl 61
	gpAge_nocbcl61<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	devExplained_NoCbcl61[b]<-summary(gpAge_nocbcl61)$dev.expl
	# fit version with just cbcl61
	gpAge_noTotProbs<-bam(g~cbcl_q61_p+s(interview_age),data=bootSamp)
	devExplained_NoTot[b]<-summary(gpAge_noTotProbs)$dev.expl
	### version with symptom count as response variable
	pgAge_full<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+cbcl_q61_p,data=bootSamp,family=nb())
	p_devExplained_full[b]<-summary(pgAge_full)$dev.expl
	# fit version with cbcl 61
	pgAge_nocbcl61<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age),data=bootSamp,family=nb())
	p_devExplained_NoCbcl61[b]<-summary(pgAge_nocbcl61)$dev.expl
	# fit version with just cbcl61
	pgAge_noG<-bam(cbcl_scr_syn_totprob_r~cbcl_q61_p+s(interview_age),data=bootSamp,family=nb())
	p_devExplained_NoG[b]<-summary(pgAge_noG)$dev.expl
}
# saveout df of dev explained for plotting
outdf=data.frame(devExplained_full,devExplained_NoCbcl61,devExplained_NoTot,p_devExplained_full,p_devExplained_NoCbcl61,p_devExplained_NoG)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/DevExplained.rds')

