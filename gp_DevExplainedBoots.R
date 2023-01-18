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
# entire subject-ID bootstrapping function from https://stackoverflow.com/questions/11919808/block-bootstrap-from-subject-list
#SubjLevelBoot <- function(x,i) {
#	unlist(lapply(i, function(n) which(x[n] == masterdf$subjectkey)))
#	print
#	do.call("rbind", lapply(i,function(n) subset(masterdf, subjectkey==x[n])))
	#model<-as.formula(f)
	#summary(model)$dev.expl
#}
#create_subject_samples <- function(data, n, subjs) {
#	  data[data$subjectkey %in% sample(subjs, n, replace = TRUE), ]
#}
#create_subject_samples <- function(data, subjects, R) {
#	resample_subjects <- sample(subjects, R, replace = TRUE)
#	data_subset <- data[data$subject %in% resample_subjects, ]
#	resample_indices <- sapply(split(data_subset, data_subset$subject), function(x) sample(1:nrow(x), nrow(x), replace = TRUE))
#	data_subset[unlist(resample_indices),]
#}
# loop over manual bootstrap
for (b in 1:10000){
	print(b)
	# get subjects to include in this bootstrap
	BootSubjs=sample(subjs,numSubjs,replace=T)
	# inefficient but interpretable loop
	# Create an empty dataframe to store the resampled observations
	resampled_df <- data.frame()
	for (j in 1:length(BootSubjs)){
		print(j)
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

