library(mgcv)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf_WithGrades.rds')
devExplained_full=rep(0,10000)
devExplained_NoGrades=rep(0,10000)
devExplained_NoTot=rep(0,10000)
e_devExplained_full=rep(0,10000)
e_devExplained_NoGrades=rep(0,10000)
e_devExplained_NoG=rep(0,10000)
grades_devExplained_full=rep(0,10000)
grades_devExplained_noG=rep(0,10000)
grades_devExplained_noTot=rep(0,10000)
# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('cbcl_scr_syn_external_r','g','subjectkey','interview_age','Grades')]
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
	#### full model
	gpAge_full<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age)+s(Grades,k=4),data=bootSamp)
	devExplained_full[b]<-summary(gpAge_full)$dev.expl
	# fit version without Grades
	gpAge_noGrades<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	devExplained_NoGrades[b]<-summary(gpAge_noGrades)$dev.expl
	# fit version with just Grades
	gpAge_noTotProbs<-bam(g~s(Grades,k=4)+s(interview_age),data=bootSamp)
	devExplained_NoTot[b]<-summary(gpAge_noTotProbs)$dev.expl
	### version with symptom count as response variable
	pgAge_full<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age)+s(Grades,k=4),data=bootSamp,family=nb())
	e_devExplained_full[b]<-summary(pgAge_full)$dev.expl
        # fit version without Grades
        pgAge_noGrades<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age),data=bootSamp,family=nb())
        e_devExplained_NoGrades[b]<-summary(pgAge_noGrades)$dev.expl
	# fit version with just Grades
	pgAge_noG<-bam(cbcl_scr_syn_external_r~s(Grades,k=4)+s(interview_age),data=bootSamp,family=nb())
	e_devExplained_NoG[b]<-summary(pgAge_noG)$dev.expl
	# variance explained in grades
	spgAge_full<-bam(Grades~s(g)+s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	grades_devExplained_full[b]<-summary(spgAge_full)$dev.expl
	spgAge_noG<-bam(Grades~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	grades_devExplained_noG[b]<-summary(spgAge_noG)$dev.expl
	spgAge_noTot<-bam(Grades~s(g)+s(interview_age),data=bootSamp)
	grades_devExplained_noTot[b]<-summary(spgAge_noTot)$dev.expl
}

# saveout df of dev explained for plotting
outdf=data.frame(devExplained_full,devExplained_NoGrades,devExplained_NoTot,e_devExplained_full,e_devExplained_NoGrades,e_devExplained_NoG,grades_devExplained_full,grades_devExplained_noG,grades_devExplained_noTot)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/DevExplainedExt.rds')

