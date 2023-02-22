library(mgcv)
# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/OutDfxc.rds')

# dev explained for all models
devExplained_full=rep(0,10000)
devExplained_n_died=rep(0,10000)
devExplained_n_injured=rep(0,10000)
devExplained_n_crime=rep(0,10000)
devExplained_n_friend=rep(0,10000)
devExplained_n_friend_injur=rep(0,10000)
devExplained_n_arrest=rep(0,10000)
devExplained_n_friend_died=rep(0,10000)
devExplained_n_mh=rep(0,10000)
devExplained_n_sib=rep(0,10000)
devExplained_n_victim=rep(0,10000)
devExplained_n_separ=rep(0,10000)
devExplained_n_law=rep(0,10000)
devExplained_n_school=rep(0,10000)
devExplained_n_move=rep(0,10000)
devExplained_n_jail=rep(0,10000)
devExplained_n_step=rep(0,10000)
devExplained_n_new_job=rep(0,10000)
devExplained_n_new_sib=rep(0,10000)
devExplained_n_g=rep(0,10000)
devExplained_n_interview_age=rep(0,10000)
devExplained_n_Grades=rep(0,10000)
devExplained_n_parentPcount=rep(0,10000)
devExplained_n_income=rep(0,10000)
devExplained_n_parental_education=rep(0,10000)
devExplained_n_sex=rep(0,10000)
devExplained_n_race=rep(0,10000)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','ple_died_y','ple_injured_y','ple_crime_y','ple_friend_y','ple_friend_injur_y','ple_arrest_y','ple_friend_died_y','ple_mh_y','ple_sib_y','ple_victim_y','ple_separ_y','ple_law_y','ple_school_y','ple_move_y','ple_jail_y','ple_step_y','ple_new_job_y','ple_new_sib_y','g','subjectkey','interview_age','Grades','parentPcount','income','parental_education','sex','race_ethnicity')]
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
	# fit models
	### full
	fullModel<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_full[b]=fullModel$deviance/fullModel$null.deviance
	### died
	Model_n_died<-bam(cbcl_scr_syn_totprob_r~ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_died[b]=Model_n_died$deviance/Model_n_died$null.deviance
	### injured
	Model_n_injured<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_injured[b]=Model_n_injured$deviance/Model_n_injured$null.deviance
	### crime
	Model_n_crime<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_crime[b]=Model_n_crime$deviance/Model_n_crime$null.deviance
	### friend
	Model_n_friend<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend[b]=Model_n_friend$deviance/Model_n_friend$null.deviance
	### friend_injured
	Model_n_friend_injured<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend_injured[b]=Model_n_friend_injured$deviance/Model_n_friend_injured$null.deviance
	### arrest
	Model_n_arrest<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_arrest[b]=Model_n_arrest$deviance/Model_n_arrest$null.deviance
	### friend_died
	Model_n_friend_died<-bam(cbcl_scr_syn_totprob_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_arrest_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g)+s(interview_age)+s(Grades)+s(parentPcount)+s(income)+s(parental_education)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend_died[b]=Model_n_friend_died$deviance/Model_n_friend_died$null.deviance
	### mental_health
	##### LEFT OFF HERE - CONVERT SPLINE VARIABLES TO NUMERIC AND PLES TO FACTORS
	

	#### gpAge_independent splines
	gpAge_full<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age)+Grades,data=bootSamp)
	devExplained_full[b]<-summary(gpAge_full)$dev.expl
	# fit version without Grades
	gpAge_noGrades<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	devExplained_NoGrades[b]<-summary(gpAge_noGrades)$dev.expl
	# fit version with just Grades
	gpAge_noTotProbs<-bam(g~Grades+s(interview_age),data=bootSamp)
	devExplained_NoTot[b]<-summary(gpAge_noTotProbs)$dev.expl
	### version with symptom count as response variable
	pgAge_full<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+Grades,data=bootSamp,family=nb())
	p_devExplained_full[b]<-summary(pgAge_full)$dev.expl
	# fit version without Grades
	pgAge_noGrades<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age),data=bootSamp,family=nb())
	p_devExplained_NoGrades[b]<-summary(pgAge_noGrades)$dev.expl
	# fit version with just Grades
	pgAge_noG<-bam(cbcl_scr_syn_totprob_r~Grades+s(interview_age),data=bootSamp,family=nb())
	p_devExplained_NoG[b]<-summary(pgAge_noG)$dev.expl
	# variance explained in grades
	bootSamp$Grades<-as.numeric(bootSamp$Grades)
	spgAge_full<-bam(Grades~s(g)+s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	grades_devExplained_full[b]<-summary(spgAge_full)$dev.expl
	spgAge_noG<-bam(Grades~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=bootSamp)
	grades_devExplained_noG[b]<-summary(spgAge_noG)$dev.expl
	spgAge_noTot<-bam(Grades~s(g)+s(interview_age),data=bootSamp)
	grades_devExplained_noTot[b]<-summary(spgAge_noTot)$dev.expl
}
# saveout df of dev explained for plotting
outdf=data.frame(devExplained_full,devExplained_NoGrades,devExplained_NoTot,p_devExplained_full,p_devExplained_NoGrades,p_devExplained_NoG,grades_devExplained_full,grades_devExplained_noG,grades_devExplained_noTot)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/DevExplained.rds')

