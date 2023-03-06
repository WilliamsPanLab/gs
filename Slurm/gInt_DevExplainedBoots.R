library(mgcv)
library(testit)
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
masterdf=masterdf[,c('cbcl_scr_syn_internal_r','ple_died_y','ple_injured_y','ple_crime_y','ple_friend_y','ple_friend_injur_y','ple_arrest_y','ple_friend_died_y','ple_mh_y','ple_sib_y','ple_victim_y','ple_separ_y','ple_law_y','ple_school_y','ple_move_y','ple_jail_y','ple_step_y','ple_new_job_y','ple_new_sib_y','g','subjectkey','interview_age','Grades','parentPcount','income','parental_education','sex','race_ethnicity')]
# loop over manual bootstrap
for (b in 1:6000){
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
	fullModel<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_full[b]=fullModel$deviance/fullModel$null.deviance
	### died
	Model_n_died<-bam(cbcl_scr_syn_internal_r~ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_died[b]=Model_n_died$deviance/Model_n_died$null.deviance
	### injured
	Model_n_injured<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_injured[b]=Model_n_injured$deviance/Model_n_injured$null.deviance
	### crime
	Model_n_crime<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_crime[b]=Model_n_crime$deviance/Model_n_crime$null.deviance
	### friend
	Model_n_friend<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend[b]=Model_n_friend$deviance/Model_n_friend$null.deviance
	### friend_injured
	Model_n_friend_injured<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend_injur[b]=Model_n_friend_injured$deviance/Model_n_friend_injured$null.deviance
	### arrest
	Model_n_arrest<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_arrest[b]=Model_n_arrest$deviance/Model_n_arrest$null.deviance
	### friend_died
	Model_n_friend_died<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_friend_died[b]=Model_n_friend_died$deviance/Model_n_friend_died$null.deviance
	### mental_health
	Model_n_mh<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_mh[b]=Model_n_mh$deviance/Model_n_mh$null.deviance
	### sibling
	Model_n_sib<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_sib[b]=Model_n_sib$deviance/Model_n_sib$null.deviance
	### victim
	Model_n_victim<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_victim[b]=Model_n_victim$deviance/Model_n_victim$null.deviance
	### separation
	Model_n_separ<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_separ[b]=Model_n_separ$deviance/Model_n_separ$null.deviance
	### law
	Model_n_law<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_law[b]=Model_n_law$deviance/Model_n_law$null.deviance
	### school
	Model_n_school<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_school[b]=Model_n_school$deviance/Model_n_school$null.deviance
	### move
	Model_n_move<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_move[b]=Model_n_move$deviance/Model_n_move$null.deviance
	### jail
	Model_n_jail<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_jail[b]=Model_n_jail$deviance/Model_n_jail$null.deviance
	### step
	Model_n_step<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_step[b]=Model_n_step$deviance/Model_n_step$null.deviance
	### new job
	Model_n_job<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_new_job[b]=Model_n_job$deviance/Model_n_job$null.deviance
	### new sib
	Model_n_new_sib<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_new_sib[b]=Model_n_new_sib$deviance/Model_n_new_sib$null.deviance
	### g
	Model_n_g<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_g[b]=Model_n_g$deviance/Model_n_g$null.deviance
	### age
	Model_n_age<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_interview_age[b]=Model_n_age$deviance/Model_n_age$null.deviance
	### grades
	Model_n_grades<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_Grades[b]=Model_n_grades$deviance/Model_n_grades$null.deviance
	### parentPcount
	Model_n_parentP<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(income,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_parentPcount[b]=Model_n_parentP$deviance/Model_n_parentP$null.deviance
	### income
	Model_n_income<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(parental_education,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_income[b]=Model_n_income$deviance/Model_n_income$null.deviance
	### parental education
	Model_n_parental_edu<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+sex+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_parental_education[b]=Model_n_parental_edu$deviance/Model_n_parental_edu$null.deviance
	### sex
	Model_n_sex<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+race_ethnicity,data=bootSamp,family=nb())
	devExplained_n_sex[b]=Model_n_sex$deviance/Model_n_sex$null.deviance
	### race
	Model_n_race<-bam(cbcl_scr_syn_internal_r~ple_died_y+ple_injured_y+ple_crime_y+ple_friend_y+ple_friend_injur_y+ple_arrest_y+ple_friend_died_y+ple_mh_y+ple_sib_y+ple_victim_y+ple_separ_y+ple_law_y+ple_school_y+ple_move_y+ple_jail_y+ple_step_y+ple_new_job_y+ple_new_sib_y+s(g,k=4)+s(interview_age,k=4)+s(Grades,k=4)+s(parentPcount,k=4)+s(income,k=4)+s(parental_education,k=4)+sex,data=bootSamp,family=nb())
	devExplained_n_race[b]=Model_n_race$deviance/Model_n_race$null.deviance
	### verify same number of terms in each model (except first because it is full/+1)
	assert((length(fullModel$var.summary)-1)==length(Model_n_died$var.summary))
	assert(length(Model_n_died$var.summary)==length(Model_n_injured$var.summary))
	assert(length(Model_n_injured$var.summary)==length(Model_n_crime$var.summary))
	assert(length(Model_n_crime$var.summary)==length(Model_n_friend$var.summary))
	assert(length(Model_n_friend$var.summary)==length(Model_n_friend_injured$var.summary))
	assert(length(Model_n_friend_injured$var.summary)==length(Model_n_arrest$var.summary))
	assert(length(Model_n_arrest$var.summary)==length(Model_n_friend_died$var.summary))
	assert(length(Model_n_friend_died$var.summary)==length(Model_n_mh$var.summary))
	assert(length(Model_n_mh$var.summary)==length(Model_n_sib$var.summary))
	assert(length(Model_n_sib$var.summary)==length(Model_n_victim$var.summary))
	assert(length(Model_n_victim$var.summary)==length(Model_n_separ$var.summary))
	assert(length(Model_n_separ$var.summary)==length(Model_n_law$var.summary))
	assert(length(Model_n_law$var.summary)==length(Model_n_school$var.summary))
	assert(length(Model_n_school$var.summary)==length(Model_n_move$var.summary))
	assert(length(Model_n_move$var.summary)==length(Model_n_jail$var.summary))
	assert(length(Model_n_jail$var.summary)==length(Model_n_step$var.summary))
	assert(length(Model_n_step$var.summary)==length(Model_n_job$var.summary))
	assert(length(Model_n_job$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_new_sib$var.summary)==length(Model_n_parental_edu$var.summary))
	assert(length(Model_n_g$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_age$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_grades$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_parentP$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_income$var.summary)==length(Model_n_new_sib$var.summary))
	assert(length(Model_n_parental_edu$var.summary)==length(Model_n_sex$var.summary))
	assert(length(Model_n_sex$var.summary)==length(Model_n_race$var.summary))
}
# saveout df of dev explained for plotting
outdf=data.frame(devExplained_full,devExplained_n_died,devExplained_n_injured,devExplained_n_crime,devExplained_n_friend,devExplained_n_friend_injur,devExplained_n_arrest,devExplained_n_friend_died,devExplained_n_mh,devExplained_n_sib,devExplained_n_victim,devExplained_n_separ,devExplained_n_law,devExplained_n_school,devExplained_n_move,devExplained_n_jail,devExplained_n_step,devExplained_n_new_job,devExplained_n_new_sib,devExplained_n_g,devExplained_n_interview_age,devExplained_n_Grades,devExplained_n_parentPcount,devExplained_n_income,devExplained_n_parental_education,devExplained_n_sex,devExplained_n_race)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/DevExplained_Int.rds')

