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

# add derivatives plotting
source(here::here('/oak/stanford/groups/leanew1/users/apines/scripts/gp', 'get_derivs_and_plot.R'))

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDF_wAdultp.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('parentPcount','g','subjectkey','interview_age','cbcl_q61_p')]
# get length of df for later
lenDF=dim(masterdf)[1]
masterdf$interview_age=as.numeric(masterdf$interview_age)
masterdf$subjectkey=as.factor(masterdf$subjectkey)
# initialize cross-boot vectors

# linear? 1xbootlength
linBoots=rep(0,10000)
# predicted for CIs? bootlengthxpredictlength
# NOTE THIS IS PREDICTED ON MASTERDF ORDER: NOT FROM LOW-TO-HIGH ON RESPONSE VARIABLE
predCIs=matrix(0,10000,lenDF)
# interaction? 1xbootlength
intrxnBootsStat=rep(0,10000)
intrxnBootsP=rep(0,10000)
# DR2 of relationship? 1xbootlength
dr2Boots=rep(0,10000)
# p of relationship? 1xbootlength
pBoots=rep(0,10000)
# loop over manual bootstrap
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
	##############################
	######## I PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
	pgAge<-bam(g~s(parentPcount)+s(interview_age),data=bootSamp)
	# use model fit for saving
	forSpline<-predict(pgAge, newdata = masterdf)
	# save out forspline to a matrix
	predCIs[b,]=forSpline
	######## II FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(g~parentPcount+s(parentPcount,m=c(2,0))+s(interview_age)+ti(parentPcount,interview_age),data=bootSamp)
	linBoots[b]=summary(pgAge)$s.pv[1]
	####### III TEST INTERACTION
	### interaction test
	gpAge_intrxn<-bam(parentPcount~s(g)+s(interview_age)+g*interview_age,data=bootSamp,family=nb())
	intrxnBootsStat[b]=summary(gpAge_intrxn)$p.coeff[4]
	intrxnBootsP[b]=summary(gpAge_intrxn)$s.pv[4]
	###### FULL VS REDUCED ANOVA: P AND DR2
	no_g_Gam<-bam(parentPcount~s(interview_age),data=bootSamp,family=nb())
	no_g_Sum<-summary(no_g_Gam)
	# g-included model for measuring difference
	gGam<-bam(parentPcount~s(g)+s(interview_age),data=bootSamp,family=nb())
	gSum<-summary(gGam)
	dif<-gSum$r.sq-no_g_Sum$r.sq
	# insert into vec
	dr2Boots[b]=dif
	# test of dif with anova.gam
	anovaRes<-anova.gam(no_g_Gam,gGam,test='Chisq')
	anovaP<-anovaRes$`Pr(>Chi)`
	anovaP2<-unlist(anovaP)
	# insert into vec
	pBoots[b]=anovaP2[2]
}
# SAVEOUT
outdf=data.frame(linBoots,intrxnBootsStat,intrxnBootsP,dr2Boots,pBoots)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gparentPBoots.rds')
outdf=data.frame(predCIs)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gparentPPredictedBoots.rds')
print('done with g~parentP fit bootstrapping!')
