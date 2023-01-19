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
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','g','subjectkey','interview_age','cbcl_q61_p')]
# get length of df for later
lenDF=dim(masterdf)[1]

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
	##############################
	######## I PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#### g as response variable
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age),data=masterdf)
	# use model fit for saving
	forSpline<-predict(pgAge, newdata = masterdf)
	# save out forspline to a matrix
	predCIs[b,]=forSpline
	######## II FORMALLY TEST FOR NON-LINEARITY
	#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
	pgAge<-bam(g~cbcl_scr_syn_totprob_r+s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age)+ti(cbcl_scr_syn_totprob_r,interview_age),data=masterdf)
	linBoots[b]=summary(pgAge)$s.pv[1]
	####### III TEST INTERACTION
	### interaction test
	gpAge_intrxn<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+g*interview_age,data=masterdf,family=nb())
	intrxnBootsStat[b]=summary(gpAge_intrxn)$p.coeff[4]
	intrxnBootsP[bb=summary(gpAge_intrxn)$s.pv[4]
	###### FULL VS REDUCED ANOVA: P AND DR2
	no_g_Gam<-bam(cbcl_scr_syn_totprob_r~s(interview_age),data=masterdf,family=nb())
	no_g_Sum<-summary(no_g_Gam)
	# g-included model for measuring difference
	gGam<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age),data=masterdf,family=nb())
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
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpBoots.rds')
outdf=data.frame(predCIs)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpPredictedBoots.rds')
print('done with g~p fit bootstrapping!')
