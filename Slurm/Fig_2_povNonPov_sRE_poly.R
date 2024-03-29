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
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','site','g','subjectkey','interview_age','Pov_v2')]
# get length of df for later
lenDF=dim(masterdf)[1]
# set site to facotr
masterdf$site=as.factor(masterdf$site)
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
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
# predicted derivatives: set to maximum value for ncol
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
iMaxVal=max(masterdf$cbcl_scr_syn_internal_r)
emaxVal=max(masterdf$cbcl_scr_syn_external_r)
somMaxVal=max(masterdf$cbcl_scr_syn_somatic_r)
anxMaxVal=max(masterdf$cbcl_scr_syn_anxdep_r)
thoMaxVal=max(masterdf$cbcl_scr_syn_thought_r)
witMaxVal=max(masterdf$cbcl_scr_syn_withdep_r)
socMaxVal=max(masterdf$cbcl_scr_syn_social_r)
attMaxVal=max(masterdf$cbcl_scr_syn_attention_r)
rulMaxVal=max(masterdf$cbcl_scr_syn_rulebreak_r)
aggMaxVal=max(masterdf$cbcl_scr_syn_aggressive_r)
# predicted values: set to maximum value for ncol
pFit=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFit=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFit=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFit=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFit=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFit=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFit=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFit=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFit=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFit=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFit=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
# initialize version of each for poverty and nonpoverty groups
pFitPov=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitPov=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitPov=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitPov=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitPov=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitPov=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitPov=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitPov=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitPov=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitPov=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitPov=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
pFitNonPov=matrix(0,nrow=10000,ncol=(pMaxVal)+1)
intFitNonPov=matrix(0,nrow=10000,ncol=(iMaxVal)+1)
extFitNonPov=matrix(0,nrow=10000,ncol=(emaxVal)+1)
somFitNonPov=matrix(0,nrow=10000,ncol=(somMaxVal)+1)
anxFitNonPov=matrix(0,nrow=10000,ncol=(anxMaxVal)+1)
thoFitNonPov=matrix(0,nrow=10000,ncol=(thoMaxVal)+1)
witFitNonPov=matrix(0,nrow=10000,ncol=(witMaxVal)+1)
socFitNonPov=matrix(0,nrow=10000,ncol=(socMaxVal)+1)
attFitNonPov=matrix(0,nrow=10000,ncol=(attMaxVal)+1)
rulFitNonPov=matrix(0,nrow=10000,ncol=(rulMaxVal)+1)
aggFitNonPov=matrix(0,nrow=10000,ncol=(aggMaxVal)+1)
# initialize difference in AIC for each cbcl model with and without poverty interaction
pDiff=rep(0,10000)
intDiff=rep(0,10000)
extDiff=rep(0,10000)
somDiff=rep(0,10000)
anxDiff=rep(0,10000)
thoDiff=rep(0,10000)
witDiff=rep(0,10000)
socDiff=rep(0,10000)
attDiff=rep(0,10000)
rulDiff=rep(0,10000)
aggDiff=rep(0,10000)
# difference in AIC for each cbcl model with psuedo poverty interaction
pDiffPseudo=rep(0,10000)
intDiffPseudo=rep(0,10000)
extDiffPseudo=rep(0,10000)
somDiffPseudo=rep(0,10000)
anxDiffPseudo=rep(0,10000)
thoDiffPseudo=rep(0,10000)
witDiffPseudo=rep(0,10000)
socDiffPseudo=rep(0,10000)
attDiffPseudo=rep(0,10000)
rulDiffPseudo=rep(0,10000)
aggDiffPseudo=rep(0,10000)
# initialize difference in adjusted r^2 for each cbcl model with and without poverty interaction
pDiffAdj=rep(0,10000)
intDiffAdj=rep(0,10000)
extDiffAdj=rep(0,10000)
somDiffAdj=rep(0,10000)
anxDiffAdj=rep(0,10000)
thoDiffAdj=rep(0,10000)
witDiffAdj=rep(0,10000)
socDiffAdj=rep(0,10000)
attDiffAdj=rep(0,10000)
rulDiffAdj=rep(0,10000)
aggDiffAdj=rep(0,10000)
# initialize difference in adjusted r^2 for each cbcl model with psuedo poverty interaction
pDiffAdjPseudo=rep(0,10000)
intDiffAdjPseudo=rep(0,10000)
extDiffAdjPseudo=rep(0,10000)
somDiffAdjPseudo=rep(0,10000)
anxDiffAdjPseudo=rep(0,10000)
thoDiffAdjPseudo=rep(0,10000)
witDiffAdjPseudo=rep(0,10000)
socDiffAdjPseudo=rep(0,10000)
attDiffAdjPseudo=rep(0,10000)
rulDiffAdjPseudo=rep(0,10000)
aggDiffAdjPseudo=rep(0,10000)
# initialize difference in adjusted r^2 for each asr model with and without poverty interaction
asrPDiffAdj=rep(0,10000)
asrintDiffAdj=rep(0,10000)
asrextDiffAdj=rep(0,10000)
asrsomDiffAdj=rep(0,10000)
asranxDiffAdj=rep(0,10000)
asrthoDiffAdj=rep(0,10000)
asrwitDiffAdj=rep(0,10000)
asrattDiffAdj=rep(0,10000)
asrrulDiffAdj=rep(0,10000)
asraggDiffAdj=rep(0,10000)
# interaction vs. no interaction difference in MAE: cbcl
pDiffMAE=rep(0,10000)
intDiffMAE=rep(0,10000)
extDiffMAE=rep(0,10000)
somDiffMAE=rep(0,10000)
anxDiffMAE=rep(0,10000)
thoDiffMAE=rep(0,10000)
witDiffMAE=rep(0,10000)
socDiffMAE=rep(0,10000)
attDiffMAE=rep(0,10000)
rulDiffMAE=rep(0,10000)
aggDiffMAE=rep(0,10000)
# maximum value in each iteration
pMax=rep(0,10000)
intMax=rep(0,10000)
extMax=rep(0,10000)
somMax=rep(0,10000)
anxMax=rep(0,10000)
thoMax=rep(0,10000)
witMax=rep(0,10000)
socMax=rep(0,10000)
attMax=rep(0,10000)
rulMax=rep(0,10000)
aggMax=rep(0,10000)
# finally, add pov classification (v2)
masterdf$poverty=0
# poverty now defined in sample construction
masterdf$poverty[masterdf$Pov_v2==1]=1
masterdf$poverty=as.factor(masterdf$poverty)
# loop over manual bootstrap
for (b in 1:1000){
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
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
        # last pre-step: need to create NULL variables. Use count of Poverty as count of membership to psuedo-groups, but randomly distribute membership
	# get unique instances of subjects, note duplcated and !duplicated are same because 2x instances of each subj
	unqSubjs=bootSamp[duplicated(bootSamp$subjectkey),]
        # get count poverty in this boot
	subjects_with_poverty = unqSubjs$subjectkey[unqSubjs$poverty == 1]
	subjects_with_pseudopoverty = sample(unqSubjs$subjectkey, length(subjects_with_poverty))
	# create a vector to assign pseudopoverty to entire subjects
	bootSamp$pseudopoverty = 0
        # randomly assign Povcount people to pseudopoverty
	bootSamp$pseudopoverty[bootSamp$subjectkey %in% subjects_with_pseudopoverty] = 1
	bootSamp$pseudopoverty=as.factor(bootSamp$pseudopoverty)
        # make a df that is just the same number of kids in poverty, but actually from the non-poverty group. The point is to see if we can recover g~p slope in a an equivalent number of non-poverty kids
        #bootSamp$pseudopoverty2=0
        #bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
	##############################
	# get max for this iteration
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
	bimax=max(bootSamp$cbcl_scr_syn_internal_r)
	bemax=max(bootSamp$cbcl_scr_syn_external_r)
	bsommax=max(bootSamp$cbcl_scr_syn_somatic_r)
	banxmax=max(bootSamp$cbcl_scr_syn_anxdep_r)
	bthomax=max(bootSamp$cbcl_scr_syn_thought_r)
	bwitmax=max(bootSamp$cbcl_scr_syn_withdep_r)
	bsocmax=max(bootSamp$cbcl_scr_syn_social_r)
	battmax=max(bootSamp$cbcl_scr_syn_attention_r)
	brulmax=max(bootSamp$cbcl_scr_syn_rulebreak_r)
	baggmax=max(bootSamp$cbcl_scr_syn_aggressive_r)
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(0:bpmax)
	eachIntcount=seq(0:bimax)
	eachExtcount=seq(0:bemax)
	eachSomcount=seq(0:bsommax)
	eachAnxcount=seq(0:banxmax)
	eachThocount=seq(0:bthomax)
	eachWitcount=seq(0:bwitmax)
	eachSoccount=seq(0:bsocmax)
	eachAttcount=seq(0:battmax)
	eachRulcount=seq(0:brulmax)
	eachAggcount=seq(0:baggmax)
	# set age to to median for predict df, also set child symptom score to median for predict df
	#####################################
        # set age to to median for predict df, +1 to max values because 0-max sequence is 1 > max in length.
        predictDFp=data.frame(eachPcount,rep(median(masterdf$interview_age),bpmax+1),rep('site13',bpmax+1))
        predictDFint=data.frame(eachIntcount,rep(median(masterdf$interview_age),bimax+1),rep('site13',bimax+1))
        predictDFext=data.frame(eachExtcount,rep(median(masterdf$interview_age),bemax+1),rep('site13',bemax+1))
        predictDFsom=data.frame(eachSomcount,rep(median(masterdf$interview_age),bsommax+1),rep('site13',bsommax+1))
        predictDFanx=data.frame(eachAnxcount,rep(median(masterdf$interview_age),banxmax+1),rep('site13',banxmax+1))
        predictDFtho=data.frame(eachThocount,rep(median(masterdf$interview_age),bthomax+1),rep('site13',bthomax+1))
        predictDFwit=data.frame(eachWitcount,rep(median(masterdf$interview_age),bwitmax+1),rep('site13',bwitmax+1))
        predictDFsoc=data.frame(eachSoccount,rep(median(masterdf$interview_age),bsocmax+1),rep('site13',bsocmax+1))
        predictDFatt=data.frame(eachAttcount,rep(median(masterdf$interview_age),battmax+1),rep('site13',battmax+1))
        predictDFrul=data.frame(eachRulcount,rep(median(masterdf$interview_age),brulmax+1),rep('site13',brulmax+1))
        predictDFagg=data.frame(eachAggcount,rep(median(masterdf$interview_age),baggmax+1),rep('site13',baggmax+1))
        # set colnames so predict can work
        colnames(predictDFp)=c('cbcl_scr_syn_totprob_r','interview_age','site')
        colnames(predictDFint)=c('cbcl_scr_syn_internal_r','interview_age','site')
        colnames(predictDFext)=c('cbcl_scr_syn_external_r','interview_age','site')
        colnames(predictDFsom)=c('cbcl_scr_syn_somatic_r','interview_age','site')
        colnames(predictDFanx)=c('cbcl_scr_syn_anxdep_r','interview_age','site')
        colnames(predictDFtho)=c('cbcl_scr_syn_thought_r','interview_age','site')
        colnames(predictDFwit)=c('cbcl_scr_syn_withdep_r','interview_age','site')
        colnames(predictDFsoc)=c('cbcl_scr_syn_social_r','interview_age','site')
        colnames(predictDFatt)=c('cbcl_scr_syn_attention_r','interview_age','site')
        colnames(predictDFrul)=c('cbcl_scr_syn_rulebreak_r','interview_age','site')
        colnames(predictDFagg)=c('cbcl_scr_syn_aggressive_r','interview_age','site')	
	####### II PREDICT POVERTY INTERACTIONS #######
	# fit models with standalone poverty term
	# cbcl
	pgAge_pov=bam(g~poly(cbcl_scr_syn_totprob_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	intgAge_pov=bam(g~poly(cbcl_scr_syn_internal_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	extgAge_pov=bam(g~poly(cbcl_scr_syn_external_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	somgAge_pov=bam(g~poly(cbcl_scr_syn_somatic_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	anxgAge_pov=bam(g~poly(cbcl_scr_syn_anxdep_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	thogAge_pov=bam(g~poly(cbcl_scr_syn_thought_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	witgAge_pov=bam(g~poly(cbcl_scr_syn_withdep_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	socgAge_pov=bam(g~poly(cbcl_scr_syn_social_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	attgAge_pov=bam(g~poly(cbcl_scr_syn_attention_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	rulgAge_pov=bam(g~poly(cbcl_scr_syn_rulebreak_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	agggAge_pov=bam(g~poly(cbcl_scr_syn_aggressive_r,2)*poverty+poly(interview_age,2)+poverty+s(site,bs="re"),data=bootSamp)
	# make new predict dataframes with poverty variable
	# cbcl
	predictDFp$poverty=1
	predictDFint$poverty=1
	predictDFext$poverty=1
	predictDFsom$poverty=1
	predictDFanx$poverty=1
	predictDFtho$poverty=1	
	predictDFwit$poverty=1
	predictDFsoc$poverty=1
	predictDFatt$poverty=1
	predictDFrul$poverty=1
	predictDFagg$poverty=1
	# get poverty fits
	# cbcl
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	# cbcl
	pFitPov[b,1:(bpmax+1)]=forFitp
	intFitPov[b,1:(bimax+1)]=forFitint
	extFitPov[b,1:(bemax+1)]=forFitext
	somFitPov[b,1:(bsommax+1)]=forFitsom
	anxFitPov[b,1:(banxmax+1)]=forFitanx
	thoFitPov[b,1:(bthomax+1)]=forFittho
	witFitPov[b,1:(bwitmax+1)]=forFitwit
	socFitPov[b,1:(bsocmax+1)]=forFitsoc
	attFitPov[b,1:(battmax+1)]=forFitatt
	rulFitPov[b,1:(brulmax+1)]=forFitrul
	aggFitPov[b,1:(baggmax+1)]=forFitagg
	# update predict df to nonpoverty, get nonpoverty fits
	# cbcl
	predictDFp$poverty=0
	predictDFint$poverty=0
	predictDFext$poverty=0
	predictDFsom$poverty=0
	predictDFanx$poverty=0
	predictDFtho$poverty=0
	predictDFwit$poverty=0
	predictDFsoc$poverty=0
	predictDFatt$poverty=0
	predictDFrul$poverty=0
	predictDFagg$poverty=0
	# get nonpoverty fits
	# cbcl
	forFitp=predict(pgAge_pov,predictDFp)
	forFitint=predict(intgAge_pov,predictDFint)
	forFitext=predict(extgAge_pov,predictDFext)
	forFitsom=predict(somgAge_pov,predictDFsom)
	forFitanx=predict(anxgAge_pov,predictDFanx)
	forFittho=predict(thogAge_pov,predictDFtho)
	forFitwit=predict(witgAge_pov,predictDFwit)
	forFitsoc=predict(socgAge_pov,predictDFsoc)
	forFitatt=predict(attgAge_pov,predictDFatt)
	forFitrul=predict(rulgAge_pov,predictDFrul)
	forFitagg=predict(agggAge_pov,predictDFagg)
	# print out fit
	# cbcl
	pFitNonPov[b,1:(bpmax+1)]=forFitp
	intFitNonPov[b,1:(bimax+1)]=forFitint
	extFitNonPov[b,1:(bemax+1)]=forFitext
	somFitNonPov[b,1:(bsommax+1)]=forFitsom
	anxFitNonPov[b,1:(banxmax+1)]=forFitanx
	thoFitNonPov[b,1:(bthomax+1)]=forFittho
	witFitNonPov[b,1:(bwitmax+1)]=forFitwit
	socFitNonPov[b,1:(bsocmax+1)]=forFitsoc
	attFitNonPov[b,1:(battmax+1)]=forFitatt
	rulFitNonPov[b,1:(brulmax+1)]=forFitrul
}
# SAVEOUT
# save out poverty and nonpoverty fits - cbcl
outdf=data.frame(pFitPov,pFitNonPov,intFitPov,intFitNonPov,extFitPov,extFitNonPov,somFitPov,somFitNonPov,anxFitPov,anxFitNonPov,thoFitPov,thoFitNonPov,witFitPov,witFitNonPov,socFitPov,socFitNonPov,attFitPov,attFitNonPov,rulFitPov,rulFitNonPov,aggFitPov,aggFitNonPov)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/gpFitBoots_cbcl_pNp_sRE_poly_1k.rds')

print('done with g~p fit bootstrapping!')
