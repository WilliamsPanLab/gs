library(mgcv)
library(pammtools)
library(patchwork)
library(gratia)
library(scales)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','g','subjectkey','interview_age','sex','income')]
# get length of df for later
lenDF=dim(masterdf)[1]
# will need to get full and reduced models for each boot, as well as a null distribution
# in addition to derivatives and fits, save F values of interaction + null distribution F values
# interactions to be tested: p*sex, p*poverty, p*sex*poverty
# convert all cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
# sex (g in for x to get co-pilot around its NSFW filter) and poverty to factors
masterdf$seg<-as.ordered(masterdf$sex)
masterdf$poverty=0
masterdf$income<-as.numeric(masterdf$income)
# poverty now defined in sample construction
masterdf$poverty[masterdf$Pov_v2==1]=1
masterdf$poverty=as.factor(masterdf$poverty)
### initialize cross-boot vectors
# predicted derivatives: set to maximum value for ncol +1, as 0-maxvalue is 1 longer than maxvalue
pMaxVal=max(masterdf$cbcl_scr_syn_totprob_r)
F_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
M_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
P_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
R_pDeriv=matrix(0,nrow=10000,ncol=(pMaxVal+1))
# predicted values: set to maximum value for ncol
F_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
M_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
P_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
R_pFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
# a lil' vector just to track max p int and ext over iterations
pMax=rep(0,10000)
# and modular fit for pov vs. psuedopov 2 (matched #) fits
povFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
pseudopovFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
FullNonpovFit=matrix(0,nrow=10000,ncol=(pMaxVal+1))
#################################################
##################### loop over manual bootstraps
#################################################
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
	# get max for this iteration
	bpmax=max(bootSamp$cbcl_scr_syn_totprob_r)
	# last pre-step: need to create NULL variables. Use count of F and count of Poverty as count of membership to psuedo-groups, but randomly distribute membership
	# get count F in this boot
	Fcount=length(which(bootSamp$sex=="F"))
	# get count poverty in this boot
	Povcount=length(which(bootSamp$poverty==1))
	# randomly assign Fcount people to pseudoseg
	bootSamp$psuedoseg=0
	bootSamp$psuedoseg[sample(1:dim(bootSamp)[1],Fcount)]=1
	# randomly assign Povcount people to pseudopoverty
	bootSamp$psuedopoverty=0
	bootSamp$psuedopoverty[sample(1:dim(bootSamp)[1],Povcount)]=1
	# make a df that is just the same number of kids in poverty, but actually from the non-poverty group. The point is to see if we can recover g~p slope in a an equivalent number of non-poverty kids
	bootSamp$pseudopoverty2=0
	bootSamp$pseudopoverty2[sample(which(bootSamp$poverty==0),Povcount)]=1
	
	#
	######## I FIT MODELS
	#
	
	#### g as response variable: REDUCED MODELS
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age,k=4),data=bootSamp)
	#### FULL MODELS
	pgAge_seg<-bam(g~s(cbcl_scr_syn_totprob_r,by=seg,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+seg+s(interview_age,k=4),data=bootSamp)
	pgAge_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+poverty+s(interview_age,k=4),data=bootSamp)
	# comparison models to evaluate AIC gain from interactions specifically. Should include main effects but minus interaction of interest
	pgAge_seg_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r,k=4)+seg+s(interview_age,k=4),data=bootSamp)
	pgAge_pov_noIntrxn=bam(g~s(cbcl_scr_syn_totprob_r,k=4)+poverty+s(interview_age,k=4),data=bootSamp)
	pgAge_seg_pov_noIntrxn<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty,k=4)+s(cbcl_scr_syn_totprob_r,by=seg,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+seg*poverty+s(interview_age,k=4),data=bootSamp)
	#### VERY FULL MODELS
	pgAge_seg_pov<-bam(g~s(cbcl_scr_syn_totprob_r,by=poverty,k=4)+s(cbcl_scr_syn_totprob_r,by=seg,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+seg*poverty+s(interview_age,k=4)+s(cbcl_scr_syn_totprob_r, by = interaction(seg, poverty),k=4),data=bootSamp)
	# fit null models for use later (no reduced, reduced is the same as real reduced)
	#### FULL NULL MODELS
	pgAge_seg_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedoseg,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+psuedoseg+s(interview_age,k=4),data=bootSamp)
	pgAge_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+psuedopoverty+s(interview_age,k=4),data=bootSamp)
	#### VERY FULL NULL MODELS
	pgAge_seg_pov_n<-bam(g~s(cbcl_scr_syn_totprob_r,by=psuedopoverty,k=4)+s(cbcl_scr_syn_totprob_r,by=psuedoseg,k=4)+s(cbcl_scr_syn_totprob_r,k=4)+psuedoseg*psuedopoverty+s(interview_age,k=4)+s(cbcl_scr_syn_totprob_r, by = interaction(psuedoseg, psuedopoverty),k=4),data=bootSamp)
	
	#
	######## II PREDICT VARIABLE OF INTEREST WITH FIT SPLINE
	#
	
	# use PREDICTED VALUES of model fit for each symptom count for saving
	eachPcount=seq(0:bpmax)
	# set age to to median for predict df, +1 because 0:max sequence is one lnger than max
	predictDFp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)))
	# set predict df to interacting factors of interest (seg, poverty)
	predictDF_segp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("F",(bpmax+1)))
	predictDF_povp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("1",(bpmax+1)))
	predictDF_segpovp=data.frame(eachPcount,rep(median(bootSamp$interview_age),(bpmax+1)),rep("F",(bpmax+1)),rep("1",(bpmax+1)))
	# set colnames so predict can work
	colnames(predictDF_segp)=c('cbcl_scr_syn_totprob_r','interview_age','seg')
	colnames(predictDF_povp)=c('cbcl_scr_syn_totprob_r','interview_age','poverty')
	colnames(predictDF_segpovp)=c('cbcl_scr_syn_totprob_r','interview_age','seg','poverty')
	# predict girl
	forfitp_F=predict(pgAge_seg,predictDF_segp)
	forDerivp_F=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=(bpmax+1))
	# predict boy
	predictDF_segp$seg=rep("M",(bpmax+1))
	forfitp_M=predict(pgAge_seg,predictDF_segp)
	forDerivp_M=derivatives(pgAge_seg,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_segp,n=(bpmax+1))
	# predict poverty
	forfitp_P=predict(pgAge_pov,predictDF_povp)
	forDerivp_P=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=(bpmax+1))
	# predict rich
	predictDF_povp$poverty=rep("0",(bpmax+1))
	forfitp_R=predict(pgAge_pov,predictDF_povp)
	forDerivp_R=derivatives(pgAge_pov,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_povp,n=(bpmax+1))
	# print out fit
	F_pFit[b,1:(bpmax+1)]=forfitp_F
	M_pFit[b,1:(bpmax+1)]=forfitp_M
	P_pFit[b,1:(bpmax+1)]=forfitp_P
	R_pFit[b,1:(bpmax+1)]=forfitp_R
	# print out derivatives
	F_pDeriv[b,1:(bpmax+1)]=forDerivp_F$derivative
	M_pDeriv[b,1:(bpmax+1)]=forDerivp_M$derivative
	P_pDeriv[b,1:(bpmax+1)]=forDerivp_P$derivative
	R_pDeriv[b,1:(bpmax+1)]=forDerivp_R$derivative
	# print out max of unconverted versions to anchor em later
	pMax[b]=bpmax
	#
	# III modular portion for poverty and equivalently-sized group from non-poverty:
	# extract poverty and pseudopov2 group
	povertyGroup=subset(bootSamp,poverty==1)
	pseudoPov2Group=subset(bootSamp,pseudopoverty2==1)
	nonpovGroup=subset(bootSamp,poverty==0)
	# find maximum cbcl_scr_syn_totprob_r value present across both subgroups, make prediction df with that range of cbcl_scr_syn_totprob_r values
	bpmax_2=min(c(max(povertyGroup$cbcl_scr_syn_totprob_r),max(pseudoPov2Group$cbcl_scr_syn_totprob_r)))
	eachPcount2=seq(0,bpmax_2)
	predictDF_2=data.frame(eachPcount2,rep(median(bootSamp$interview_age),(bpmax_2+1)))
	# fit g~p to poverty group
	pov_gp=gam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age),data=povertyGroup)
	# fit g~p to same-size nonpoverty group
	nonpov_gp=gam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age),data=pseudoPov2Group)
	# fit g~p to the full nonpoverty group
	fnonpov_gp=gam(g~s(cbcl_scr_syn_totprob_r,k=4)+s(interview_age),data=nonpovGroup)
	# make prediction df to use for both models
	colnames(predictDF_2)=c('cbcl_scr_syn_totprob_r','interview_age')
	# extract fit derivatives
	povFit[b,1:(bpmax_2+1)]=derivatives(pov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative
	pseudopovFit[b,1:(bpmax_2+1)]=derivatives(nonpov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative
	FullNonpovFit[b,1:(bpmax_2+1)]=derivatives(fnonpov_gp,term="s(cbcl_scr_syn_totprob_r)",data=predictDF_2,n=(bpmax_2+1))$derivative

}
# SAVEOUT
# save out version with all F stats and AICs, include max values for all iterations as well
# save out fits
outdf=data.frame(F_pFit,M_pFit,P_pFit,R_pFit)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpFits.rds')
# save out derivatives
outdf=data.frame(F_pDeriv,M_pDeriv,P_pDeriv,R_pDeriv)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpDerivs.rds')
print('done with g~p fit bootstrapping!')
# save out modular poverty and equivlanetly-sized poverty fits
outdf=data.frame(povFit,pseudopovFit,FullNonpovFit)
saveRDS(o`utdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/F3_gpPovNonPov.rds')
