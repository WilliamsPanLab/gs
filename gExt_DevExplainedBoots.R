library(mgcv)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
devExplained_full=rep(0,1000)
devExplained_NoCbcl61=rep(0,1000)
devExplained_NoTot=rep(0,1000)
# garner num subjs for bootstrapping
numSubjs=dim(masterdf)[1]
# loop over manual bootstrap
for (b in 1:10000){
	print(b)
	BootInd<-sample(seq(1,numSubjs),numSubjs,replace=T)
	# bootstrap sample
	bootSamp=masterdf[BootInd,]
	#### gpAge_independent splines
	gpAge_full<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age)+cbcl_q61_p,data=bootSamp)
	devExplained_full[b]<-summary(gpAge_full)$dev.expl
	# fit version with cbcl 61
	gpAge_nocbcl61<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age),data=bootSamp)
	devExplained_NoCbcl61[b]<-summary(gpAge_nocbcl61)$dev.expl
	# fit version with just cbcl61
	gpAge_noTotProbs<-bam(g~cbcl_q61_p+s(interview_age),data=bootSamp)
	devExplained_NoTot[b]<-summary(gpAge_noTotProbs)$dev.expl
}

# saveout df of dev explained for plotting
outdf=data.frame(devExplained_full,devExplained_NoCbcl61,devExplained_NoTot)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/DevExplainedExt.rds')

