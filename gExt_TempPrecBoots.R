library(mgcv)
# load in cross-TP data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/OutDFTmpPrec.rds')
dr2adj_G=rep(0,10000)
dr2adj_PP=rep(0,10000)
dr2adj_Gr=rep(0,10000)
# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_external_r','g','subjectkey','interview_age','Grades')]
# loop over manual bootstrap
for (b in 1:10000){
	print(b)
	# get subjects to include in this bootstrap
	BootSubjs=sample(seq(1,numSubjs),numSubjs,replace=T)	
	# bootstrap sample
	bootSamp=masterdf[BootSubjs,]
	# dr2 for G
	pgAge_full<-summary(bam(cbcl_scr_syn_external_r.y~+s(cbcl_scr_syn_external_r.x)+s(parentPcount.x)+s(g.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	# reduced
	pgAge_reduced<-summary(bam(cbcl_scr_syn_external_r.y~s(cbcl_scr_syn_external_r.x)+s(parentPcount.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	dr2adj_G[b]=pgAge_full$r.sq-pgAge_reduced$r.sq
	# dr2 for parentP
	pgAge_full<-summary(bam(cbcl_scr_syn_external_r.y~s(cbcl_scr_syn_external_r.x)+s(parentPcount.x)+s(g.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	# reduced
	pgAge_reduced<-summary(bam(cbcl_scr_syn_external_r.y~s(cbcl_scr_syn_external_r.x)+s(g.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	dr2adj_PP[b]=pgAge_full$r.sq-pgAge_reduced$r.sq
	# dr2 for Grades
	pgAge_full<-summary(bam(cbcl_scr_syn_external_r.y~s(cbcl_scr_syn_external_r.x)+s(parentPcount.x)+s(g.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	# reduced
	pgAge_reduced<-summary(bam(cbcl_scr_syn_external_r.y~s(cbcl_scr_syn_external_r.x)+s(parentPcount.x)+s(g.x)+s(interview_age.x)+Grades.x,data=bootSamp,family=nb()))
	dr2adj_Gr[b]=pgAge_full$r.sq-pgAge_reduced$r.sq
}
# saveout df of dev explained for plotting
outdf=data.frame(dr2adj_G,dr2adj_PP,dr2adj_Gr)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/dr2_Ext.rds')

