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

masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

print('fitting externalizing ~ te(age,g,) with REs for subj and family')
### externalizing
if (!file.exists("/scratch/users/apines/gp/g_ext_Age_te.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_external_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	# save
	rdata_file = file("/scratch/users/apines/gp/g_ext_Age_te.rds", blocking = TRUE)
	saveRDS(mixedEfModel, file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/g_ext_Age_te.rds. Loading')
	mixedEfModel=readRDS("/scratch/users/apines/gp/g_ext_Age_te.rds")
}
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExtAge_te.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()
# print confidence intervals for supplemental figure
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExtAge_te_ci.png',width=2400,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel,ci=T)
dev.off()

print('fitting externalizing ~ te(age,g,)+cbcl61 with REs for subj and family')
# with cbcl61
if (!file.exists("/scratch/users/apines/gp/g_ext_Age_te_cbcl61.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_external_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf,family=nb())
	# save
	rdata_file = file("/scratch/users/apines/gp/g_ext_Age_te_cbcl61.rds", blocking = TRUE)
	saveRDS(mixedEfModel, file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/g_ext_Age_te_cbcl61.rds. Loading')
	mixedEfModel=readRDS("/scratch/users/apines/gp/g_ext_Age_te_cbcl61.rds")
}
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExtAge_te_cbcl61.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

print('testing externalizing~g*age interaction with REs for subj and family')
# interaction
if (!file.exists("/scratch/users/apines/gp/g_ext_Age_ti.rds")){
	gGam_intrxn<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age)+ti(g,interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	print('externalizing ~ ti(g,age) included model')
	rdata_file = file("/scratch/users/apines/gp/g_ext_Age_ti.rds", blocking = TRUE)
        saveRDS(gGam_intrxn, file=rdata_file)
        close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/g_ext_Age_ti.rds. Loading')
	gGam_intrxn=readRDS("/scratch/users/apines/gp/g_ext_Age_ti.rds")
}
summary(gGam_intrxn)

### FULL VS REDUCED P + DR2
print('calculating effect size for g within externalizing symptoms')
# fit as independent splines to get p and dr2 for Ext w/ and w/o g
if (!file.exists("/scratch/users/apines/gp/noG_Gam_ext.rds")){
	no_g_Gam<-bam(cbcl_scr_syn_external_r~s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/noG_Gam_ext.rds", blocking = TRUE)
        saveRDS(no_g_Gam, file=rdata_file)
        close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/noG_Gam_ext.rds. Loading')
	no_g_Gam=readRDS("/scratch/users/apines/gp/noG_Gam_ext.rds")
}

no_g_Sum<-summary(no_g_Gam)
# g-included model for measuring difference
if (!file.exists("/scratch/users/apines/gp/G_Gam_ext.rds")){
	gGam<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/G_Gam_ext.rds", blocking = TRUE)
        saveRDS(gGam, file=rdata_file)
        close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/G_Gam_ext.rds. Loading')
	gGam=readRDS("/scratch/users/apines/gp/G_Gam_ext.rds")
}

gSum<-summary(gGam)
dif<-gSum$r.sq-no_g_Sum$r.sq
print('difference in r^2: g vs. no g in externalizing model')
print(dif)
# test of dif with anova.gam
anovaRes<-anova.gam(no_g_Gam,gGam,test='Chisq')
anovaP<-anovaRes$`Pr(>Chi)`
anovaP2<-unlist(anovaP)
print('chi-sq p value: g vs. no g in externalizing model')
print(anovaP2[2])

######## I SCATTERPLOTS ON TWO VARIABLES OF INTEREST WITH THEIR FIT SPLINE
#### g as response variable
if (!file.exists("/scratch/users/apines/gp/ExtgAge.rds")){
	pgAge<-bam(g~s(cbcl_scr_syn_external_r)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf)
        rdata_file = file("/scratch/users/apines/gp/ExtgAge.rds", blocking = TRUE)
	saveRDS(pgAge, file=rdata_file)
	close(rdata_file)
} else {
	print(' found ExtgAge.rds. Loading.')
        pgAge=readRDS('/scratch/users/apines/gp/ExtgAge.rds')
}
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/ExtgAge_derivs_pg.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(pgAge,'cbcl_scr_syn_external_r')
dev.off()
# use model fit for geom_smooth plotting
forSpline<-predict(pgAge, newdata = masterdf)
plotdf<-data.frame(masterdf$cbcl_scr_syn_external_r,forSpline,masterdf$g,as.factor(masterdf$eventname))
colnames(plotdf)<-c('external','predicted','g','eventname')
# make scatter plot showing g~total problems with their spline fit
thePlot<-ggplot(plotdf,aes(external,g))+
geom_point(alpha=.1,aes(color=eventname))+
geom_smooth(aes(y=predicted),size=2,se=F,color='black')+
theme_classic(base_size=24)+
xlab('Externalizing Problems Score')+
scale_color_discrete(name="Vist",breaks=c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),labels=c("Baseline","2 Years"))+
guides(colour = guide_legend(override.aes = list(alpha = 1)))
# make the image to be printed
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/ExtgAge_Scatter_pg.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()

############ II FORMALLY TEST FOR NON-LINEARITY
#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
ExtgAge<-bam(g~s(cbcl_scr_syn_external_r,m=c(2,0))+s(interview_age)+s(subjectkey,bs='re'),data=masterdf)
summary(ExtgAge)

print('done with externalizing')
