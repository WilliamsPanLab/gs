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
# convert family id to factor
masterdf$rel_family_id=as.factor(masterdf$rel_family_id)
# FYI
print('dimensions of dataframe')
dim(masterdf)

#### PRINT TENSORS
print('data loaded, fitting te(age,g) for internalizing symptoms, REs for participant and family')
### internalizing
if (!file.exists("/scratch/users/apines/gp/g_int_Age_te.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_internal_r~te(interview_age,g)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
	# save
	rdata_file = file("/scratch/users/apines/gp/g_int_Age_te.rds", blocking = TRUE)
	saveRDS(mixedEfModel, file=rdata_file)
	close(rdata_file)
} else {
	print('found g_int_Age_te.rds. Loading')
	mixedEfModel=readRDS("/scratch/users/apines/gp/g_int_Age_te.rds")
}
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_te.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()
print('tensor printed. Printing confidence intervals')
# print confidence intervals for supplemental figure
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_te_ci.png',width=2400,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel,ci=T)
dev.off()

print('fitting te(age,g) for internalizing symptoms with cbcl61 as covariate, REs for participant and family')
# FIT WITH CBCL61
if (!file.exists("/scratch/users/apines/gp/g_int_Age_te_cbcl61.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_internal_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
	# save
	rdata_file = file("/scratch/users/apines/gp/g_int_Age_te_cbcl61.rds", blocking = TRUE)
	saveRDS(mixedEfModel, file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/g_int_Age_te_cbcl61.rds. Loading')
	mixedEfModel=readRDS("/scratch/users/apines/gp/g_int_Age_te_cbcl61.rds")
}
# print png of tensor
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_te_cbcl61.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

print('calculating effect size for g within internalizing symptoms')

# fit as independent splines to get p and dr2 for int w/ and without g
###### FULL VS REDUCED ANOVA: P AND DR2

if (!file.exists("/scratch/users/apines/gp/noG_Gam_int.rds")){
	no_g_Gam<-bam(cbcl_scr_syn_internal_r~s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/noG_Gam_int.rds", blocking = TRUE)
        saveRDS(no_g_Gam, file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/noG_Gam_int.rds. Loading')
	no_g_Gam=readRDS("/scratch/users/apines/gp/noG_Gam_int.rds")
}
no_g_Sum<-summary(no_g_Gam)
# g-included model for measuring difference: note gGam gets used later for comparison with ti tensor
if (!file.exists("/scratch/users/apines/gp/G_Gam_int.rds")){
	gGam<-bam(cbcl_scr_syn_internal_r~s(g)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/G_Gam_int.rds", blocking = TRUE)
        saveRDS(gGam, file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/G_Gam_int.rds. Loading')
	gGam=readRDS("/scratch/users/apines/gp/G_Gam_int.rds")
}

gSum<-summary(gGam)
dif<-gSum$r.sq-no_g_Sum$r.sq
print('difference in r^2: g vs. no g in internalizing model')
print(dif)
# test of dif with anova.gam
anovaRes<-anova.gam(no_g_Gam,gGam,test='Chisq')
anovaP<-anovaRes$`Pr(>Chi)`
anovaP2<-unlist(anovaP)
print('chi-sq p value: g vs. no g in internalizing model')
print(anovaP2[2])

print('testing for interaction of age*g on internalizing symptoms')
# FIT AS INDEPENDENT SPLINES TO TEST FOR INTERACTION: try anova.gam with and without ti interaction (NOTE NO CBCL61)
if (!file.exists("/scratch/users/apines/gp/g_int_Age_ti.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_internal_r~s(interview_age)+s(g)+ti(interview_age,g)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
	rdata_file=file("/scratch/users/apines/gp/g_int_Age_ti.rds", blocking = TRUE)
	saveRDS(mixedEfModel,file=rdata_file)
	close(rdata_file)
} else {
	print('found /scratch/users/apines/gp/g_int_Age_ti.rds. Loading')
	mixedEfModel=readRDS("/scratch/users/apines/gp/g_int_Age_ti.rds")
}
# print png of tensor
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_ti_nocbcl61.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()
# ti-included model for measuring difference
tiSum<-summary(mixedEfModel)
print(tiSum)
notiSum<-summary(gGam)
dif<-tiSum$r.sq-gGam$r.sq
print('difference in r^2: ti vs. no ti(age,g) in internalizing model')
print(dif)
# test of dif with anova.gam
anovaRes<-anova.gam(gGam,mixedEfModel,test='Chisq')
anovaP<-anovaRes$`Pr(>Chi)`
anovaP2<-unlist(anovaP)
print('chi-sq p value: ti vs. no ti(age,g) in internalizing model')
print(anovaP2[2])

######## I SCATTERPLOTS ON TWO VARIABLES OF INTEREST WITH THEIR FIT SPLINE
#### g as response variable
if (!file.exists("/scratch/users/apines/gp/IntgAge.rds")){
	pgAge<-bam(g~s(cbcl_scr_syn_internal_r)+s(interview_age)+ti(cbcl_scr_syn_internal_r,interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf)
        rdata_file = file("/scratch/users/apines/gp/IntgAge.rds", blocking = TRUE)
	saveRDS(pgAge, file=rdata_file)
	close(rdata_file)
} else {
	print(' found IntgAge.rds. Loading.')
        pgAge=readRDS('/scratch/users/apines/gp/IntgAge.rds')
}
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/IntgAge_derivs_pg.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(pgAge,'cbcl_scr_syn_internal_r')
dev.off()
# use model fit for geom_smooth plotting
forSpline<-predict(pgAge, newdata = masterdf)
plotdf<-data.frame(masterdf$cbcl_scr_syn_internal_r,forSpline,masterdf$g,as.factor(masterdf$eventname))
colnames(plotdf)<-c('internal','predicted','g','eventname')
# make scatter plot showing g~total problems with their spline fit
thePlot<-ggplot(plotdf,aes(internal,g))+
geom_point(alpha=.1,aes(color=eventname))+
geom_smooth(aes(y=predicted),size=2,se=F,color='black')+
theme_classic(base_size=24)+
xlab('Internalizing Problems Score')+
scale_color_discrete(name="Vist",breaks=c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),labels=c("Baseline","2 Years"))+
guides(colour = guide_legend(override.aes = list(alpha = 1)))
# make the image to be printed
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/IntgAge_Scatter_pg.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()

############ II FORMALLY TEST FOR NON-LINEARITY
#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
IntgAge<-bam(g~s(cbcl_scr_syn_internal_r,m=c(2,0))+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re')+ti(cbcl_scr_syn_internal_r,interview_age),data=masterdf)
summary(IntgAge)

print('done with internalizing')
