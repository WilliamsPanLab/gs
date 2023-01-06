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

#### gpAge_te
if (!file.exists("/scratch/users/apines/gp/gpAge_te.rds")){
	gpAge_te<-bam(cbcl_scr_syn_totprob_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	# resilient version
	rdata_file = file("/scratch/users/apines/gp/gpAge_te.rds", blocking = TRUE)
	saveRDS(gpAge_te, file=rdata_file)
	close(rdata_file)
} else {
	print('found gpAge_te.rds. Loading')
	gpAge_te=readRDS("/scratch/users/apines/gp/gpAge_te.rds")
}
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_te.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(gpAge_te)
dev.off()
# print confidence intervals for supplemental figure
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_te_ci.png',width=2400,height=800)
# prinout gg_tensor
gg_tensor(gpAge_te,ci=T)
dev.off()

#### gpAge_independent splines
if (!file.exists("/scratch/users/apines/gp/gpAge.rds")){
	print(dim(masterdf))
	gpAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/gpAge.rds", blocking = TRUE)
	saveRDS(gpAge, file=rdata_file)
	close(rdata_file)
} else {
	print('found gpAge.rds. Loading')
	gpAge=readRDS("/scratch/users/apines/gp/gpAge.rds")
}
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_derivs_gp.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(gpAge,'g')
dev.off()

##############################
######## I SCATTERPLOTS ON TWO VARIABLES OF INTEREST WITH THEIR FIT SPLINE
#### g as response variable
if (!file.exists("/scratch/users/apines/gp/pgAge.rds")){
	print('fitting p~g on masterdf')
	print('master df dims')
	print(dim(masterdf))
	pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf)
	rdata_file = file("/scratch/users/apines/gp/pgAge.rds", blocking = TRUE)
	saveRDS(pgAge, file=rdata_file)
	close(rdata_file)
} else {
	print(' found pgAge.rds. Loading.')
	pgAge=readRDS('/scratch/users/apines/gp/pgAge.rds')
}
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_derivs_pg.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(pgAge,'cbcl_scr_syn_totprob_r')
dev.off()
# use model fit for geom_smooth plotting
forSpline<-predict(pgAge, newdata = masterdf)
print('length of fit:')
print(length(forSpline))
plotdf<-data.frame(masterdf$cbcl_scr_syn_totprob_r,forSpline,masterdf$g,as.factor(masterdf$eventname))
colnames(plotdf)<-c('totprob','predicted','g','eventname')
# make scatter plot showing g~total problems with their spline fit
thePlot<-ggplot(plotdf,aes(totprob,g))+
geom_point(alpha=.1,aes(color=eventname))+
geom_smooth(aes(y=predicted),size=2,se=F,color='black')+
theme_classic(base_size=24)+
xlab('Total Problems Score')+
scale_color_discrete(name="Vist",breaks=c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),labels=c("Baseline","2 Years"))+
guides(colour = guide_legend(override.aes = list(alpha = 1)))
# make the image to be printed
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_Scatter_pg.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()
######## II FORMALLY TEST FOR NON-LINEARITY
#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age)+s(subjectkey,bs='re')+ti(cbcl_scr_syn_totprob_r,interview_age),data=masterdf)
summary(pgAge)



############################## USE THIS CHUNK ON ti() MODEL IF ti() FITS: OTHERWISE, KEEP ON NO-INTERACTIONS MODEL
### step one: extract max fitted values for each tensor plot to prospectively set color range of plots
print('Max and Min of gpAge_te (no cbcl61) fitted values')
print(max(gpAge_te$fitted.values))
print(min(gpAge_te$fitted.values))
#gg_tensor(fixedEfModel)+scale_fill_gradient2(
#    low = muted("blue"), 
#    mid = "white", 
#    high = muted("red"), 
#    midpoint = 0, limits=c(-2,2)
#)
###
### AND USE max(fixedEfModel$fitted.values) OF ALL CBCL61 TE FITS AND ALL NO CBCL61 TE FITS FOR SAME COLOR SCALE, AT LEAST CONISTENTLY ACROSS CBCL61-NOCBCL61 WITHIN SYMPTOM TYPE CATEGORIES
### 
### PLOT TENSORS WITH
### low = muted("blue"), 
#     mid = "white", 
#     high = muted("red"), 
#     midpoint = 0, limits=c(-2,2)
#)
####

print('testing for interaction between g and age on totprobs')
### interaction test
gpAge_intrxn<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+ti(g,interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
print('model with interaction (ti) for age and g')
print(summary(gpAge_intrxn))


print('full vs reduced model for effect size and anova of g~totprob')

###### FULL VS REDUCED ANOVA: P AND DR2
if (!file.exists("/scratch/users/apines/gp/noG_Gam_p.rds")){
        no_g_Gam<-bam(cbcl_scr_syn_totprob_r~s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
	rdata_file = file("/scratch/users/apines/gp/noG_Gam_p.rds", blocking = TRUE)
        saveRDS(no_g_Gam, file=rdata_file)
        close(rdata_file)
} else {
        print('found /scratch/users/apines/gp/noG_Gam_p.rds. Loading')
        no_g_Gam=readRDS("/scratch/users/apines/gp/noG_Gam_p.rds")
}
no_g_Sum<-summary(no_g_Gam)

# g-included model for measuring difference
gGam<-gpAge
gSum<-summary(gGam)
dif<-gSum$r.sq-no_g_Sum$r.sq
print('difference in r^2: g vs. no g in totprobs model')
print(dif)
# test of dif with anova.gam
anovaRes<-anova.gam(no_g_Gam,gGam,test='Chisq')
anovaP<-anovaRes$`Pr(>Chi)`
anovaP2<-unlist(anovaP)
print('chi-sq p value: g vs. no g in totprobs model')
print(anovaP2[2])

print('fitting te() with cbcl61 included + visualizing')
# fit version with cbcl 61
if (!file.exists("/scratch/users/apines/gp/gpAge_te_cbcl61.rds")){
	mixedEfModel<-bam(cbcl_scr_syn_totprob_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf,family=nb())
	# save
	rdata_file = file("/scratch/users/apines/gp/gpAge_te_cbcl61.rds", blocking = TRUE)
        saveRDS(mixedEfModel, file=rdata_file)
        close(rdata_file)
} else {
	print('found gpAge_te_cbcl61.rds. Loading')
	mixedEfModel=readRDS('/scratch/users/apines/gp/gpAge_te_cbcl61.rds')
}
######## plot tensor 
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_te_cbcl61.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

print('Max and Min of gpAge_te (with cbcl61) fitted values')
print(max(mixedEfModel$fitted.values))
print(min(mixedEfModel$fitted.values))

print('done with g~p fits')
