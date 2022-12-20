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

#### gpAge_te
gpAge_te<-bam(cbcl_scr_syn_totprob_r~te(interview_age,g)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
# resilient version
rdata_file = file("/scratch/users/apines/gp/gpAge_te.rds", blocking = TRUE)
saveRDS(gpAge_te, file=rdata_file)
close(rdata_file)
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
gpAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
rdata_file = file("/scratch/users/apines/gp/gpAge.rds", blocking = TRUE)
saveRDS(gpAge, file=rdata_file)
close(rdata_file)
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_derivs_gp.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(gpAge,'g')
dev.off()

############################## USE THIS CHUNK ON ti() MODEL IF ti() FITS: OTHERWISE, KEEP ON NO-INTERACTIONS MODEL
############################## USE THIS CHUNK ON ti() MODEL IF ti() FITS: OTHERWISE, KEEP ON NO-INTERACTIONS MODEL
######## I SCATTERPLOTS ON TWO VARIABLES OF INTEREST WITH THEIR FIT SPLINE
#### g as response variable
pgAge<-bam(g~s(cbcl_scr_syn_totprob_r)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf)
rdata_file = file("/scratch/users/apines/gp/pgAge.rds", blocking = TRUE)
saveRDS(gpAge, file=rdata_file)
close(rdata_file)
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_derivs_pg.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(pgAge,'cbcl_scr_syn_totprob_r')
dev.off()
# use model fit for geom_smooth plotting
forSpline<-predict(pgAge, data = masterdf)
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
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_Scater_pg.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()
######## II FORMALLY TEST FOR NON-LINEARITY - SWITCH TO TI() IF THAT IS CHOSEN MODEL, KEEP HERE IF NOT
#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
pgAge<-bam(g~s(cbcl_scr_syn_totprob_r,m=c(2,0))+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf)
summary(pgAge)
############################## USE THIS CHUNK ON ti() MODEL IF ti() FITS: OTHERWISE, KEEP ON NO-INTERACTIONS MODEL
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
############################## USE THIS CHUNK ON ti() MODEL IF ti() FITS: OTHERWISE, KEEP ON NO-INTERACTIONS MODEL

print('testing for interaction between g and age on totprobs')
### interaction test
gpAge_intrxn<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+ti(g,interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
print('model with interaction (ti) for age and g')
print(summary(gpAge_intrxn))


print('full vs reduced model for effect size and anova of g~totprob')
###### FULL VS REDUCED ANOVA: P AND DR2
no_g_Gam<-bam(cbcl_scr_syn_totprob_r~s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
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

print('generating scatter plot of externalizing ~ g (residualized)')
### SCATTER OF EXTERNAL ALONG G after controlling for age
# use no g gam to control for age in scatterplot of g~ext (scatter on residuals)
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gp_ageRegressed.png',width=700,height=700)
# prinout gg_tensor
visreg(gGam,'g')
dev.off()

print('fitting te() with cbcl61 included + visualizing')
# fit version with cbcl 61
mixedEfModel<-bam(cbcl_scr_syn_totprob_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
######## plot tensor 
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_te_cbcl61.png',width=800,height=800)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()
# save
rdata_file = file("/scratch/users/apines/gp/gpAge_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)

