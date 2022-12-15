library(mgcv)
library(pammtools)
library(visreg)

masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
# convert family id to factor
masterdf$rel_family_id=as.factor(masterdf$rel_family_id)

### externalizing

#### LOAD IN THE MODEL FROM RDS INSTEAD

mixedEfModel<-bam(cbcl_scr_syn_external_r~te(interview_age,g)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/scratch/users/apines/gp/g_ext_Age_te.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExtAge_te.png',width=1200,height=1200)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

# with cbcl61
mixedEfModel<-bam(cbcl_scr_syn_external_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/scratch/users/apines/gp/g_ext_Age_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExtAge_te_cbcl61.png',width=1200,height=1200)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

### FULL VS REDUCED P + DR2
# fit as independent splines to get p and dr2 for Ext w/ and w/o g
no_g_Gam<-bam(cbcl_scr_syn_external_r~s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
no_g_Sum<-summary(no_g_Gam)
# g-included model for measuring difference
gGam<-bam(cbcl_scr_syn_external_r~s(g)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
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

### SCATTER OF EXTERNAL ALONG G after controlling for age
# use no g gam to control for age in scatterplot of g~ext (scatter on residuals)
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gExt_ageRegressed.png',width=700,height=700)
# prinout gg_tensor
visreg(gGam,'g')
dev.off()
