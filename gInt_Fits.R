library(mgcv)
library(pammtools)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')

### internalizing
mixedEfModel<-gam(cbcl_scr_syn_internal_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_int_Age_te.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_te.png',width=1500,height=1500)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

# FIT WITH CBCL61
mixedEfModel<-gam(cbcl_scr_syn_internal_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_int_Age_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)
# print png of tensor
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_te_cbcl61.png',width=1500,height=1500)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off()

# fit as independent splines to get p and dr2 for int w/ and without g
###### FULL VS REDUCED ANOVA: P AND DR2
no_g_Gam<-gam(cbcl_scr_syn_internal_r~s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
no_g_Sum<-summary(no_g_Gam)
# g-included model for measuring difference
gGam<-gam(cbcl_scr_syn_internal_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
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


# FIT AS INDEPENDENT SPLINES TO TEST FOR INTERACTION: try anova.gam with and without ti interaction (NOTE NO CBCL61)
mixedEfModel<-gam(cbcl_scr_syn_internal_r~s(interview_age)+s(g)+ti(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# print png of tensor
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gIntAge_ti_nocbcl61.png',width=1500,height=1500)
# prinout gg_tensor
gg_tensor(mixedEfModel)
dev.off() 
reduced_mixedEfModel<-gam(cbcl_scr_syn_internal_r~s(interview_age)+s(g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# ti-included model for measuring difference
tiSum<-summary(mixedEfModel)
notiSum<-summary(reduced_mixedEfModel)
dif<-tiSum$r.sq-reduced_mixedEfModel$r.sq
print('difference in r^2: ti vs. no ti(age,g) in internalizing model')
print(dif)
# test of dif with anova.gam
anovaRes<-anova.gam(reduced_mixedEfModel,mixedEfModel,test='Chisq')
anovaP<-anovaRes$`Pr(>Chi)`
anovaP2<-unlist(anovaP)
print('chi-sq p value: ti vs. no ti(age,g) in internalizing model')
print(anovaP2[2])

