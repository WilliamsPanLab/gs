library(mgcv)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
# convert family id to factor
masterdf$rel_family_id=as.factor(masterdf$rel_family_id)

#### gpAge_independent splines
nb_gpAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf,family=nb())
####### plot residuals
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_nb_residuals_qq.png',width=700,height=700)
# prinout gg_derivs
qq.gam(nb_gpAge)
dev.off()

# gaussian
#### gpAge_independent splines
gaussian_gpAge<-bam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re')+s(rel_family_id,bs='re'),data=masterdf)
####### plot residuals
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_gaus_residuals_qq.png',width=700,height=700)
# prinout gg_derivs
qq.gam(gaussian_gpAge)
dev.off()

# print AIC and BIC
AIC(nb_gpAge,gaussian_gpAge)
BIC(nb_gpAge,gaussian_gpAge)
