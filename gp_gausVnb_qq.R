library(mgcv)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')

#### gpAge_independent splines
gpAge<-gam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
####### plot residuals
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_nb_residuals_qq.png',width=1500,height=1500)
# prinout gg_derivs
qq.gam(gpAge)
dev.off()

# gaussian
#### gpAge_independent splines
gpAge<-gam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf)
####### plot residuals
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_gaus_residuals_qq.png',width=1500,height=1500)
# prinout gg_derivs
qq.gam(gpAge)
dev.off()
