library(mgcv)
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')
gpAge<-gam(cbcl_scr_syn_totprob_r~s(g)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf,family=nb())
gpAge_te<-gam(cbcl_scr_syn_totprob_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# resilient version
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/gpAge_te.rds", blocking = TRUE)
saveRDS(gpAge_te, file=rdata_file)
close(rdata_file)
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/gpAge.rds", blocking = TRUE)
saveRDS(gpAge, file=rdata_file)
close(rdata_file)
################
# Comparative p, int, and ext with and without controlling for cbcl_q61_p
#########
### P
#mixedEfModel<-gam(cbcl_scr_syn_totprob_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
#rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/gpAge_te.rds", blocking = TRUE)
#saveRDS(gpAge_te, file=rdata_file)

mixedEfModel<-gam(cbcl_scr_syn_totprob_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/gpAge_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)

### internalizing
mixedEfModel<-gam(cbcl_scr_syn_internal_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_int_Age_te.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)

mixedEfModel<-gam(cbcl_scr_syn_internal_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_int_Age_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)

### externalizing
mixedEfModel<-gam(cbcl_scr_syn_external_r~te(interview_age,g)+s(subjectkey,bs='re'),data=masterdf,family=nb())
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_ext_Age_te.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)

mixedEfModel<-gam(cbcl_scr_syn_external_r~te(interview_age,g)+cbcl_q61_p+s(subjectkey,bs='re'),data=masterdf)
# save
rdata_file = file("/oak/stanford/groups/leanew1/users/apines/data/gp/g_ext_Age_te_cbcl61.rds", blocking = TRUE)
saveRDS(mixedEfModel, file=rdata_file)
close(rdata_file)
