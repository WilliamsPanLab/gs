library(mgcv)
library(ggplot2)
library(tidyr)
library(dplyr)

# load in dataframe
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')

# ensure family ID as factor
masterdf$rel_family_id=as.factor(masterdf$rel_family_id)

# version where only those with familial representation are in df
famDF=masterdf[!is.na(masterdf$rel_family_id),]

# load in model
Model=readRDS('/scratch/users/apines/gp/pgAge.rds')

#### WITH G AS RESPONSE

# make "predict" version for plotting RE's
newData <- tidyr::expand(famDF, nesting(subjectkey,rel_family_id), cbcl_scr_syn_totprob_r=unique(cbcl_scr_syn_totprob_r))

# fit at mean age
newData$interview_age=11

# predicted p
predP <- bind_cols(newData, as.data.frame(predict(Model,newdata=newData,se.fit=TRUE)))

# cut out real values of famDF
famDFvals<-famDF[,c("subjectkey","g","cbcl_scr_syn_totprob_r","rel_family_id")]

# Plot it out
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_REs.png',width=11000,height=11000)
# prinout random effects
ggplot(predP, aes(x = cbcl_scr_syn_totprob_r, y = fit, group = subjectkey, color=rel_family_id))+geom_line()+geom_point(data=famDFvals,aes(y=g))+facet_wrap(~ subjectkey)
dev.off()

################

### family-level random effects
# select just a few families for viz
index <- (famDF$rel_family_id == 1) | (famDF$rel_family_id == 2) | (famDF$rel_family_id == 3) | (famDF$rel_family_id == 4) | (famDF$rel_family_id == 5) | (famDF$rel_family_id == 6) | (famDF$rel_family_id == 7) | (famDF$rel_family_id == 8) | (famDF$rel_family_id == 9) | (famDF$rel_family_id == 10) | (famDF$rel_family_id == 11) | (famDF$rel_family_id == 12)
famDFvals_select=famDF[index,]

newData <- tidyr::expand(famDFvals_select, nesting(rel_family_id,subjectkey), cbcl_scr_syn_totprob_r=unique(cbcl_scr_syn_totprob_r))

# fit at mean age
newData$interview_age=11

# predicted p
predP <- bind_cols(newData, as.data.frame(predict(Model,newdata=newData,se.fit=TRUE)))

# cut out real values of famDF
famDFvals_select<-famDFvals_select[,c("subjectkey","g","cbcl_scr_syn_totprob_r","rel_family_id")]

### PLOT FAMILY INSTEAD OF SUBJ-LEVEL RE
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_famREs.png',width=1500,height=1500)
# prinout random effects
ggplot(predP, aes(x = cbcl_scr_syn_totprob_r, y = fit, group = subjectkey, color=subjectkey))+geom_line()+geom_point(data=famDFvals_select,aes(y=g))+facet_wrap(~ rel_family_id)
dev.off()

