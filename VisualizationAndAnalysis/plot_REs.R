library(mgcv)
library(ggplot2)
library(tidyr)
library(dplyr)

# load in dataframe
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDf.rds')

# load in model
Model=readRDS('/scratch/users/apines/gp/pgAge.rds')

#### WITH G AS RESPONSE

# make "predict" version for plotting RE's
newData <- tidyr::expand(masterdf, nesting(subjectkey), cbcl_scr_syn_totprob_r=unique(cbcl_scr_syn_totprob_r))

# fit at mean age
newData$interview_age=11

# predicted p
predP <- bind_cols(newData, as.data.frame(predict(Model,newdata=newData,se.fit=TRUE)))

# sparsify dataset of interest to remove possiblities of errors
masterdf<-masterdf[,c("subjectkey","g","cbcl_scr_syn_totprob_r")]

# Plot it out
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/gpAge_REs.png',width=13000,height=13000)
# prinout random effects
ggplot(predP, aes(x = cbcl_scr_syn_totprob_r, y = fit, group = subjectkey))+geom_line()+geom_point(data=masterdf,aes(y=g))+facet_wrap(~ subjectkey)
dev.off()

################
