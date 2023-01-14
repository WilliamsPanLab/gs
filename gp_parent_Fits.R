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
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/mixedEfDF_wAdultp.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# could be placed in preproc script
masterdf$interview_age=as.numeric(masterdf$interview_age)
masterdf$subjectkey=as.factor(masterdf$subjectkey)

##############################
######## I SCATTERPLOTS ON TWO VARIABLES OF INTEREST WITH THEIR FIT SPLINE
#### g as response variable
if (!file.exists("/scratch/users/apines/gp/pgAge_parent.rds")){
	print('fitting p~g on masterdf')
	print('master df dims')
	print(dim(masterdf))
	pgAge<-bam(g~s(parentP)+s(interview_age)+s(subjectkey,bs='re'),data=masterdf)
	rdata_file = file("/scratch/users/apines/gp/pgAge_parent.rds", blocking = TRUE)
	saveRDS(pgAge, file=rdata_file)
	close(rdata_file)
} else {
	print(' found pgAge.rds. Loading.')
	pgAge=readRDS('/scratch/users/apines/gp/pgAge_parent.rds')
}
####### plot derivs
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_derivs_pg_parent.png',width=800,height=800)
# prinout gg_derivs
get_derivs_and_plot(pgAge,'parentP')
dev.off()
# use model fit for geom_smooth plotting
forSpline<-predict(pgAge, newdata = masterdf)
print('length of fit:')
print(length(forSpline))
plotdf<-data.frame(masterdf$parentP,forSpline,masterdf$g,as.factor(masterdf$eventname))
colnames(plotdf)<-c('parentP','predicted','g','eventname')
# make scatter plot showing g~total problems with their spline fit
thePlot<-ggplot(plotdf,aes(parentP,g))+
geom_point(alpha=.1,aes(color=eventname))+
geom_smooth(aes(y=predicted),size=2,se=F,color='black')+
theme_classic(base_size=24)+
xlab('parentP')+
scale_color_discrete(name="Vist",breaks=c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),labels=c("Baseline","2 Years"))+
guides(colour = guide_legend(override.aes = list(alpha = 1)))
# make the image to be printed
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/pgAge_Scatter_pg_parent.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()
######## II FORMALLY TEST FOR NON-LINEARITY
#### uses this proposed test https://stats.stackexchange.com/questions/449641/is-there-a-hypothesis-test-that-tells-us-whether-we-should-use-gam-vs-glm
pgAge<-bam(g~s(parentP,m=c(2,0))+s(interview_age)+s(subjectkey,bs='re')+ti(parentP,interview_age),data=masterdf)
summary(pgAge)

print('ti model: parent P')
pgAge<-bam(g~s(parentP)+s(interview_age)+s(subjectkey,bs='re')+ti(parentP,interview_age),data=masterdf)
summary(pgAge)

