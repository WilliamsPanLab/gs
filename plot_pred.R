test=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gIntPredictedBoots.rds')
Intquantiles=apply( test , 2 , quantile , probs = c(0.05,0.95), na.rm = TRUE )
plotdf=data.frame(t(Intquantiles))
colnames(plotdf)<-c('low','high')


# make scatter plot showing g~total problems with their spline fit
thePlot<-ggplot(plotdf,aes(internal,g))+
geom_point(alpha=.1,aes(color=eventname))+
geom_smooth(aes(y=predicted),size=2,se=F,color='black')+
theme_classic(base_size=24)+
xlab('Internalizing Problems Score')+
scale_color_discrete(name="Vist",breaks=c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),labels=c("Baseline","2 Years"))+
guides(colour = guide_legend(override.aes = list(alpha = 1)))
# make the image to be printed
png('/oak/stanford/groups/leanew1/users/apines/figs/gp/IntgAge_Scatter_pg.png',width=800,height=800)
ggMarginal(thePlot,groupFill=T)
dev.off()


