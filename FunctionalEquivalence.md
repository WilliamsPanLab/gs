FunctionalEquivalence
================
2023-02-10

``` r
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
###  ∆∆∆ Demonstration of functionalEquivalence in bcpa g ~ pcag, p ~ p count, and parentP ~parentP count
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
```

``` r
#### LOAD libraries
library(rapportools)
```

    ## 
    ## Attaching package: 'rapportools'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, median, sd, var

    ## The following objects are masked from 'package:base':
    ## 
    ##     max, mean, min, range, sum

``` r
###########∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆##############
### This chunk processes mental health data ###

### LOAD in cbcl data
cbcl=read.delim('~/Downloads/Package_1205735/abcd_cbcl01.txt')
cbcls=read.delim('/Users/panlab/Downloads/Package_1205735/abcd_cbcls01.txt')
# subset timepoints
cbclsBV=subset(cbcls,eventname=='baseline_year_1_arm_1')
cbcls2=subset(cbcls,eventname=='2_year_follow_up_y_arm_1')
# subset timepoints
cbclBV=subset(cbcl,eventname=='baseline_year_1_arm_1')
cbcl2=subset(cbcl,eventname=='2_year_follow_up_y_arm_1')
# merge with other cbcl
cbclsBV=merge(cbclsBV,cbclBV,by=c('subjectkey','eventname'))
cbcls2=merge(cbcls2,cbcl2,by=c('subjectkey','eventname'))

# initialize master df
masterdf<-merge(cbcls,cbcl,by=c('subjectkey','eventname','interview_age','src_subject_id'))
cbcldim<-dim(masterdf)
print(cbcldim)
```

    ## [1] 39767   218

``` r
### LOAD in grades, ∆∆∆ will need to correct for incongruency between tp1 measure (decent granularity) and tp2 measure (high granularity) ∆∆∆
gradesInfoBV=readRDS('~/Downloads/DEAP-data-download-13.rds')
# extract baseline
gradesInfoBV=subset(gradesInfoBV,event_name=='baseline_year_1_arm_1')
gradesInfoBV$Grades<-as.numeric(gradesInfoBV$ksads_back_grades_in_school_p)
# convert ndar value to R
gradesInfoBV$Grades[gradesInfoBV$Grades==-1]=NA
# convert ndar colnames to other ndar colnames
gradesInfoBV$eventname=gradesInfoBV$event_name
gradesInfoBV$subjectkey=gradesInfoBV$src_subject_id
# for tp2, the key is 1 = A's, 2 = B's, 3 = C's, 4 = D's, 5 = F's, -1 = NA
gradesInfoY2=read.delim('~/Downloads/Package_1207225/abcd_saag01.txt')
gradesInfoY2=subset(gradesInfoY2,eventname=='2_year_follow_up_y_arm_1')
gradesInfoY2$sag_grade_type<-as.numeric(gradesInfoY2$sag_grade_type)
# key: 1=100-97,2=96-93,3=92-90,4=89-87,5=86-83,6=82-80,7=79-77,8=76-73,9=72-70,10=69-67,11=66-65,12=0-65,-1=NA,777= no answer
gradesInfoY2$sag_grade_type[gradesInfoY2$sag_grade_type==-1]=NA
gradesInfoY2$sag_grade_type[gradesInfoY2$sag_grade_type==777]=NA
# now convert to be equivalent with timepoint 1 grades measure
ind12=gradesInfoY2$sag_grade_type==12
ind11=gradesInfoY2$sag_grade_type==11
ind10=gradesInfoY2$sag_grade_type==10
ind9=gradesInfoY2$sag_grade_type==9
ind8=gradesInfoY2$sag_grade_type==8
ind7=gradesInfoY2$sag_grade_type==7
ind6=gradesInfoY2$sag_grade_type==6
ind5=gradesInfoY2$sag_grade_type==5
ind4=gradesInfoY2$sag_grade_type==4
ind3=gradesInfoY2$sag_grade_type==3
ind2=gradesInfoY2$sag_grade_type==2
ind1=gradesInfoY2$sag_grade_type==1
#### Set indices to low-res versions
# < 65 becomes failing
gradesInfoY2$sag_grade_type[ind12]=5
# 66-69 = Ds
gradesInfoY2$sag_grade_type[ind11]=4
gradesInfoY2$sag_grade_type[ind10]=4
# 70-79 = Cs
gradesInfoY2$sag_grade_type[ind7]=3
gradesInfoY2$sag_grade_type[ind8]=3
gradesInfoY2$sag_grade_type[ind9]=3
# 80-89 = Bs
gradesInfoY2$sag_grade_type[ind4]=2
gradesInfoY2$sag_grade_type[ind5]=2
gradesInfoY2$sag_grade_type[ind6]=2
# 90+ = As
gradesInfoY2$sag_grade_type[ind1]=1
gradesInfoY2$sag_grade_type[ind2]=1
gradesInfoY2$sag_grade_type[ind3]=1
gradesInfoY2$Grades<-gradesInfoY2$sag_grade_type

###### ∆∆∆∆∆∆∆ create grades info from both of em
NeededColNames=c('subjectkey','eventname','Grades')
gradesInfo<-rbind(gradesInfoBV[,NeededColNames],gradesInfoY2[,NeededColNames])
gradesInfo$Grades<-as.ordered(gradesInfo$Grades)
###### ∆∆∆∆∆∆∆

# merge and count losses
masterdf<-merge(masterdf,gradesInfo,by=c('subjectkey','eventname'))
gradesdim=dim(masterdf)
print(gradesdim)
```

    ## [1] 22289   219

``` r
dif=cbcldim[1]-gradesdim[1]
print(paste0(dif,' rows lost from grades merge, note loss of rows due to no 1 year timepoint'))
```

    ## [1] "17478 rows lost from grades merge, note loss of rows due to no 1 year timepoint"

``` r
### LOAD in ASR data
asr=read.delim('~/Downloads/Package_1207917/pasr01.txt',na.strings=c("","NA"))
masterdf<-merge(masterdf,asr,by=c('subjectkey','eventname','interview_age'))
asrdim=dim(masterdf)
print(asrdim)
```

    ## [1] 22289   364

``` r
dif=gradesdim[1]-asrdim[1]
print(paste0(dif,' rows lost from asr merge'))
```

    ## [1] "0 rows lost from asr merge"

``` r
# load in a DEAP file for rel_family_ID
DEAP=readRDS('~/Downloads/DEAP-data-download-13.rds')
DEAP$subjectkey<-DEAP$src_subject_id
DEAP$eventname=DEAP$event_name
DEAP=DEAP[,c('rel_family_id','subjectkey','eventname')]
masterdf<-merge(masterdf,DEAP,by=c('subjectkey','eventname'))
deapdim=dim(masterdf)
print(deapdim)
```

    ## [1] 22288   365

``` r
dif=asrdim[1]-deapdim[1]
print(paste0(dif,' rows lost from deap familyID merge'))
```

    ## [1] "1 rows lost from deap familyID merge"

``` r
### CLEAN data
# subjectkey as factor
masterdf$subjectkey<-as.factor(masterdf$subjectkey)
# convert cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r<-as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r<-as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r<-as.numeric(masterdf$cbcl_scr_syn_external_r)
# remove instances of NA tot probs
masterdf=masterdf[!is.na(masterdf$cbcl_scr_syn_totprob_r),]
newDim=dim(masterdf)
print(paste0(newDim[1],' after removing NAs for totprob_r, ',(deapdim[1]- newDim[1]),' lost after removing'))
```

    ## [1] "19951 after removing NAs for totprob_r, 2337 lost after removing"

``` r
# and for is empty
masterdf=masterdf[!is.empty(masterdf$cbcl_scr_syn_totprob_r),]
newDim2=dim(masterdf)
print(paste0(newDim2[1],' after removing isempty for totprob_r, ',(newDim[1]- newDim2[1]),' lost after removing'))
```

    ## [1] "18935 after removing isempty for totprob_r, 1016 lost after removing"

``` r
###########∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆##############
###   This chunk processes cognitive data    ###

#### LOAD in cognitive data
nihCog=read.delim('~/Downloads/Package_1206930/abcd_tbss01.txt')
othCog=read.delim('~/Downloads/Package_1206930/abcd_ps01.txt')
littleMan=read.delim('~/Downloads/Package_1206931/lmtp201.txt')

# merge in
masterdf<-merge(masterdf,nihCog,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(masterdf, nihCog, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated
    ## in the result

``` r
newDim3=dim(masterdf)
print(paste0(newDim3[1],' after merging nih toolbox, ',(newDim2[1]- newDim3[1]),' lost after removing'))
```

    ## [1] "18935 after merging nih toolbox, 0 lost after removing"

``` r
masterdf<-merge(masterdf,othCog,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(masterdf, othCog, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'src_subject_id.x', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'src_subject_id.y' are duplicated in the result

``` r
newDim4=dim(masterdf)
print(paste0(newDim4[1],' after merging other cognitive measures, ',(newDim3[1]- newDim4[1]),' lost after removing'))
```

    ## [1] "18935 after merging other cognitive measures, 0 lost after removing"

``` r
masterdf<-merge(masterdf,littleMan,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(masterdf, littleMan, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'src_subject_id.x', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'collection_id.x', 'dataset_id.x', 'src_subject_id.y',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated
    ## in the result

``` r
newDim5=dim(masterdf)
print(paste0(newDim5[1],' after merging little man, ',(newDim4[1]- newDim5[1]),' lost after removing'))
```

    ## [1] "18935 after merging little man, 0 lost after removing"

``` r
# clean age
masterdf$interview_age<-as.numeric(masterdf$interview_age)
masterdf$interview_age<-as.numeric(masterdf$interview_age)/12
```

``` r
###########∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆##############
##This chunk preps for cognition factorization##
###########∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆∆##############

# use thompson 2019 recreation of non nih-tb measures
ind_pea_ravlt = c(which(names(masterdf)=="pea_ravlt_sd_trial_i_tc"),which(names(masterdf)=="pea_ravlt_sd_trial_ii_tc"),
    which(names(masterdf)=="pea_ravlt_sd_trial_iii_tc"),which(names(masterdf)=="pea_ravlt_sd_trial_iv_tc"),
    which(names(masterdf)=="pea_ravlt_sd_trial_v_tc")); names(masterdf)[ind_pea_ravlt];
```

    ## [1] "pea_ravlt_sd_trial_i_tc"   "pea_ravlt_sd_trial_ii_tc" 
    ## [3] "pea_ravlt_sd_trial_iii_tc" "pea_ravlt_sd_trial_iv_tc" 
    ## [5] "pea_ravlt_sd_trial_v_tc"

``` r
# set numbers to numeric
masterdf$pea_ravlt_sd_trial_i_tc=as.numeric(masterdf$pea_ravlt_sd_trial_i_tc)
masterdf$pea_ravlt_sd_trial_ii_tc=as.numeric(masterdf$pea_ravlt_sd_trial_ii_tc)
masterdf$pea_ravlt_sd_trial_iii_tc=as.numeric(masterdf$pea_ravlt_sd_trial_iii_tc)
masterdf$pea_ravlt_sd_trial_iv_tc=as.numeric(masterdf$pea_ravlt_sd_trial_vi_tc)
masterdf$pea_ravlt_sd_trial_v_tc=as.numeric(masterdf$pea_ravlt_sd_trial_v_tc)

# total correct across trials
masterdf$pea_ravlt_ld = masterdf$pea_ravlt_sd_trial_i_tc + masterdf$pea_ravlt_sd_trial_ii_tc + masterdf$pea_ravlt_sd_trial_iii_tc + masterdf$pea_ravlt_sd_trial_iv_tc + masterdf$pea_ravlt_sd_trial_v_tc

# change to numeric
masterdf$nihtbx_picvocab_uncorrected<-as.numeric(masterdf$nihtbx_picvocab_uncorrected)
masterdf$nihtbx_flanker_uncorrected<-as.numeric(masterdf$nihtbx_flanker_uncorrected)
masterdf$nihtbx_list_uncorrected<-as.numeric(masterdf$nihtbx_list_uncorrected)
masterdf$nihtbx_cardsort_uncorrected<-as.numeric(masterdf$nihtbx_cardsort_uncorrected)
masterdf$nihtbx_pattern_uncorrected<-as.numeric(masterdf$nihtbx_pattern_uncorrected)
masterdf$nihtbx_picture_uncorrected<-as.numeric(masterdf$nihtbx_picture_uncorrected)
masterdf$nihtbx_reading_uncorrected<-as.numeric(masterdf$nihtbx_reading_uncorrected)
masterdf$pea_wiscv_tss<-as.numeric(masterdf$pea_wiscv_tss)
masterdf$lmt_scr_perc_correct<-as.numeric(masterdf$lmt_scr_perc_correct)


# for isolating PCA dataframe
pcVars=c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","pea_ravlt_ld","lmt_scr_perc_correct")

# test for completeness before running PCA. Better move to calculate ONLY in the sample that we are running analyses on (more technically accurate than PC structure slightly misaligned with sample of interest)
# get other vars of interest to check for complete cases
KidVarsOfInt=c('Grades','cbcl_scr_syn_totprob_r','cbcl_scr_syn_external_r','cbcl_scr_syn_internal_r')
# asr columns of interest to gauge completeness of
ColsOfInt=asr[,c(11:141)]
ASRVarsOfInt=colnames(ColsOfInt)

# only use subjects with both timepoints as complete cases
subjs=unique(masterdf$subjectkey)
for (s in subjs){
  # if there are less than two complete cases of the variables of interest
  if (sum(complete.cases(masterdf[masterdf$subjectkey==s,c(pcVars,KidVarsOfInt,ASRVarsOfInt)]))<2){
    subjs=subjs[subjs!=s]
  }
}
# convert masterdf to df with complete observations for cognition
masterdf=masterdf[masterdf$subjectkey %in% subjs,]

newDim6=dim(masterdf)
print(paste0(newDim6[1],' after retaining only subjs with vars of int at BOTH timepoints, ',(newDim5[1]- newDim6[1]),' lost after removing'))
```

    ## [1] "11452 after retaining only subjs with vars of int at BOTH timepoints, 7483 lost after removing"

``` r
print(dim(masterdf))
```

    ## [1] 11452   597

``` r
# finish cleaning data for sherlock runs: one family member per family to facilitate random sample
masterdf$id_fam = NULL
# default value of family size (# of children in abcd study)
masterdf$fam_size = 1

# counter index
ind=0

# set each instance of multiple family members to a family ID, as ind
set.seed(1)
for(f in 1:length(unique(masterdf$rel_family_id))){
  # calculate family size
  famsize=sum(masterdf$rel_family_id == unique(masterdf$rel_family_id)[f]) / 2
  masterdf$fam_size[masterdf$rel_family_id == unique(masterdf$rel_family_id)[f]] = famsize
  # note that each  person is represented twice at this point:
  # divide by 2 to take number of visits to number of people, if there's more than 2x visits per family ID, izza family
  # this logic gets hairy. Starting from outside in: > 1 is family size >1, /2 is divided by 2 for two visits, [f] is unique familyID, rel_family_id is place in column of masterdf
  if(famsize>1){
    # remove one from instances where family-id = this relative family id (sequence for siblings, 1:size(Family))
    #print(paste0('family size ',famsize))
    # keep one sib
    kept=sample(seq(1,famsize),1)
    #print(paste0('kept ',kept))
    # use to select one subject id
    famIDs=unique(masterdf$subjectkey[masterdf$rel_family_id == unique(masterdf$rel_family_id)[f]])
    # chosen sib
    keeper=famIDs[kept]
    left=famIDs[-c(kept)]
    # leave rest
    masterdf=masterdf[masterdf$subjectkey!=left,] 
    #print(paste0('left ',left))
    # calc index of family
    ind=ind+1   
    # set index of family
    masterdf$id_fam[masterdf$rel_family_id == unique(masterdf$rel_family_id)[f]] = ind
  } 
}
```

    ## Warning in `!=.default`(masterdf$subjectkey, left): longer object length is not
    ## a multiple of shorter object length

    ## Warning in is.na(e1) | is.na(e2): longer object length is not a multiple of
    ## shorter object length

    ## Warning in `!=.default`(masterdf$subjectkey, left): longer object length is not
    ## a multiple of shorter object length

    ## Warning in is.na(e1) | is.na(e2): longer object length is not a multiple of
    ## shorter object length

``` r
# make family ID for those with families represented in ABCD
masterdf$rel_family_id=masterdf$id_fam

newDim7=dim(masterdf)
print(paste0(newDim7[1],' after retaining only one subjs per family, ',(newDim6[1]- newDim7[1]),' lost after removing'))
```

    ## [1] "10101 after retaining only one subjs per family, 1351 lost after removing"

``` r
#       NOW 
# THAT'S WHAT I CALL PCAPREP
#       271

pcVars=c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","pea_ravlt_ld","lmt_scr_perc_correct","pea_wiscv_tss","nihtbx_list_uncorrected","nihtbx_cardsort_uncorrected")
# only pca vars, only timepoint 1 for eval of overlap with bpca
masterdfbv=masterdf[masterdf$eventname=="baseline_year_1_arm_1",]
masterdfbv=masterdfbv[complete.cases(masterdfbv[,pcVars]),]
# and for a pcadf version (pure numeric)
pcaDf<-masterdfbv[complete.cases(masterdfbv[,pcVars]),pcVars]
```

``` r
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
#### G: BPCA ~= PCA
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 

# compare to thompson PC1
PC1=readRDS('/Users/panlab/Downloads/DEAP-data-download.rds')
PC1$subjectkey=PC1$src_subject_id
PC1$eventname<-PC1$event_name
# only baseline is meaningful
PC1<-PC1[PC1$event_name=='baseline_year_1_arm_1',]

# derive pcs
Y = as.matrix(scale(pcaDf[,pcVars]))
# equiv for binding scores to IDs and eventnames
pcVarsAndIDs=c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","pea_ravlt_ld","lmt_scr_perc_correct","subjectkey","eventname")
Yextended=masterdfbv[,pcVarsAndIDs]
ncomp = 3
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
y.pca$loadings
```

    ## 
    ## Loadings:
    ##                             RC1   RC2   RC3  
    ## nihtbx_picvocab_uncorrected 0.765       0.171
    ## nihtbx_flanker_uncorrected  0.210 0.722      
    ## nihtbx_pattern_uncorrected        0.795 0.102
    ## nihtbx_picture_uncorrected        0.137 0.856
    ## nihtbx_reading_uncorrected  0.815 0.116      
    ## pea_ravlt_ld                0.326 0.148 0.700
    ## lmt_scr_perc_correct        0.503 0.294      
    ## pea_wiscv_tss               0.535       0.379
    ## nihtbx_list_uncorrected     0.524 0.195 0.418
    ## nihtbx_cardsort_uncorrected 0.217 0.710 0.219
    ## 
    ##                  RC1   RC2   RC3
    ## SS loadings    2.264 1.850 1.642
    ## Proportion Var 0.226 0.185 0.164
    ## Cumulative Var 0.226 0.411 0.576

``` r
# assign scores to subjs
Yextended$g<-y.pca$scores[,1]
# merge in cog data
masterdfbv$g<-Yextended$g
# mergin'
compareG_df<-merge(PC1,masterdfbv,by=c('src_subject_id'))
```

    ## Warning in merge.data.frame(PC1, masterdfbv, by = c("src_subject_id")):
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.x',
    ## 'dataset_id.x', 'src_subject_id.y', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y' are duplicated in the result

``` r
# compare
compareG_df$g<-as.numeric(compareG_df$g)
compareG_df$neurocog_pc1.bl<-as.numeric(compareG_df$neurocog_pc1.bl)
print(cor.test(compareG_df$g,compareG_df$neurocog_pc1.bl))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  compareG_df$g and compareG_df$neurocog_pc1.bl
    ## t = 250.52, df = 4973, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9604948 0.9645782
    ## sample estimates:
    ##       cor 
    ## 0.9625911

``` r
library(ggplot2)
plotdf<-compareG_df[,c('g','neurocog_pc1.bl')]
ComparePlot<-ggplot(data=plotdf,aes(x=g,y=neurocog_pc1.bl))+geom_smooth(method='lm',color='black')+geom_point(alpha=.05)+theme_classic()
print(ComparePlot)
```

    ## `geom_smooth()` using formula 'y ~ x'

    ## Warning: Removed 14 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 14 rows containing missing values (geom_point).

![](FunctionalEquivalence_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
### Can add cormat among loadings for all 3 components it looks like they are functionally equivalent all the way down
```

``` r
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
#### P: P FACTOR  ~= SYMPTOM COUNT
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
## write a version parallel to michelini et al 2019
#### low freq removal
##The following CBCL items were removed because of low frequency: “Drinks alcohol without parents' approval”, “Sexual problems”, “Smokes, chews, or sniffs tobacco”, “Truancy, skips school”, “Uses drugs for non-medical #purposes (don't include alcohol or tobacco)”. 
## find items
lf1<-grep('Drinks alcohol without parents',cbcl[1,])
lf2<-grep('Sexual problem',cbcl[1,])
lf3<-grep('Smokes, chews, or sniffs tobacco',cbcl[1,])
lf4<-grep('Truancy, skips school',cbcl[1,])
lf5<-grep('non medical purpose',cbcl[1,])
# check if these meet criteria in tp2
lf1_nonzero=sum(as.numeric(cbcl2[,lf1])>0)
lf2_nonzero=sum(as.numeric(cbcl2[,lf2])>0)
lf3_nonzero=sum(as.numeric(cbcl2[,lf3])>0)
lf4_nonzero=sum(as.numeric(cbcl2[,lf4])>0)
lf5_nonzero=sum(as.numeric(cbcl2[,lf5])>0)
# add tp1 values
lf1_nonzero=lf1_nonzero+sum(as.numeric(cbclBV[,lf1])>0)
lf2_nonzero=lf2_nonzero+sum(as.numeric(cbclBV[,lf2])>0)
lf3_nonzero=lf3_nonzero+sum(as.numeric(cbclBV[,lf3])>0)
lf4_nonzero=lf4_nonzero+sum(as.numeric(cbclBV[,lf4])>0)
lf5_nonzero=lf5_nonzero+sum(as.numeric(cbclBV[,lf5])>0)

# num rows
numRows=dim(cbcl2)[1]+dim(cbclBV)[1]

# is prop greater 
lf1_nonzero/numRows
```

    ## [1] 0.001390758

``` r
lf2_nonzero/numRows
```

    ## [1] 0.002422611

``` r
lf3_nonzero/numRows
```

    ## [1] 0.0004486317

``` r
lf4_nonzero/numRows
```

    ## [1] 0.006191117

``` r
lf5_nonzero/numRows
```

    ## [1] 0.000807537

``` r
# looks like truancy meets their inclusion criteria with tp2 included

# remove those meeting criteria
cbclBV=cbclBV[-c(lf1,lf2,lf3,lf5)]
cbcl2=cbcl2[-c(lf1,lf2,lf3,lf5)]
# and col names to index later
CBCLcolnames=cbcl[1,]
CBCLcolnames=CBCLcolnames[-c(lf1,lf2,lf3,lf5)]
```

``` r
# polycors to eval if tp1 vars still meet criterion
library(polycor)

# isolate cbcl qs
cbclBV_qs=cbclBV[,10:124]
# fix character specification for variables
cbclBV_qs <- as.data.frame(lapply(cbclBV_qs, as.ordered))
# polycormat
cbclBV_qs_cormat=suppressWarnings(hetcor(cbclBV_qs))
# ok piecemeal go through and find rows with >.75 cors
for (i in 1:dim(cbclBV_qs_cormat$correlations)[1]){
  # 1 is for diagonal
  if (sum(cbclBV_qs_cormat$correlations[i,]>.75)>1){
    print('Correlated Items:')
    # + 9 because 10th col is first col
    print(CBCLcolnames[1,i+9])
    BVoverCorrelateds<-colnames(CBCLcolnames[1,which(cbclBV_qs_cormat$correlations[i,]>.75)+9])
    print(unlist(CBCLcolnames[BVoverCorrelateds]))
    print('-------')
    }
}
```

    ## [1] "Correlated Items:"
    ## [1] "Can't concentrate, can't pay attention for long No puede concentrarse o prestar atención por mucho tiempo"
    ##                                                                                                  cbcl_q08_p 
    ## "Can't concentrate, can't pay attention for long No puede concentrarse o prestar atención por mucho tiempo" 
    ##                                                                                                  cbcl_q10_p 
    ##    "Can't sit still, restless, or hyperactive No puede quedarse quieto(a); es inquieto(a) o hiperactivo(a)" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Can't sit still, restless, or hyperactive No puede quedarse quieto(a); es inquieto(a) o hiperactivo(a)"
    ##                                                                                                  cbcl_q08_p 
    ## "Can't concentrate, can't pay attention for long No puede concentrarse o prestar atención por mucho tiempo" 
    ##                                                                                                  cbcl_q10_p 
    ##    "Can't sit still, restless, or hyperactive No puede quedarse quieto(a); es inquieto(a) o hiperactivo(a)" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Cruelty, bullying, or meanness to others Es cruel, abusador(a), y malo(a) con los demás"
    ##                                                                                cbcl_q16_p 
    ## "Cruelty, bullying, or meanness to others Es cruel, abusador(a), y malo(a) con los demás" 
    ##                                                                                cbcl_q97_p 
    ##                                                        "Threatens people Amenaza a otros" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Destroys their own things"
    ##                                                                 cbcl_q20_p 
    ##                                                "Destroys their own things" 
    ##                                                                cbcl_q106_p 
    ## "Vandalism Comete actos de vandalismo, como romper ventanas u otras cosas" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Disobedient at home  Desobedece en casa"
    ##                                                                                                       cbcl_q22_p 
    ##                                                                        "Disobedient at home  Desobedece en casa" 
    ##                                                                                                       cbcl_q28_p 
    ## "Breaks rules at home, school or elsewhere  No respeta/rompe las reglas en casa, en la escuela, o en otro lugar" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Disobedient at school Desobedece en la escuela"
    ##                                                                                                       cbcl_q23_p 
    ##                                                                 "Disobedient at school Desobedece en la escuela" 
    ##                                                                                                       cbcl_q28_p 
    ## "Breaks rules at home, school or elsewhere  No respeta/rompe las reglas en casa, en la escuela, o en otro lugar" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Doesn't get along with other kids No se lleva bien con otros niños(as)/jóvenes"
    ##                                                                       cbcl_q25_p 
    ## "Doesn't get along with other kids No se lleva bien con otros niños(as)/jóvenes" 
    ##                                                                       cbcl_q48_p 
    ##               "Not liked by other kids No le cae bien a otros niños(as)/jóvenes" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Breaks rules at home, school or elsewhere  No respeta/rompe las reglas en casa, en la escuela, o en otro lugar"
    ##                                                                                                       cbcl_q22_p 
    ##                                                                        "Disobedient at home  Desobedece en casa" 
    ##                                                                                                       cbcl_q23_p 
    ##                                                                 "Disobedient at school Desobedece en la escuela" 
    ##                                                                                                       cbcl_q28_p 
    ## "Breaks rules at home, school or elsewhere  No respeta/rompe las reglas en casa, en la escuela, o en otro lugar" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Gets teased a lot Los demás se burlan de él/ella a menudo"
    ##                                                         cbcl_q38_p 
    ##        "Gets teased a lot Los demás se burlan de él/ella a menudo" 
    ##                                                         cbcl_q48_p 
    ## "Not liked by other kids No le cae bien a otros niños(as)/jóvenes" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Not liked by other kids No le cae bien a otros niños(as)/jóvenes"
    ##                                                                       cbcl_q25_p 
    ## "Doesn't get along with other kids No se lleva bien con otros niños(as)/jóvenes" 
    ##                                                                       cbcl_q38_p 
    ##                      "Gets teased a lot Los demás se burlan de él/ella a menudo" 
    ##                                                                       cbcl_q48_p 
    ##               "Not liked by other kids No le cae bien a otros niños(as)/jóvenes" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Nausea, feels sick Náuseas, ganas de vomitar"
    ##                                    cbcl_q56c_p 
    ## "Nausea, feels sick Náuseas, ganas de vomitar" 
    ##                                    cbcl_q56f_p 
    ##             "Stomachaches Dolores de estómago" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Stomachaches Dolores de estómago"
    ##                                    cbcl_q56c_p 
    ## "Nausea, feels sick Náuseas, ganas de vomitar" 
    ##                                    cbcl_q56f_p 
    ##             "Stomachaches Dolores de estómago" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Threatens people Amenaza a otros"
    ##                                                                                cbcl_q16_p 
    ## "Cruelty, bullying, or meanness to others Es cruel, abusador(a), y malo(a) con los demás" 
    ##                                                                                cbcl_q97_p 
    ##                                                        "Threatens people Amenaza a otros" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ## [1] "Vandalism Comete actos de vandalismo, como romper ventanas u otras cosas"
    ##                                                                 cbcl_q20_p 
    ##                                                "Destroys their own things" 
    ##                                                                cbcl_q106_p 
    ## "Vandalism Comete actos de vandalismo, como romper ventanas u otras cosas" 
    ## [1] "-------"

``` r
#### composite creation
####The following composites were created: Attacks/threatens (“Physically attacks people”, “Threatens people”); Destroys (“Destroys his/her own things”, “Destroys things belonging to his/her family or others”, “Vandalism”); Disobeys rules (“Disobedient at home”, “Disobedient at school”, “Breaks rules at home, school or elsewhere”); Steals (“Steals at home”, “Steals outside ###the home”); Peer problems (“Doesn't get along with other kids”, “Not liked by other kids”); Distracted/Hyperactive (“Can't concentrate, can't pay attention for long”, “Inattentive or easily distracted”, “Can't sit still, restless, or hyperactive”); Hallucinations (“Hears sound or voices that aren't there”, “Sees things that aren't there”); Sex play (“Plays with own sex ###parts in public”, “Plays with own sex parts too much”); Weight problems (“Overeating”, “Overweight”)

# merge timepoints
mergedTPs=cbclBV
# retain subjs and 
mergedTPsSubjs<-mergedTPs$subjectkey
mergedTPsEventName<-mergedTPs$eventname
mergedTPs <- as.data.frame(lapply(mergedTPs, as.numeric))
```

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

    ## Warning in lapply(mergedTPs, as.numeric): NAs introduced by coercion

``` r
# inattentive composite
mergedTPs$cbcl_inattentive=round((mergedTPs$cbcl_q08_p+mergedTPs$cbcl_q10_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`cbcl_q08_p`,`cbcl_q10_p`))
# aggressive composite
mergedTPs$cbcl_aggresive=round((mergedTPs$cbcl_q16_p+mergedTPs$cbcl_q97_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`cbcl_q16_p`,`cbcl_q97_p`))
# anarchic composite
mergedTPs$cbcl_anarchic=round((mergedTPs$cbcl_q22_p+mergedTPs$cbcl_q23_p+mergedTPs$cbcl_q28_p)/3)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`cbcl_q22_p`,`cbcl_q23_p`,`cbcl_q28_p`))
# unpopular composite
mergedTPs$cbcl_unpopular=round((mergedTPs$cbcl_q25_p+mergedTPs$cbcl_q48_p+mergedTPs$cbcl_q38_p)/3)
# remove constits
mergedTPs=subset(mergedTPs, select = -c(`cbcl_q25_p`,`cbcl_q48_p`,`cbcl_q38_p`))
# tummy composite
mergedTPs$cbcl_tummy=round((mergedTPs$cbcl_q56f_p+mergedTPs$cbcl_q56c_p)/2)

####### ∆∆∆∆∆∆∆∆∆∆
## run pca on cbcl
pcaDf_p<-mergedTPs[,10:122]
pcaDf_p$SubjsNames<-mergedTPsSubjs
pcaDf_p$EventNames<-mergedTPsEventName
pcaDf_p=subset(pcaDf_p, select = -c(`eventname`,`timept`,`collection_title`))
# remove NA vars
pcaDf_p_Complete<-pcaDf_p[complete.cases(pcaDf_p),]
pcaDf_p_CompleteSubjs<-pcaDf_p_Complete$SubjsNames
pcaDf_p_CompleteEventNames<-pcaDf_p_Complete$EventNames
# isolate numeric
pcaDf_p_num<-pcaDf_p_Complete[,1:110]
# convert to numeric for pca
pcaDf_p_num <- as.data.frame(lapply(pcaDf_p_num, as.numeric))
# derive pcs
ncomp = 1
y.pca = psych::principal(pcaDf_p_num, rotate="geominT", nfactors=ncomp, scores=TRUE)
y.pca$loadings
```

    ## 
    ## Loadings:
    ##                  PC1  
    ## cbcl_q01_p       0.462
    ## cbcl_q03_p       0.592
    ## cbcl_q04_p       0.571
    ## cbcl_q05_p       0.444
    ## cbcl_q06_p       0.140
    ## cbcl_q07_p       0.390
    ## cbcl_q09_p       0.582
    ## cbcl_q11_p       0.462
    ## cbcl_q12_p       0.491
    ## cbcl_q13_p       0.420
    ## cbcl_q14_p       0.451
    ## cbcl_q15_p       0.218
    ## cbcl_q17_p       0.443
    ## cbcl_q18_p       0.223
    ## cbcl_q19_p       0.589
    ## cbcl_q20_p       0.524
    ## cbcl_q21_p       0.535
    ## cbcl_q24_p       0.346
    ## cbcl_q26_p       0.487
    ## cbcl_q27_p       0.539
    ## cbcl_q29_p       0.355
    ## cbcl_q30_p       0.337
    ## cbcl_q31_p       0.366
    ## cbcl_q32_p       0.304
    ## cbcl_q33_p       0.529
    ## cbcl_q34_p       0.453
    ## cbcl_q35_p       0.510
    ## cbcl_q36_p       0.361
    ## cbcl_q37_p       0.402
    ## cbcl_q39_p       0.341
    ## cbcl_q40_p       0.224
    ## cbcl_q41_p       0.638
    ## cbcl_q42_p       0.385
    ## cbcl_q43_p       0.526
    ## cbcl_q44_p       0.250
    ## cbcl_q45_p       0.550
    ## cbcl_q46_p       0.409
    ## cbcl_q47_p       0.384
    ## cbcl_q49_p       0.265
    ## cbcl_q50_p       0.494
    ## cbcl_q51_p       0.276
    ## cbcl_q52_p       0.365
    ## cbcl_q53_p       0.323
    ## cbcl_q54_p       0.367
    ## cbcl_q55_p       0.186
    ## cbcl_q56a_p      0.283
    ## cbcl_q56b_p      0.287
    ## cbcl_q56c_p      0.338
    ## cbcl_q56d_p      0.176
    ## cbcl_q56e_p      0.216
    ## cbcl_q56f_p      0.328
    ## cbcl_q56g_p      0.181
    ## cbcl_q56h_p      0.181
    ## cbcl_q57_p       0.428
    ## cbcl_q58_p       0.377
    ## cbcl_q59_p       0.144
    ## cbcl_q60_p       0.163
    ## cbcl_q61_p       0.473
    ## cbcl_q62_p       0.439
    ## cbcl_q63_p       0.363
    ## cbcl_q64_p       0.388
    ## cbcl_q65_p       0.402
    ## cbcl_q66_p       0.464
    ## cbcl_q67_p       0.234
    ## cbcl_q68_p       0.543
    ## cbcl_q69_p       0.456
    ## cbcl_q70_p       0.220
    ## cbcl_q71_p       0.483
    ## cbcl_q72_p       0.154
    ## cbcl_q74_p       0.458
    ## cbcl_q75_p       0.309
    ## cbcl_q76_p       0.363
    ## cbcl_q77_p       0.223
    ## cbcl_q78_p       0.606
    ## cbcl_q79_p       0.208
    ## cbcl_q80_p       0.458
    ## cbcl_q81_p       0.385
    ## cbcl_q82_p       0.323
    ## cbcl_q83_p       0.400
    ## cbcl_q84_p       0.475
    ## cbcl_q85_p       0.443
    ## cbcl_q86_p       0.604
    ## cbcl_q87_p       0.633
    ## cbcl_q88_p       0.547
    ## cbcl_q89_p       0.453
    ## cbcl_q90_p       0.397
    ## cbcl_q91_p       0.322
    ## cbcl_q92_p       0.218
    ## cbcl_q93_p       0.455
    ## cbcl_q94_p       0.450
    ## cbcl_q95_p       0.602
    ## cbcl_q96_p       0.213
    ## cbcl_q98_p            
    ## cbcl_q100_p      0.427
    ## cbcl_q101_p      0.159
    ## cbcl_q102_p      0.383
    ## cbcl_q103_p      0.539
    ## cbcl_q104_p      0.527
    ## cbcl_q106_p      0.299
    ## cbcl_q107_p      0.160
    ## cbcl_q108_p      0.153
    ## cbcl_q109_p      0.501
    ## cbcl_q110_p           
    ## cbcl_q111_p      0.421
    ## cbcl_q112_p      0.489
    ## cbcl_inattentive 0.591
    ## cbcl_aggresive   0.408
    ## cbcl_anarchic    0.602
    ## cbcl_unpopular   0.541
    ## cbcl_tummy       0.344
    ## 
    ##                   PC1
    ## SS loadings    18.466
    ## Proportion Var  0.168

``` r
# assign scores to subjs
subjPvalues=data.frame(pcaDf_p_CompleteSubjs,pcaDf_p_CompleteEventNames,y.pca$scores[,1])
colnames(subjPvalues)<-c('subjectkey','eventname','p')

# compare 'em!
testEquivDf=merge(masterdf,subjPvalues,by=c('subjectkey','eventname'))
```

    ## Warning in merge.data.frame(masterdf, subjPvalues, by = c("subjectkey", :
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.x',
    ## 'dataset_id.x', 'src_subject_id.y', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y' are duplicated in the result

``` r
print(cor.test(testEquivDf$p,testEquivDf$cbcl_scr_syn_totprob_r))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  testEquivDf$p and testEquivDf$cbcl_scr_syn_totprob_r
    ## t = 615.2, df = 5052, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9930184 0.9937453
    ## sample estimates:
    ##       cor 
    ## 0.9933918

``` r
plotdf<-testEquivDf[,c('p','cbcl_scr_syn_totprob_r')]
ComparePlot<-ggplot(data=plotdf,aes(x=p,y=cbcl_scr_syn_totprob_r))+geom_smooth(method='lm',col='black')+geom_point(alpha=.05)+theme_classic()
print(ComparePlot)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](FunctionalEquivalence_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
#### PARENT P: COUNT ~= FACTOR
################# ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ ∆∆∆ 
### derive adult p
# columns of interest to gauge completeness of
ColsOfInt=asr[,c(11:141)]
# retain complete cases
completeInd=ColsOfInt[rowSums(is.na(ColsOfInt)) == 0,]
ASRcolnames=asr[1,]

### following Micheli approach of composite creation

# polycors to eval if tp1 vars still meet criterion
library(polycor)
# subset asr into timepoints
asrBV=subset(asr,eventname=='baseline_year_1_arm_1')
asr2=subset(asr,eventname=='2_year_follow_up_y_arm_1')

# save subjIDs
asrBVSubjs=asrBV$subjectkey
asr2Subjs=asr2$subjectkey
# isolate asr qs
asrBV_qs=asrBV[,11:141]
asr2_qs=asr2[,11:141]
# fix character specification for variables
asrBV_qs <- as.data.frame(lapply(asrBV_qs, as.ordered))
asr2_qs <- as.data.frame(lapply(asr2_qs, as.ordered))
# polycormat - supress wave of warnings
asrBV_qs_cormat=suppressWarnings(hetcor(asrBV_qs))
asr2_qs_cormat=suppressWarnings(hetcor(asr2_qs))
# ok piecemeal go through and find rows with >.75 cors
for (i in 1:dim(asrBV_qs_cormat$correlations)[1]){
  # a single 1 is expected for diagonal
  if (sum(asrBV_qs_cormat$correlations[i,]>.75)>1){
    # if it is > .75 in both tps
    if (sum(asr2_qs_cormat$correlations[i,]>.75)>1){
      print('Correlated Items:')
      # + 10 because 11th col is first col
      BVoverCorrelateds<-colnames(ASRcolnames[1,which(asrBV_qs_cormat$correlations[i,]>.75)+10])
      year2overCorrelateds<-colnames(ASRcolnames[1,which(asr2_qs_cormat$correlations[i,]>.75)+10])
      intersection=intersect(BVoverCorrelateds,year2overCorrelateds)
      print(unlist(ASRcolnames[intersection]))
      print('-------')
      }
    }
}
```

    ## [1] "Correlated Items:"
    ##                                                                             asr_q20_p 
    ##                            "I damage or destroy my things Destruyo mis propias cosas" 
    ##                                                                             asr_q21_p 
    ## "I damage or destroy things belonging to others Destruyo las cosas de otras personas" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                             asr_q20_p 
    ##                            "I damage or destroy my things Destruyo mis propias cosas" 
    ##                                                                             asr_q21_p 
    ## "I damage or destroy things belonging to others Destruyo las cosas de otras personas" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                                     asr_q36_p 
    ## "I accidentally get hurt a lot, accident-prone Me lastimo accidentalmente con mucha frecuencia, soy propenso(a) a accidentes" 
    ##                                                                                                                     asr_q62_p 
    ##                                                         "I am poorly coordinated or clumsy Tengo mala coordinación o torpeza" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                           asr_q40_p 
    ## "I hear sounds and voices that other people think aren't there Oigo sonidos o voces que otros creen que no existen" 
    ##                                                                                                           asr_q70_p 
    ##                        "I see things that other people think aren't there Veo cosas que otros creen que no existen" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                           asr_q45_p 
    ##                 "I am nervous or tense Soy nervioso(a), o tenso(a)" 
    ##                                                           asr_q50_p 
    ## "I am too fearful or anxious Soy demasiado miedoso(a) o ansioso(a)" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                           asr_q45_p 
    ##                 "I am nervous or tense Soy nervioso(a), o tenso(a)" 
    ##                                                           asr_q50_p 
    ## "I am too fearful or anxious Soy demasiado miedoso(a) o ansioso(a)" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                 asr_q55_p 
    ## "My moods swing between elation and depression Mi humor cambia entre euforia y depresión" 
    ##                                                                                 asr_q87_p 
    ##       "My moods or feeling change suddenly Tengo bruscos cambios de humor o sentimientos" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                     asr_q56c_p                     asr_q56g_p 
    ##    "Nausea, feel sick Náuseas" "Vomiting, throwing up Vómito" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                     asr_q56c_p                     asr_q56g_p 
    ##    "Nausea, feel sick Náuseas" "Vomiting, throwing up Vómito" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                                     asr_q36_p 
    ## "I accidentally get hurt a lot, accident-prone Me lastimo accidentalmente con mucha frecuencia, soy propenso(a) a accidentes" 
    ##                                                                                                                     asr_q62_p 
    ##                                                         "I am poorly coordinated or clumsy Tengo mala coordinación o torpeza" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                           asr_q40_p 
    ## "I hear sounds and voices that other people think aren't there Oigo sonidos o voces que otros creen que no existen" 
    ##                                                                                                           asr_q70_p 
    ##                        "I see things that other people think aren't there Veo cosas que otros creen que no existen" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                             asr_q84_p 
    ##              "I do things that other people think are strange Hago cosas que otras personas piensan que son extrañas" 
    ##                                                                                                             asr_q85_p 
    ## "I have thoughts that other people would think are strange Tengo ideas que otras personas pensarían que son extrañas" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                             asr_q84_p 
    ##              "I do things that other people think are strange Hago cosas que otras personas piensan que son extrañas" 
    ##                                                                                                             asr_q85_p 
    ## "I have thoughts that other people would think are strange Tengo ideas que otras personas pensarían que son extrañas" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                 asr_q55_p 
    ## "My moods swing between elation and depression Mi humor cambia entre euforia y depresión" 
    ##                                                                                 asr_q87_p 
    ##       "My moods or feeling change suddenly Tengo bruscos cambios de humor o sentimientos" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                                             asr_q114_p 
    ## "I fail to pay my debts or meet other financial responsibilities No pago mis deudas ni me hago cargo de responsabilidades financieras" 
    ##                                                                                                                             asr_q117_p 
    ##                        "I have trouble managing my money or credit card Me cuesta trabajo manejar el dinero o las tarjetas de crédito" 
    ## [1] "-------"
    ## [1] "Correlated Items:"
    ##                                                                                                                             asr_q114_p 
    ## "I fail to pay my debts or meet other financial responsibilities No pago mis deudas ni me hago cargo de responsabilidades financieras" 
    ##                                                                                                                             asr_q117_p 
    ##                        "I have trouble managing my money or credit card Me cuesta trabajo manejar el dinero o las tarjetas de crédito" 
    ## [1] "-------"

``` r
#### merge timepoints
mergedTPs=rbind(asrBV_qs,asr2_qs)
# retain subjsIDs and 
mergedTPsSubjs<-c(asrBVSubjs,asr2Subjs)
mergedTPsEventName<-rep()
mergedTPs <- as.data.frame(lapply(mergedTPs, as.numeric))

# destroyer composite
mergedTPs$asr_destroyer=round((mergedTPs$asr_q20_p+mergedTPs$asr_q21_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q20_p`,`asr_q21_p`))

# hallucinations composite
mergedTPs$asr_halluc=round((mergedTPs$asr_q40_p+mergedTPs$asr_q70_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q40_p`,`asr_q70_p`))

# Odd composite
mergedTPs$asr_odd=round((mergedTPs$asr_q84_p+mergedTPs$asr_q85_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q84_p`,`asr_q85_p`))

# sick composite
mergedTPs$asr_sick=round((mergedTPs$asr_q56c_p+mergedTPs$asr_q56g_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q56c_p`,`asr_q56g_p`))

# clumsy composite
mergedTPs$asr_clumsy=round((mergedTPs$asr_q36_p+mergedTPs$asr_q62_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q36_p`,`asr_q62_p`))

# anxious composite
mergedTPs$asr_anx=round((mergedTPs$asr_q45_p+mergedTPs$asr_q50_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q45_p`,`asr_q50_p`))

# swinger composite
mergedTPs$asr_swing=round((mergedTPs$asr_q55_p+mergedTPs$asr_q87_p)/2)
# remove constituent variables
mergedTPs=subset(mergedTPs, select = -c(`asr_q55_p`,`asr_q87_p`))

# insolvent composite
mergedTPs$asr_insolvent=round((mergedTPs$asr_q114_p+mergedTPs$asr_q117_p)/2)
# remove constits
mergedTPs=subset(mergedTPs, select = -c(`asr_q114_p`,`asr_q117_p`))
#####################

## run pca on asr
pcaDf_p<-mergedTPs
pcaDf_p$SubjsNames<-mergedTPsSubjs
pcaDf_p$eventname<-c(rep('baseline_year_1_arm_1',dim(asrBV)[1]),rep('2_year_follow_up_y_arm_1',dim(asr2)[1]))
# remove NA vars
pcaDf_p_Complete<-pcaDf_p[complete.cases(pcaDf_p),]
pcaDf_p_CompleteSubjs<-pcaDf_p_Complete$SubjsNames
pcaDf_p_CompleteEventNames<-pcaDf_p_Complete$eventname
pcaDf_p=subset(pcaDf_p, select = -c(`eventname`,`SubjsNames`))
# isolate numeric
pcaDf_p_num<-pcaDf_p_Complete[,1:123]
# convert to numeric for pca
pcaDf_p_num <- as.data.frame(lapply(pcaDf_p_num, as.numeric))
# derive pcs
pcaMat_p_complete = as.matrix(scale(pcaDf_p_num))
ncomp = 1
y.pca = psych::principal(pcaMat_p_complete, rotate="geomin", nfactors=ncomp, scores=TRUE)
y.pca$loadings
```

    ## 
    ## Loadings:
    ##               PC1   
    ## asr_q01_p      0.473
    ## asr_q02_p     -0.204
    ## asr_q03_p      0.401
    ## asr_q04_p     -0.170
    ## asr_q05_p      0.333
    ## asr_q06_p      0.152
    ## asr_q07_p      0.210
    ## asr_q08_p      0.528
    ## asr_q09_p      0.571
    ## asr_q10_p      0.389
    ## asr_q11_p      0.406
    ## asr_q12_p      0.567
    ## asr_q13_p      0.561
    ## asr_q14_p      0.465
    ## asr_q15_p           
    ## asr_q16_p      0.321
    ## asr_q17_p      0.390
    ## asr_q18_p      0.193
    ## asr_q19_p      0.235
    ## asr_q22_p      0.472
    ## asr_q23_p      0.308
    ## asr_q24_p      0.415
    ## asr_q25_p      0.393
    ## asr_q26_p      0.176
    ## asr_q27_p      0.400
    ## asr_q28_p      0.409
    ## asr_q29_p      0.319
    ## asr_q30_p      0.406
    ## asr_q31_p      0.396
    ## asr_q32_p      0.360
    ## asr_q33_p      0.492
    ## asr_q34_p      0.405
    ## asr_q35_p      0.589
    ## asr_q37_p      0.243
    ## asr_q38_p      0.298
    ## asr_q39_p      0.219
    ## asr_q41_p      0.503
    ## asr_q42_p      0.438
    ## asr_q43_p      0.325
    ## asr_q44_p      0.546
    ## asr_q46_p      0.466
    ## asr_q47_p      0.573
    ## asr_q48_p      0.451
    ## asr_q49_p      0.214
    ## asr_q51_p      0.408
    ## asr_q52_p      0.513
    ## asr_q53_p      0.595
    ## asr_q54_p      0.565
    ## asr_q56a_p     0.428
    ## asr_q56b_p     0.362
    ## asr_q56d_p     0.248
    ## asr_q56e_p     0.257
    ## asr_q56f_p     0.375
    ## asr_q56h_p     0.427
    ## asr_q56i_p     0.422
    ## asr_q57_p      0.200
    ## asr_q58_p      0.291
    ## asr_q59_p      0.550
    ## asr_q60_p      0.520
    ## asr_q61_p      0.404
    ## asr_q63_p      0.365
    ## asr_q64_p      0.532
    ## asr_q65_p      0.418
    ## asr_q66_p      0.403
    ## asr_q67_p      0.469
    ## asr_q68_p      0.460
    ## asr_q69_p      0.465
    ## asr_q71_p      0.505
    ## asr_q72_p      0.440
    ## asr_q73_p           
    ## asr_q74_p      0.241
    ## asr_q75_p      0.365
    ## asr_q76_p      0.429
    ## asr_q77_p      0.362
    ## asr_q78_p      0.524
    ## asr_q79_p      0.207
    ## asr_q80_p           
    ## asr_q81_p      0.411
    ## asr_q82_p      0.184
    ## asr_q83_p      0.437
    ## asr_q86_p      0.544
    ## asr_q88_p     -0.158
    ## asr_q89_p      0.452
    ## asr_q90_p      0.183
    ## asr_q91_p      0.313
    ## asr_q92_p      0.260
    ## asr_q93_p      0.308
    ## asr_q94_p      0.251
    ## asr_q95_p      0.477
    ## asr_q96_p      0.309
    ## asr_q97_p      0.251
    ## asr_q98_p      0.120
    ## asr_q99_p      0.345
    ## asr_q100_p     0.461
    ## asr_q101_p     0.280
    ## asr_q102_p     0.567
    ## asr_q103_p     0.653
    ## asr_q104_p     0.325
    ## asr_q105_p     0.471
    ## asr_q106_p          
    ## asr_q107_p     0.545
    ## asr_q108_p     0.512
    ## asr_q109_p          
    ## asr_q110_p     0.224
    ## asr_q111_p     0.418
    ## asr_q112_p     0.576
    ## asr_q113_p     0.468
    ## asr_q115_p     0.570
    ## asr_q116_p     0.595
    ## asr_q118_p     0.514
    ## asr_q119_p     0.367
    ## asr_q120_p     0.267
    ## asr_q121_p     0.322
    ## asr_q122_p     0.314
    ## asr_q123_p    -0.332
    ## asr_destroyer  0.241
    ## asr_halluc     0.281
    ## asr_odd        0.476
    ## asr_sick       0.413
    ## asr_clumsy     0.404
    ## asr_anx        0.583
    ## asr_swing      0.626
    ## asr_insolvent  0.467
    ## 
    ##                   PC1
    ## SS loadings    20.092
    ## Proportion Var  0.163

``` r
subjPvalues=data.frame(pcaDf_p_CompleteSubjs,pcaDf_p_CompleteEventNames,y.pca$scores[,1])
colnames(subjPvalues)<-c('subjectkey','eventname','parentP')

### ∆∆∆ Save out first OutDF
OutDF=merge(masterdf,subjPvalues,by=c('subjectkey','eventname'))
```

    ## Warning in merge.data.frame(masterdf, subjPvalues, by = c("subjectkey", :
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.x',
    ## 'dataset_id.x', 'src_subject_id.y', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y' are duplicated in the result

``` r
dimOutDF=dim(OutDF)

# exclude subjs without data for both timepoints
OutDFBV=subset(OutDF,eventname=='baseline_year_1_arm_1')
OutDF2Y=subset(OutDF,eventname=='2_year_follow_up_y_arm_1')
# intersection of subjs in both
BothTPsubjs=intersect(OutDFBV$subjectkey,OutDF2Y$subjectkey)
# index out intersection from non tp-split df
OutDF=OutDF[OutDF$subjectkey %in% BothTPsubjs,]
outDf2dim=dim(OutDF)
print(outDf2dim)
```

    ## [1] 10076   600

``` r
dif=dimOutDF[1]-outDf2dim[1]
print(paste0(dif,' rows lost from only using subjs with both timepoints'))
```

    ## [1] "25 rows lost from only using subjs with both timepoints"

``` r
# make count version
ASRdfNum<-as.data.frame(lapply(asr[-1,11:141],as.numeric))
ASRtotal=rowSums(ASRdfNum)
# and subtract reverse score items because they were included in sum above, and modeling "happiness" as symmetric to "symptoms" seems like a strong assumption
# reverse scored = face validity AND loading in expected direction
ASRtotal=ASRtotal-ASRdfNum$asr_q02_p
ASRtotal=ASRtotal-ASRdfNum$asr_q04_p
ASRtotal=ASRtotal-ASRdfNum$asr_q15_p
ASRtotal=ASRtotal-ASRdfNum$asr_q73_p
ASRtotal=ASRtotal-ASRdfNum$asr_q80_p
ASRtotal=ASRtotal-ASRdfNum$asr_q88_p
ASRtotal=ASRtotal-ASRdfNum$asr_q106_p
ASRtotal=ASRtotal-ASRdfNum$asr_q109_p
ASRtotal=ASRtotal-ASRdfNum$asr_q123_p

# merge in (first row is colnames)
asr$parentPcount=c(NA,ASRtotal)
# fix asr age for merge
asr$interview_age=as.numeric(asr$interview_age)/12
```

    ## Warning: NAs introduced by coercion

``` r
# set subjectkey to factor for merge
asr$subjectkey<-as.factor(asr$subjectkey)

# merge
OutDF=merge(OutDF,asr,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(OutDF, asr, by = c("subjectkey", "eventname", :
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.x',
    ## 'dataset_id.x', 'src_subject_id.y', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'src_subject_id.y' are
    ## duplicated in the result

``` r
### confirm tangentivity of sidequest
print(cor.test(OutDF$parentP,OutDF$parentPcount))
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  OutDF$parentP and OutDF$parentPcount
    ## t = 924.7, df = 10074, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.9939289 0.9943838
    ## sample estimates:
    ##       cor 
    ## 0.9941608

``` r
plotdf<-OutDF[,c('parentP','parentPcount')]
ComparePlot<-ggplot(data=plotdf,aes(x=parentPcount,y=parentP))+geom_smooth(method='lm',col='black')+geom_point(alpha=.05)+theme_classic()
print(ComparePlot)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](FunctionalEquivalence_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
