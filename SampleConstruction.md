SampleConstruction
================
Adam
2023-01-25

``` r
# libraries
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
# load in data
cbcl=read.delim('~/Downloads/Package_1205735/abcd_cbcl01.txt')
cbcls=read.delim('/Users/panlab/Downloads/Package_1205735/abcd_cbcls01.txt')
# subset timepoints
cbclsBV=subset(cbcls,eventname=='baseline_year_1_arm_1')
cbcls2=subset(cbcls,eventname=='2_year_follow_up_y_arm_1')
#### add non summary items to cbcl for schoolwork q
# subset timepoints
cbclBV=subset(cbcl,eventname=='baseline_year_1_arm_1')
cbcl2=subset(cbcl,eventname=='2_year_follow_up_y_arm_1')
# merge with other cbcl
cbclsBV=merge(cbclsBV,cbclBV,by=c('subjectkey','eventname'))
cbcls2=merge(cbcls2,cbcl2,by=c('subjectkey','eventname'))

# load in ASR
asr=read.delim('~/Downloads/Package_1207917/pasr01.txt',na.strings=c("","NA"))

# load in clinicalish data
clinc1=read.delim('/Users/panlab/Downloads/Package_1205908/abcd_ksad501.txt')
kBV=subset(clinc1,eventname=='baseline_year_1_arm_1')
k1=subset(clinc1,eventname=='1_year_follow_up_y_arm_1')
k2=subset(clinc1,eventname=='2_year_follow_up_y_arm_1')
clinc1_p=read.delim('/Users/panlab/Downloads/Package_1205908/abcd_ksad01.txt')

# read in grades, annoying distinction between tp1 measure (decent granulrity) and tp2 measure (high granularity)
gradesInfoBV=readRDS('~/Downloads/DEAP-data-download-13.rds')
gradesInfoBV=subset(gradesInfoBV,event_name=='baseline_year_1_arm_1')
gradesInfoBV$Grades<-as.numeric(gradesInfoBV$ksads_back_grades_in_school_p)
gradesInfoBV$Grades[gradesInfoBV$Grades==-1]=NA
gradesInfoBV$eventname=gradesInfoBV$event_name
gradesInfoBV$subjectkey=gradesInfoBV$src_subject_id
# for this one, the key is 1 = A's, 2 = B's, 3 = C's, 4 = D's, 5 = F's, -1 = NA
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

################################### ∆∆∆∆∆

# initialize master df
masterdf<-merge(cbcls,cbcl,by=c('subjectkey','eventname','interview_age','src_subject_id'))
masterdf<-merge(masterdf,clinc1,by=c('subjectkey','eventname','interview_age','src_subject_id'))
masterdf<-merge(masterdf,clinc1_p,by=c('subjectkey','eventname','interview_age','src_subject_id'))
```

    ## Warning in merge.data.frame(masterdf, clinc1_p, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated
    ## in the result

``` r
masterdf<-merge(masterdf,gradesInfo,by=c('subjectkey','eventname'))
```

    ## Warning in merge.data.frame(masterdf, gradesInfo, by = c("subjectkey",
    ## "eventname")): column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated
    ## in the result

``` r
# load in a DEAP file for rel_family_ID
DEAP=readRDS('~/Downloads/DEAP-data-download-13.rds')
DEAP$subjectkey<-DEAP$src_subject_id
DEAP$eventname=DEAP$event_name
DEAP=DEAP[,c('rel_family_id','subjectkey','eventname')]

#### clean data

# convert to numeric
masterdf$interview_age<-as.numeric(masterdf$interview_age)
# cbcl sum indep. of the q of interest
masterdf$cbcl_scr_syn_totprob_r<-as.numeric(masterdf$cbcl_scr_syn_totprob_r)-as.numeric(masterdf$cbcl_q61_p)
masterdf$cbcl_q61_p<-as.ordered(masterdf$cbcl_q61_p)
masterdf$subjectkey<-as.factor(masterdf$subjectkey)
# clean data
masterdf$cbcl_scr_syn_totprob_r<-as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r<-as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r<-as.numeric(masterdf$cbcl_scr_syn_external_r)
masterdf$subjectkey<-as.factor(masterdf$subjectkey)
masterdf$cbcl_q61_p<-as.ordered(masterdf$cbcl_q61_p)
# remove instances of NA tot probs
masterdf=masterdf[!is.na(masterdf$cbcl_scr_syn_totprob_r),]
# and for q61
masterdf=masterdf[!is.na(masterdf$cbcl_q61_p),]
# and for is empty
masterdf=masterdf[!is.empty(masterdf$cbcl_scr_syn_totprob_r),]
masterdf=masterdf[!is.empty(masterdf$cbcl_q61_p),]
# r can't take the hint
masterdf=masterdf[!masterdf$cbcl_scr_syn_totprob_r=='',]
masterdf=masterdf[!masterdf$cbcl_q61_p=='',]

#### now load in cognitive data
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
masterdf<-merge(masterdf,othCog,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(masterdf, othCog, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y',
    ## 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y' are duplicated in the result

``` r
masterdf<-merge(masterdf,littleMan,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(masterdf, littleMan, by = c("subjectkey",
    ## "eventname", : column names 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y',
    ## 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'src_subject_id.x',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'src_subject_id.y' are
    ## duplicated in the result

``` r
masterdf<-merge(masterdf,DEAP,by=c('subjectkey','eventname'))
```

    ## Warning in merge.data.frame(masterdf, DEAP, by = c("subjectkey", "eventname")):
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.y', 'dataset_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y', 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'src_subject_id.x', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'src_subject_id.y' are duplicated in the result

``` r
# clean age
masterdf$interview_age<-as.numeric(masterdf$interview_age)/12
```

``` r
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
```

``` r
#### calculate PCs

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

# for isolating PCA df
pcVars=c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","pea_ravlt_ld","lmt_scr_perc_correct")

# get other vars of interest to check for complete cases
KidVarsOfInt=c('Grades','cbcl_scr_syn_totprob_r','cbcl_scr_syn_external_r','cbcl_scr_syn_internal_r')


# columns of interest to gauge completeness of
ColsOfInt=asr[,c(11:141)]
ASRVarsOfInt=colnames(ColsOfInt)
# clean age in ASR
asr$interview_age<-as.numeric(asr$interview_age)/12
```

    ## Warning: NAs introduced by coercion

``` r
# merge with master
masterdf=merge(asr,masterdf,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(asr, masterdf, by = c("subjectkey", "eventname", :
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'collection_id.x', 'dataset_id.x', 'interview_date.x',
    ## 'sex.x', 'collection_title.x', 'collection_id.y', 'dataset_id.y',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.x',
    ## 'dataset_id.x', 'interview_date.x', 'sex.x', 'collection_title.x',
    ## 'collection_id.y', 'dataset_id.y', 'src_subject_id.x', 'interview_date.y',
    ## 'sex.y', 'collection_title.y', 'collection_id.y', 'dataset_id.y',
    ## 'src_subject_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are
    ## duplicated in the result

``` r
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

print(dim(masterdf))
```

    ## [1] 11448  2516

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

# make family ID for those with families represented in ABCD
masterdf$rel_family_id=masterdf$id_fam

print(dim(masterdf))
```

    ## [1] 10100  2518

``` r
# pea_wiscv_tss, nihtbx_list_uncorrected, and nihtbx_cardsort_uncorrected taken out for lack of longitudinal coverage
pcaDf<-masterdf[,pcVars]
```

``` r
# derive Cognitive PCS

# derive pcs
Y = as.matrix(scale(pcaDf[complete.cases(pcaDf[,pcVars]),pcVars]))
# equiv for binding scores to IDs and eventnames
pcVarsAndIDs=c("nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_pattern_uncorrected","nihtbx_picture_uncorrected","nihtbx_reading_uncorrected","pea_ravlt_ld","lmt_scr_perc_correct","subjectkey","eventname")
Yextended=masterdf[complete.cases(masterdf[,pcVarsAndIDs]),pcVarsAndIDs]
ncomp = 3
y.pca = psych::principal(Y, rotate="varimax", nfactors=ncomp, scores=TRUE)
y.pca$loadings
```

    ## 
    ## Loadings:
    ##                             RC1    RC2    RC3   
    ## nihtbx_picvocab_uncorrected  0.752  0.189       
    ## nihtbx_flanker_uncorrected   0.203  0.820       
    ## nihtbx_pattern_uncorrected   0.163  0.843       
    ## nihtbx_picture_uncorrected   0.604  0.250       
    ## nihtbx_reading_uncorrected   0.709  0.202  0.175
    ## pea_ravlt_ld                 0.763              
    ## lmt_scr_perc_correct                       0.983
    ## 
    ##                  RC1   RC2   RC3
    ## SS loadings    2.088 1.527 1.013
    ## Proportion Var 0.298 0.218 0.145
    ## Cumulative Var 0.298 0.516 0.661

``` r
# assign scores to subjs
Yextended$g<-y.pca$scores[,1]
# merge in cog data

masterdf$g<-Yextended$g
```

``` r
# save out all timepoints df for bootstrapping both-tp-fits
saveRDS(masterdf,'~/DfWithGrades.rds')
print(dim(masterdf))
```

    ## [1] 10100  2519

``` r
#### ∆∆∆ Note this is without adult P
```

``` r
# derive adult p

# columns of interest to gauge completeness of
ColsOfInt=asr[,c(11:141)]
# retain complete cases
completeInd=ColsOfInt[rowSums(is.na(ColsOfInt)) == 0,]
# merge with master
masterdf=merge(asr,masterdf,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(asr, masterdf, by = c("subjectkey", "eventname", :
    ## column names 'src_subject_id.x', 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.x',
    ## 'dataset_id.x', 'interview_date.x', 'sex.x', 'collection_title.x',
    ## 'collection_id.y', 'dataset_id.y', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'collection_id.x', 'dataset_id.x', 'src_subject_id.y',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'src_subject_id.x', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'collection_id.y', 'dataset_id.y', 'src_subject_id.y',
    ## 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated in the result

``` r
ASRcolnames=asr[1,]

#######MICHELI APPROACH TO FACTOR DELINEATION OF ASR
# composite creation

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
# polycormat
asrBV_qs_cormat=hetcor(asrBV_qs)
```

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q02_p and asr_q32_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q04_p and asr_q12_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q04_p and asr_q32_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q04_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q14_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q14_p and asr_q78_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q15_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q25_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q26_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q27_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q28_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q29_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q30_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q32_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q35_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q36_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q52_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56c_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56h_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q60_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q64_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q67_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q74_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q84_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q85_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q89_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q96_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q101_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q23_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q25_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q28_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q30_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q35_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q41_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q48_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56h_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q25_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q28_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q30_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q35_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q41_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q48_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56h_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q25_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q25_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q25_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q28_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q30_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q30_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q30_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q30_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q37_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q35_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q35_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q37_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q41_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q41_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q41_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q42_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q42_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q48_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q45_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q45_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q46_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q46_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q47_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q47_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q48_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q48_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q48_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q56f_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q55_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q55_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q61_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q65_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56h_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56h_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q58_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q62_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q62_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q62_p and asr_q106_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q62_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q65_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q65_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q65_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q68_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q69_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q69_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q73_p and asr_q94_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q77_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q78_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q81_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q83_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q86_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q93_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q94_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q94_p and asr_q106_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q94_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q94_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q95_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q98_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q99_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q100_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q100_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q102_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q102_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q103_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q104_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q104_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q105_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q105_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q107_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q108_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q108_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q112_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q118_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q121_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in hetcor.data.frame(asrBV_qs): the correlation matrix has been adjusted
    ## to make it positive-definite

``` r
asr2_qs_cormat=hetcor(asr2_qs)
```

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q06_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q01_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q03_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q04_p and asr_q63_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q04_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q05_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q05_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q05_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q05_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q08_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q09_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q32_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q106_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q06_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q08_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q09_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q10_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q11_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q12_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q13_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q14_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q14_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q14_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q15_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q21_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q17_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q18_p and asr_q49_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q22_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q24_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q25_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q26_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q27_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q28_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q29_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q30_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q32_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q33_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q35_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q36_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q51_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q52_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56b_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56c_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56h_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q60_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q62_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q64_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q67_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q84_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q85_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q89_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q90_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q96_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q101_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q21_p and asr_q122_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q26_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q33_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q51_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56c_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q60_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q74_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q22_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q23_p and asr_q98_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q26_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q31_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q33_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q38_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q43_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q51_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56c_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q60_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q74_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q24_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q26_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q27_p and asr_q88_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q28_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q29_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q30_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q31_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q39_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q32_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q33_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q35_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q36_p and asr_q40_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q38_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q52_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q64_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q39_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q42_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q45_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q46_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q47_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q50_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q52_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q53_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q54_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q55_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56c_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56e_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56f_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56h_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q56i_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q58_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q59_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q63_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q64_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q67_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q68_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q69_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q71_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q75_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q78_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q80_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q81_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q83_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q84_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q85_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q87_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q115_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q40_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q42_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q42_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q42_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q44_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q49_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q56a_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q56b_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q43_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q60_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q44_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q45_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q45_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q45_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q46_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q46_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q46_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q47_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q47_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q47_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q56d_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q66_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q49_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q50_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q51_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q51_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q53_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q53_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q54_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56a_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56b_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56c_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56c_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56d_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56d_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56d_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56d_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q73_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56e_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q80_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56f_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q73_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q56i_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q58_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q58_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q59_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q60_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q60_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q60_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q61_p and asr_q98_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q73_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q63_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q72_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q66_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q69_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q69_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q69_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q71_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q74_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q76_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q77_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q89_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q72_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q73_p and asr_q74_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q74_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q74_p and asr_q88_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q74_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q74_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q74_p and asr_q123_p produced a warning:
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q106_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q76_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q77_p and asr_q80_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q77_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q77_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q77_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q78_p and asr_q79_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q78_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q78_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q79_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q86_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q89_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q95_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q80_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q81_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q81_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q83_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q83_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q86_p and asr_q92_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q86_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q89_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q88_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q89_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q89_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q93_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q98_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q99_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q100_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q102_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q103_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q104_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q105_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q107_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q108_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q92_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q93_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q95_p and asr_q109_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q98_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q98_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q98_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q99_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q100_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q102_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q103_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q104_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q105_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q106_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q107_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q108_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q110_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q114_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q109_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q111_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q112_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q113_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q116_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q117_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q118_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q119_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q120_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in FUN(X[[i]], ...): polychoric correlation between variables asr_q110_p and asr_q121_p produced warnings:
    ##    NaNs produced
    ##    NaNs produced
    ##    NaNs produced

    ## Warning in hetcor.data.frame(asr2_qs): the correlation matrix has been adjusted
    ## to make it positive-definite

``` r
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

OutDF=merge(masterdf,subjPvalues,by=c('subjectkey','eventname'))
```

    ## Warning in merge.data.frame(masterdf, subjPvalues, by = c("subjectkey", : column
    ## names 'src_subject_id.x', 'collection_id.x', 'dataset_id.x', 'interview_date.x',
    ## 'sex.x', 'collection_title.x', 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y',
    ## 'collection_id.x', 'dataset_id.x', 'src_subject_id.y', 'interview_date.x',
    ## 'sex.x', 'collection_title.x', 'collection_id.y', 'dataset_id.y',
    ## 'src_subject_id.x', 'interview_date.y', 'sex.y', 'collection_title.y',
    ## 'collection_id.y', 'dataset_id.y', 'src_subject_id.y', 'interview_date.y',
    ## 'sex.y', 'collection_title.y' are duplicated in the result

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
OutDF=merge(OutDF,asr,by=c('subjectkey','eventname','interview_age'))
```

    ## Warning in merge.data.frame(OutDF, asr, by = c("subjectkey", "eventname", :
    ## column names 'collection_id.x', 'dataset_id.x', 'interview_date.x', 'sex.x',
    ## 'collection_title.x', 'src_subject_id.x', 'collection_id.x', 'dataset_id.x',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.x',
    ## 'dataset_id.x', 'interview_date.x', 'sex.x', 'collection_title.x',
    ## 'collection_id.y', 'dataset_id.y', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'collection_id.x', 'dataset_id.x', 'src_subject_id.y',
    ## 'interview_date.x', 'sex.x', 'collection_title.x', 'collection_id.y',
    ## 'dataset_id.y', 'src_subject_id.x', 'interview_date.y', 'sex.y',
    ## 'collection_title.y', 'collection_id.y', 'dataset_id.y', 'src_subject_id.y',
    ## 'interview_date.y', 'sex.y', 'collection_title.y', 'collection_id.y',
    ## 'dataset_id.y', 'interview_date.y', 'sex.y', 'collection_title.y' are duplicated
    ## in the result

``` r
print(dim(OutDF))
```

    ## [1] 10100  2811

``` r
saveRDS(OutDF,'~/OutDfFull.rds')

# convert to one row per subj for temporal precedence analyses
OutDFBV=subset(OutDF,eventname=='baseline_year_1_arm_1')
OutDF2Y=subset(OutDF,eventname=='2_year_follow_up_y_arm_1')
OutDFTmpPrec<-merge(OutDFBV,OutDF2Y,by='subjectkey')
print(dim(OutDFTmpPrec))
```

    ## [1] 5038 5621

``` r
saveRDS(OutDFTmpPrec,'~/OutDFTmpPrec.rds')
```
