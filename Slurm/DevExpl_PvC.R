# read in master df and calculate deviance in g explained by cbcl scales vs asr scales
library(mgcv)

# load in data
masterdf=readRDS('/oak/stanford/groups/leanew1/users/apines/data/gp/gp_masterdf.rds')
# FYI
print('dimensions of dataframe')
dim(masterdf)

# garner num subjs for bootstrapping
subjs=unique(masterdf$subjectkey)
numSubjs=length(subjs)
# cut df to just variables of interest to speed stuff up # add cbcl subscales
masterdf=masterdf[,c('parentPcount','cbcl_scr_syn_totprob_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_somatic_r','cbcl_scr_syn_anxdep_r','cbcl_scr_syn_thought_r','cbcl_scr_syn_withdep_r','cbcl_scr_syn_social_r','cbcl_scr_syn_attention_r','cbcl_scr_syn_rulebreak_r','cbcl_scr_syn_aggressive_r','ASRAnxDepr','ASRWithdrawn','ASRSomatic','ASRThought','ASRAttn','ASRAggr','ASRIntrusive','ASRRulB','ASRInt','ASRExt','g','subjectkey','interview_age')]
# get length of df for later
lenDF=dim(masterdf)[1]
# convert cbcl scores to numeric
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
# convert parentPcount to numeric
masterdf$parentPcount=as.numeric(masterdf$parentPcount)
masterdf$cbcl_scr_syn_totprob_r=as.numeric(masterdf$cbcl_scr_syn_totprob_r)
masterdf$cbcl_scr_syn_internal_r=as.numeric(masterdf$cbcl_scr_syn_internal_r)
masterdf$cbcl_scr_syn_external_r=as.numeric(masterdf$cbcl_scr_syn_external_r)
masterdf$cbcl_scr_syn_somatic_r=as.numeric(masterdf$cbcl_scr_syn_somatic_r)
masterdf$cbcl_scr_syn_anxdep_r=as.numeric(masterdf$cbcl_scr_syn_anxdep_r)
masterdf$cbcl_scr_syn_thought_r=as.numeric(masterdf$cbcl_scr_syn_thought_r)
masterdf$cbcl_scr_syn_withdep_r=as.numeric(masterdf$cbcl_scr_syn_withdep_r)
masterdf$cbcl_scr_syn_social_r=as.numeric(masterdf$cbcl_scr_syn_social_r)
masterdf$cbcl_scr_syn_attention_r=as.numeric(masterdf$cbcl_scr_syn_attention_r)
masterdf$cbcl_scr_syn_rulebreak_r=as.numeric(masterdf$cbcl_scr_syn_rulebreak_r)
masterdf$cbcl_scr_syn_aggressive_r=as.numeric(masterdf$cbcl_scr_syn_aggressive_r)
masterdf$ASR_anxdep=as.numeric(masterdf$ASRAnxDepr)
masterdf$ASR_withdep=as.numeric(masterdf$ASRWithdrawn)
masterdf$ASR_somatic=as.numeric(masterdf$ASRSomatic)
masterdf$ASR_thought=as.numeric(masterdf$ASRThought)
masterdf$ASR_attention=as.numeric(masterdf$ASRAttn)
masterdf$ASR_aggressive=as.numeric(masterdf$ASRAggr)
masterdf$ASR_rulebreak=as.numeric(masterdf$ASRRulB)
masterdf$ASR_intrusive=as.numeric(masterdf$ASRIntrusive)
masterdf$ASRInt=as.numeric(masterdf$ASRInt)
masterdf$ASRExt=as.numeric(masterdf$ASRExt)
# initialize a vector for each of 10k bootstraps that will hold deviance explained by a single variable
devExplTotProb=rep(0,10000)
devExplInternal=rep(0,10000)
devExplExternal=rep(0,10000)
devExplSomatic=rep(0,10000)
devExplAnxDep=rep(0,10000)
devExplThought=rep(0,10000)
devExplWithDep=rep(0,10000)
devExplSocial=rep(0,10000)
devExplAttention=rep(0,10000)
devExplRuleBreak=rep(0,10000)
devExplAggressive=rep(0,10000)
devExplParentPcount=rep(0,10000)
devExplParentInternal=rep(0,10000)
devExplParentExternal=rep(0,10000)
devExplParentSomatic=rep(0,10000)
devExplParentAnx=rep(0,10000)
devExplParentThought=rep(0,10000)
devExplParentWith=rep(0,10000)
devExplParentIntr=rep(0,10000)
devExplParentAttn=rep(0,10000)
devExplParentRule=rep(0,10000)
devExplParentAgg=rep(0,10000)
# set seed
set.seed(1)
for (b in 1:2000) {
  print(b)
  # get subjects to include in this bootstrap
  BootSubjs <- sample(subjs, numSubjs, replace = TRUE)
  
  # Create an empty dataframe to store the resampled observations
  bootSamp <- data.frame()
  
  tryCatch({
    for (j in 1:length(BootSubjs)) {
      subject_obs <- masterdf[masterdf$subjectkey == BootSubjs[j], ]
      bootSamp <- rbind(bootSamp, subject_obs)
    }
    
    # fit each model explaining g with a single scale
    TotProbMod <- bam(g ~ s(cbcl_scr_syn_totprob_r), data = bootSamp)
    InternalMod <- bam(g ~ s(cbcl_scr_syn_internal_r), data = bootSamp)
    ExternalMod <- bam(g ~ s(cbcl_scr_syn_external_r), data = bootSamp)
    SomaticMod <- bam(g ~ s(cbcl_scr_syn_somatic_r), data = bootSamp)
    AnxDepMod <- bam(g ~ s(cbcl_scr_syn_anxdep_r), data = bootSamp)
    ThoughtMod <- bam(g ~ s(cbcl_scr_syn_thought_r), data = bootSamp)
    WithDepMod <- bam(g ~ s(cbcl_scr_syn_withdep_r), data = bootSamp)
    SocialMod <- bam(g ~ s(cbcl_scr_syn_social_r), data = bootSamp)
    AttentionMod <- bam(g ~ s(cbcl_scr_syn_attention_r), data = bootSamp)
    RuleBreakMod <- bam(g ~ s(cbcl_scr_syn_rulebreak_r), data = bootSamp)
    AggressiveMod <- bam(g ~ s(cbcl_scr_syn_aggressive_r), data = bootSamp)
    ParentPcountMod <- bam(g ~ s(parentPcount), data = bootSamp)
    ParentInternalMod <- bam(g ~ s(ASRInt), data = bootSamp)
    ParentExternalMod <- bam(g ~ s(ASRExt), data = bootSamp)
    ParentSomaticMod <- bam(g ~ s(ASR_somatic), data = bootSamp)
    ParentAnxMod <- bam(g ~ s(ASR_anxdep), data = bootSamp)
    ParentThoughtMod <- bam(g ~ s(ASR_thought), data = bootSamp)
    ParentWithMod <- bam(g ~ s(ASR_withdep), data = bootSamp)
    ParentIntrMod <- bam(g ~ s(ASR_intrusive), data = bootSamp)
    ParentAttnMod <- bam(g ~ s(ASR_attention), data = bootSamp)
    ParentRuleMod <- bam(g ~ s(ASR_rulebreak), data = bootSamp)
    ParentAggMod <- bam(g ~ s(ASR_aggressive), data = bootSamp)
    
    # get deviance explained by each model
    devExplTotProb[b] <- summary(TotProbMod)$dev.expl
    devExplInternal[b] <- summary(InternalMod)$dev.expl
    devExplExternal[b] <- summary(ExternalMod)$dev.expl
    devExplSomatic[b] <- summary(SomaticMod)$dev.expl
    devExplAnxDep[b] <- summary(AnxDepMod)$dev.expl
    devExplThought[b] <- summary(ThoughtMod)$dev.expl
    devExplWithDep[b] <- summary(WithDepMod)$dev.expl
    devExplSocial[b] <- summary(SocialMod)$dev.expl
    devExplAttention[b] <- summary(AttentionMod)$dev.expl
    devExplRuleBreak[b] <- summary(RuleBreakMod)$dev.expl
    devExplAggressive[b] <- summary(AggressiveMod)$dev.expl
    devExplParentPcount[b] <- summary(ParentPcountMod)$dev.expl
    devExplParentInternal[b] <- summary(ParentInternalMod)$dev.expl
    devExplParentExternal[b] <- summary(ParentExternalMod)$dev.expl
    devExplParentSomatic[b] <- summary(ParentSomaticMod)$dev.expl
    devExplParentAnx[b] <- summary(ParentAnxMod)$dev.expl
    devExplParentThought[b] <- summary(ParentThoughtMod)$dev.expl
    devExplParentWith[b] <- summary(ParentWithMod)$dev.expl
    devExplParentIntr[b] <- summary(ParentIntrMod)$dev.expl
    devExplParentAttn[b] <- summary(ParentAttnMod)$dev.expl
    devExplParentRule[b] <- summary(ParentRuleMod)$dev.expl
    devExplParentAgg[b] <- summary(ParentAggMod)$dev.expl
  }, error = function(e) {
    # Assign "NA" to positions in case of error
    devExplTotProb[b] <- NA
    devExplInternal[b] <- NA
    devExplExternal[b] <- NA
    devExplSomatic[b] <- NA
    devExplAnxDep[b] <- NA
    devExplThought[b] <- NA
    devExplWithDep[b] <- NA
    devExplSocial[b] <- NA
    devExplAttention[b] <- NA
    devExplRuleBreak[b] <- NA
    devExplAggressive[b] <- NA
    devExplParentPcount[b] <- NA
    devExplParentInternal[b] <- NA
    devExplParentExternal[b] <- NA
    devExplParentSomatic[b] <- NA
    devExplParentAnx[b] <- NA
    devExplParentThought[b] <- NA
    devExplParentWith[b] <- NA
    devExplParentIntr[b] <- NA
    devExplParentAttn[b] <- NA
    devExplParentRule[b] <- NA
    devExplParentAgg[b] <- NA
  })
}
# save results
outdf=data.frame(devExplTotProb,devExplInternal,devExplExternal,devExplSomatic,devExplAnxDep,devExplThought,devExplWithDep,devExplSocial,devExplAttention,devExplRuleBreak,devExplAggressive,devExplParentPcount,devExplParentInternal,devExplParentExternal,devExplParentSomatic,devExplParentAnx,devExplParentThought,devExplParentWith,devExplParentIntr,devExplParentAttn,devExplParentRule,devExplParentAgg)
saveRDS(outdf,'/oak/stanford/groups/leanew1/users/apines/data/gp/PvC_gdevExplBoots.rds')
