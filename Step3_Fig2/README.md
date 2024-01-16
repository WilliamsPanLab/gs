# Figure 2: Parental and demographic analyses

Here, the goal is to submit 10k iteration of non-linear firts of g to p on sherlock for parental mental health. That script [is here](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/Fig_2_parents.R).
As prior, we'll take masterdf from the sample construction step and run 10,000 [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) to gauge model consistency over different permutations of participants. Because we have repeated measures (2x), we need to take care to remove entire participants for each bootstrap, not just observations.

## Bootstrapping/computational steps using masterdf from sample construction

1. Port masterdf over to Sherlock, or whatever your compute cluster name is, with scp. 
2. module load R/4.1 on cluster and open R (terminal). Subsequent steps in Fig_2_parents.R, which is called to be "sbatched" (slurm equivalent to qsub on SGE) by [this script](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/sbatch_Fig2_Parents.sh)
4. Load needed master dataframe
5. Glean number of subjects for bootstrapping purposes later
6. Convert all Child Behavioral Checklist and Adult Self Report scores to numeric
7. Initialize output vectors: initializing the values will make it slightly faster. Initialize them as 0s for
   
  A) Maximum values of symptom scores for each iteration
  B) predictied derivatives, note we get an entire range of values for each iteration.
  C) predicted fits children below the poverty line (P) and above (R).
  D) predicted fits for poverty vs. nonpoverty vs. nonpoverty size-matched (same n as poverty subsample to preclude possibility that n differences drive group differences)
  E) difference-in-AIC between full and reduced models across bootstraps. Diff Pseudo is for permuted poverty labels, so we can see the true difference in AIC vs. a "null" difference in AIC for manual significance testing.

9. Now we can actually begin bootstrap. We will run this 10k times as for figure 1.
10. Perform the bootstrap for each iteration. Recall that we need entire participants in/out for each iteration. We shouldn't be pulling one observation from one participant out and retaining the other inside of a bootstrap. We'll make a new dataframe by first randomly sampling participants, and then looping through each randomly sampled PTID and populating the bootstrap dataframe with that participant's data iteratively
11. Assign same number of participants actually below poverty line to a "pseudopoverty" group to test null permutations vs. observed.
12. Now we can derive our maximum values for each iteration
13. We'll also get counts of how many participants are in F group for sex (female) and poverty group (income below 25k per year). As prior, sex is labeled as seg because github co-pilot doesn't autocomplete certain words.
14. Randomly assign n people to a pseudopoverty group for null models: should be same n as in actual poverty group. For psueodopoverty2, this group is explicitly selected from the non-poverty group.
15. (I) Make some GAMs that do and don't account for poverty: we will get our predicted values by poverty from these models later.
16. (II) Predict variables of interest with the splines fit to said models. First we'll make a dataframe with each symptom level sampled with seq(0:bpmax), the median age of the bootstrap, and then we'll cycle out sex and poverty for different predictions.
17. Set column names for prediction data frames.
18. Get fits for kids, below poverty line, and above poverty line. Get derivatives of these fits for saving out.
19. Populate derivatives and fit vectors with values for this bootstrap.
21. Fit models on each group.
22. Recover derivatives of g~p (cbcl_scr_syn_totprob_r) across spectrum of p for each model, print out to vectors.
23. Save data (after looping).

## Upon completion, we should have gpFitBoots_cbclasr.rds, gpFitBoots_asr_pNp.RDS (poverty/nonpoverty), gpFitBoots_cbcl_pNp.rds, gpDiffBoots_cbcl.rds, gpDiffBoots_asr.rds, gpDiffBoots_cbclPseudo.Rds, and gpDiffBoots_asrPseudo.Rds  for use in the subsequent plotting script

## Plotting steps (assuming above .rds files are derived and locally available). Script referenced is available [here](https://github.com/WilliamsPanLab/gp/blob/master/Step3_Fig2/Fig2.md)

1. Prepare yourself, this is a long .md. Some of the code is redundant rather than clever functions for easier editing (definitely not because of a cleverness deficit of the author)
2. Load needed libraries.
3. Chunk 2 establishes some functions that we'll use throughout the markdown. One is to plot bootstraps with a fair amount of ggplot specifications, and the other is to find the furthest extent of symptoms reported in each bootstrap. For context: the significance testing we've employed uses the total number of bootstraps to derive p-values. If certain levels of symptoms are not included in some bootstraps, this changes how we can consider statistical significance. To maximize statistical certainty, only symptom ranges that were included in all bootstraps are considered. Therefore, we need this function `find_furthest_nonzero` to determine what the highest symptom count that was included across all all 10k iterations was. Finally, we set a color palette that will be used throughout.
4. Load in masterdf, glean further evidence of nonlinearity + justification for use of splines with linear vs. nonlinear AIC, as in Step 2 (macro step 2, not step 2 in this .md file).
5. Calculate clinical/subclinical boundaries. as prior. Parental mental health subscales come after CBCL subscales.
6. The next chunk loads in the bootstrapped fits poverty vs. nonpoverty fits for g~symptoms. This is printed out for every subscale as three separate chunks: poverty and nonpoverty fits co-plotted, then poverty-specific derivatives and then non-poverty specific derivatives. Note 95% confidence is obtained using pos/neg counts > 9,5000, and 99% confidence with 9,900 subbed in. Both are depicted in the manuscript.
7. True AIC difference from including poverty interaction evaluated versus null models for manual p-values.
8. Next, parental symptom bootstrap fits and derivatives are evaluated across bootstraps, as prior. Note splitting of first 1:2000 values as .rds 1, 2001:4000 as .rds 2, as prior, due to multihreading. If you run your [bootstraps](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/Fig_2_parents.R) in one shot this is not neccessary.
9. Deviance explained in g as a product of child and parental symptoms across bootstraps is loaded in and plotted next.
10. AIC of adding parental symptom scores is evaluated next.
11. Simple plot of child p-factor in poverty and non-poverty conditions is plotted next.
12. Some basic stats (median problems below and above poverty line, t-test of the difference) are derived next, independently in both timepoints of the dataframe.


