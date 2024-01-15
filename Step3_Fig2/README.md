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

1. Prepare yourself, this is a long .md. Some of the code is redundant rather than created as functions for easier editing.
2. Load needed libraries.
3. Chunk 2 establishes some functions that we'll use throughout the markdown. One is to plot bootstraps with a fair amount of ggplot specifications, and the other is to find the furthest extent of symptoms reported in each bootstrap. For context: the significance testing we've employed uses the total number of bootstraps to derive p-values. If certain levels of symptoms are not included in some bootstraps, this changes how we can consider statistical significance. To maximize statistical certainty, only symptom ranges that were included in all bootstraps are considered. Therefore, we need this function `find_furthest_nonzero` to determine what the highest symptom count that was included across all all 10k iterations was. Finally, we set a color palette that will be used throughout.
4. Load in masterdf, glean further evidence of nonlinearity + justification for use of splines with linear vs. nonlinear AIC, as in Step 2 (macro step 2, not step 2 in this .md file).
5. Calculate clinical/subclinical boundaries. as prior.
6. Plot child p-factor vs. parental p-factor, print correlation
7. The next chunk loads in the bootstrapped fits for adult self-report (`F3_gpFits.rds`) and determines borderline and full clinical cutoffs based off of the t-scores of data from masterdf (obtained in the sample construction step).
8. The next chunk gets the median value of girl fits at each symptom count out to where there is bootstrap coverage. There is only boostrap coverage to the # of symptoms that was consistently represented in every bootstrap.
9. Isolate boy fits, obtain median, and truncate at the maximum value seen across girl fits for equivalence. Boys have slightly more symptoms, so we are using the lowest common denominator (max girl symptom #) here.
10. Combine both into a dataframe and plot 'em.
11. Extract derivatives of these fits. Here we calculate derivatives directly from the fits for girls and boys, plop 'em in derivative_matrix.
12. Now we'll do our manual significance testing at 99% confidence. Specifically, we'll find where derivatives (slopes) are positive or negative across all 10,000 iterations. This is a way of testing if the 99% confidence interval for the true slope excludes 0. We'll then plot significant locations (spans of symptoms where slope is significant) using geom_raster. This is equivalent to the approach used in [this paper](https://www.sciencedirect.com/science/article/pii/S1878929320300360), but manual bootstrapping was needed in this situation (hence sherlock/slurm). We'll plot out the significant slopes for all subscales.
13. Get data for children above below the poverty line.
14. Create grades plots: box and whiskers, and directly from masterdf rather than sherlock derivatives.
15. Create deviance explained plots. These are fairly straightforward as well: we're just plotting the deviance explained in each bootstrap as a function of whether that model sought to explain deviance in pscyhopathology with g, grades, parental P, or some combination of those three. Cross-sectional is first, which is same-timepoint deviance explained, and longitudinal is second, which is deviance explained in timepoint 2 above and beyond the measurement of the predicted variable at timepoint 1.
16. A fairly basic g-Grades correspondence is established as a sanity check for a supplemental figure.
17. Load in the poverty vs. nonpoverty vs. n-matched nonpoverty fit derivatives. We will compare them in the symptom range consistently sampled in the poverty analysis, as the poverty group has the least coverage. Note we are comparing the derivatives between these three groups, and note we are comparing at 90% confidence rather than 99.
18. The last step is a temporal precedence analyses. Here, we use lavaan to re-establish what was shown in [this paper](https://pubmed.ncbi.nlm.nih.gov/34332330/). The general finding is that symptoms tend to precede cognitive impairment, except for internalizing. More information on lavaan can [be found here](https://cran.r-project.org/web/packages/lavaan/lavaan.pdf).

