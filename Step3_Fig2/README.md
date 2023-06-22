# Figure 1: Non-linear fits to asr (adult self report, parental) subscales

Here, the goal is to submit 10k iterations of non-linear fits of g to parental mental health subscales subscales on sherlock (our computing cluser, slurm-based). So we'll take masterdf from the sample construction step and run 10,000 [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) to gauge model consistency over different permutations of participants. Because we have repeated measures (2x), we need to take care to remove entire participants for each bootstrap, not just observations.

These analyses take place over three scripts: sample construction (done already), [Fig2_Parents.R](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/Fig_2_parents.R), and then plotting the values derived from the bootstrapping ([Fig1.md/.rmd](https://github.com/WilliamsPanLab/gp/blob/master/Figures/code/Fig2.md))

## Bootstrapping/computational steps (assuming sample construction is done)

1. Port masterdf over to Sherlock with scp. Subsequent steps in Fig_2_parents.R
2. module load R/4.1 on sherlock and open R (terminal)
3. Load needed master dataframe
4. Glean number of subjects for bootstrapping purposes later
5. Convert all Adult Self Report and Child Behavioral Checklist scores of interest to numeric. As a note, child behavioral checklists are included so we can derive a set of models predicting g (i.e., outcome variable, not out-of-sample or longitudinal prediction) with both parental and child factors in there.
6. Initialize output vectors: we're running a looooong for-loop but initializing the values will make it slightly faster. Initialize them as 0s for 
  A) Linear vs. non-linear stat tests
  B) True maximum values for each subscale in the ABCD study (after QC)
  C) predictied derivatives, note we get an entire range of values for each iteration!
  D) predicted fits (from which the derivatives are derived)
  E) Maximum values for each subscale for each iteration
7. Now we can actually begin bootstrap. We will run this 10k instead of 1k times, so that we can calculate significance to a greater degree of accuracy.
8. Perform the bootstrap for each iteration. Recall that we need entire participants in/out for each iteration. We shouldn't be pulling one observation from one participant out and retaining the other inside of a bootstrap. We'll make a new dataframe by first randomly sampling participants, and then looping through each randomly sampled PTID and populating the bootstrap dataframe with that participant's data iteratively
9. Now we can derive our maximum values for each iteration
10. (labeled I in-script). Formally test for non-linearity. This fits a linear and a purely non-linear spline. Usually splines combine linear and non-linear fits, but here we distill solely the non-linear component to significance test whether or not there is a significant non-linear component to the fit. 
11. (labeled II in-script) Now we derive the spline fits for this bootstrap for each CBCL subscale with bam.
12. To interrogate each model fit, we examine the predicted value for each # of symptoms for each subscale. Each # of symptoms is reconstructed with seq by 1 to the maximum value of that bootstrap iteration
13. We need to make a pseudodataframe with the median age populated throughout to see what the rest of our model fit looks like given a constant age, as age is a major factor in cognitive scores in youth.
14. After fixing column names, we use the predict function, which takes in the fit model and the psuedodataframe
15. Now we printout the fit. Note we are only printing the the maximum value of each subscale in this iteration: we can't really extrapolate past that well. This will lead to 0's in iteration where participants with the highest cbcl subscale scores are not selected, which we will deal with in post-processing.
16. Now we will derive the derivatives of each fit, that is the pointwise slope of the fit across the range of number of reported symptoms.
17. Same spiel as the fit: print out to maximum bootstrap-specific value
18. saveout all the data as .rds files

## Upon completion, we should have gpBoots_asr.rds, gpDerivBoots_asr.rds, and gpFitBoots_asr.rds for use in the subsequent plotting script

## Plotting steps (assuming above .rds files are derived and locally available)

1. Chunk 1: load needed libraries.
2. Chunk 2 establishes some functions that we'll use throughout the markdown. One is to plot bootstraps with a fair amount of ggplot specifications, and the other is to find the furthest extent of symptoms reported in each bootstrap. For context: the significance testing we've employed uses the total number of bootstraps to derive p-values. If certain levels of symptoms are not included in some bootstraps, this changes how we can consider statistical significance. To maximize statistical certainty, only symptom ranges that were included in all bootstraps are considered. Therefore, we need this function `find_furthest_nonzero` to determine what the highest symptom count that was included across all all 10k iterations was. Finally, we set a color palette that will be used throughout.
3. Plot a basic child p-factor to parental p-factor relationship. ggMarginal is added to show the difference in the distribution of child symptoms and parental symptoms. Print out the correlation. :BEGIN UPDATES HERE AFTER SHERLOCK IS BACK ONLINE
4. The next chunk loads in the bootstrapped fits (`gParentpFitBoots.rds`) and plots the p-factor. We'll extract subfactors from the big all-subfactor-bootstrap dataframe. The maximum extent of each subfactor is the maximum value for that subfactor in masterdf. They are extracted accordingly: see the previous bootstrapping script (above) for reference. Next, we derive some plot elements, such as the threshold for borderline clinical and clinical thresholds from t-scores ([see this paper for reference](https://www.nature.com/articles/s41380-022-01522-w)).
5. The next chunk does the same, but extends this processing to internalizing and externalizing subscales (while extracting the information for subscales for future plotting).
6. Now we load in bootstrapped derivatives (`gParentpDerivBoots.rds`). The structure of which values belong to which subscale is the same as that used for the fits in steps 4-5. Now we'll do our manual significance testing at 99% confidence. Specifically, we'll find where derivatives (slopes) are positive or negative across all 10,000 iterations. This is a way of testing if the 99% confidence interval for the true slope excludes 0. We'll then plot significant locations (spans of symptoms where slope is significant) using geom_raster. This is equivalent to the approach used in [this paper](https://www.sciencedirect.com/science/article/pii/S1878929320300360), but manual bootstrapping was needed in this situation (hence sherlock/slurm). All credit for this idea to the esteemed Dr. Bart Larsen. We'll plot out the significant slopes for all subscales.
7. The last step is to summarize these derivatives. Obviously the models are non-linear so they don't have a single slope, but folks will probably prefer something immediately interpretable to something nuanced and in-the-weeds. Tertile splits (bottom third of the symptom range, top third) will allow us to order the different fits relative to each other as well. To summarize, we will split into the subclinical and the clinical range of symptoms. Note scales are ordered by their first tertile median slope (first of these plots) across both plots, and the color scale is the same as the rest of figure 1.
