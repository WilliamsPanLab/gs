# Figure 1: Linear fits to cbcl subscales

Here, the goal is to submit 10k iterations of linear fits of g to cbcl subscales on sherlock (our computing cluser, slurm-based) across symtpom tertiles. So we'll take masterdf from the sample construction step and run 10,000 [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) to gauge model consistency over different permutations of participants. Because we have repeated measures (2x), we need to take care to remove entire participants for each bootstrap, not just observations.

These analyses take place over three scripts: [sample construction](https://github.com/WilliamsPanLab/gs/tree/master/Step1_SampleConstruction) (this is the previous step), [Fig1_Boots.R](https://github.com/WilliamsPanLab/gs/blob/master/Slurm/Fig1_Boots_CvSC.R), and then plotting the values derived from the bootstrapping ([Fig1.md/.rmd](https://github.com/WilliamsPanLab/gs/blob/master/Step2_Fig1/Figure1.md))

## Bootstrapping/computational steps (assuming sample construction is done)

1. Port masterdf over to your computer cluster (in Stanford's case, Sherlock) with the scp. Subsequent steps are in Fig1_Boots_CvSC.R, which is called to be sbatched (slurm equivalent of SGE qsub) by [this script](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/sbatch_Fig1.sh)
2. module load R/4.1 on sherlock and open R (terminal)
3. Load needed master dataframe
4. Glean number of subjects for bootstrapping purposes later
5. Convert all Child Behavioral Checklist scores of interest to numeric
6. Initialize output vectors: we're running a looooong for-loop but initializing the values will make it slightly faster. Initialize them as 0s for 
  A) Difference in betas for g~symptoms for total problems, internalizing, externalizing
  B) Actual betas derived for the three relationships
  C) Betas in the high symptom tertile for three relaitonships
  D) Betas in the low symptom tertile for three relaitonships
7. Derive the non-bootstrapped, full data value betas (full, high and low symptoms) for our best guess of the centerpoint across the three symptom classes.
8. Use these to get the true beta differences across tertiles
9. Now we can actually begin bootstrap. We will run this 10k instead of 1k times, so that we can calculate significance to a greater degree of accuracy.
10. Perform the bootstrap for each iteration. Recall that we need entire participants in/out for each iteration. We shouldn't be pulling one observation from one participant out and retaining the other inside of a bootstrap. We'll make a new dataframe by first randomly sampling participants, and then looping through each randomly sampled PTID and populating the bootstrap dataframe with that participant's data iteratively
11. Now we get the number of participants above and below our symptom tertile thresholds (derived from full sample) for accurate resampling of pseudo-tertiles (matched n).
12. Randomly assign to pseudotertiles based on n above and below (non-overlapping)
13. Derive beta coefficients for full bootstrap resample as well as in pseudotertiles
14. Calculate difference in pseudotertile betas for null distribution
15. Get actual tertile beta coefficients
16. Save out true values as 10,001st in beta differences
17. saveout all the data as .rds files

## Upon completion, we should have gp_CvSC_diffs.rds for g~symptom betas across full samples, tertiles, pseudotertiles, and real and permuted beta differences for our three symptom classes. Full symptom class granularity is commented out but available.

## Plotting steps (assuming above .rds files are derived and locally available)

1. load needed libraries.
2. Establish some functions that we'll use throughout both markdowns. One is to plot bootstraps with a fair amount of ggplot specifications, and the other is to find the furthest extent of symptoms reported in each bootstrap. For context: the significance testing we've employed uses the total number of bootstraps to derive p-values. If certain levels of symptoms are not included in some bootstraps, this changes how we can consider statistical significance. To maximize statistical certainty, only symptom ranges that were included in all bootstraps are considered. Therefore, we need this function `find_furthest_nonzero` to determine what the highest symptom count that was included across all all 10k iterations was. Finally, we set a color palette that will be used throughout.
3. Caclulate AIC in linear vs. nonlinear models as further confirmation of nonlinear relationship between g~symptoms across subscales.
4. Plot our reference linear model. This is equivalent to the reported correlation. The linear fit is also extracted to plot as a reference on future plots. `basic` is the name of the first plot, which is a hexplot. We also want to construct equivalent plots across symtpom tertiles.
5. Calculate linear model stats for reporting.
6. Bonus chunk not in publication: create a .gift demonstrating the dependence of the g~symptom fit on the symptom range being sampled.
7. Load in bootstrapped fits.
8. For total problems, internalizing, and externalizing, plot both the distributions of low-symptom and high-symptom tertile fits as well as the true difference in betas vs. the null difference in betas (null distribution from psuedotertile differences).
9. Significance testing of pseudotertiles is at end. If there were bootstrap iterations where the difference in slopes between pseudotertiles exceeded the true difference, the commented out lines (533-536) would calculate the p-value. Because true > null in all instances, we know p is <0.0001 in this operationalization.
