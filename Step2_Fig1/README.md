# Figure 1: Linear fits to cbcl subscales

Here, the goal is to submit 10k iterations of linear fits of g to cbcl subscales on sherlock (our computing cluser, slurm-based) across symtpom tertiles. So we'll take masterdf from the sample construction step and run 10,000 [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) to gauge model consistency over different permutations of participants. Because we have repeated measures (2x), we need to take care to remove entire participants for each bootstrap, not just observations.

These analyses take place over three scripts: [sample construction](https://github.com/WilliamsPanLab/gs/tree/master/Step1_SampleConstruction) (this is the previous step), [Fig1_Boots.R](https://github.com/WilliamsPanLab/gs/blob/master/Slurm/Fig1_Boots_CvSC.R), and then plotting the values derived from the bootstrapping ([Fig1.md/.rmd](https://github.com/WilliamsPanLab/gs/blob/master/Step2_Fig1/Figure1.md))

## Bootstrapping/computational steps (assuming sample construction is done)

1. Port masterdf over to your computer cluster (in Stanford's case, Sherlock) with the scp. Subsequent steps are in Fig1_Boots_CvSC.R, which is called to be sbatched (slurm equivalent of SGE qsub) by [this script](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/sbatch_Fig1.sh)
2. module load R/4.1 on sherlock and open R (terminal)
3. Load needed master dataframe
5. Glean number of subjects for bootstrapping purposes later
6. Convert all Child Behavioral Checklist scores of interest to numeric
7. Initialize output vectors: we're running a looooong for-loop but initializing the values will make it slightly faster. Initialize them as 0s for 
  A) Linear vs. non-linear stat tests
  B) True maximum values for each subscale in the ABCD study (after QC)
  C) predictied derivatives, note we get an entire range of values for each iteration!
  D) predicted fits (from which the derivatives are derived)
  E) Maximum values for each subscale for each iteration
8. Now we can actually begin bootstrap. We will run this 10k instead of 1k times, so that we can calculate significance to a greater degree of accuracy.
9. Perform the bootstrap for each iteration. Recall that we need entire participants in/out for each iteration. We shouldn't be pulling one observation from one participant out and retaining the other inside of a bootstrap. We'll make a new dataframe by first randomly sampling participants, and then looping through each randomly sampled PTID and populating the bootstrap dataframe with that participant's data iteratively
10. Now we can derive our maximum values for each iteration
11. (labeled I in-script). Formally test for non-linearity. This fits a linear and a purely non-linear spline. Usually splines combine linear and non-linear fits, but here we distill solely the non-linear component to significance test whether or not there is a significant non-linear component to the fit. 
12. (labeled II in-script) Now we derive the spline fits for this bootstrap for each CBCL subscale with bam.
13. To interrogate each model fit, we examine the predicted value for each # of symptoms for each subscale. Each # of symptoms is reconstructed with seq by 1 to the maximum value of that bootstrap iteration
14. We need to make a pseudodataframe with the median age populated throughout to see what the rest of our model fit looks like given a constant age, as age is a major factor in cognitive scores in youth.
15. After fixing column names, we use the predict function, which takes in the fit model and the psuedodataframe
16. Now we printout the fit. Note we are only printing the the maximum value of each subscale in this iteration: we can't really extrapolate past that well. This will lead to 0's in iteration where participants with the highest cbcl subscale scores are not selected, which we will deal with in post-processing.
17. Now we will derive the derivatives of each fit, that is the pointwise slope of the fit across the range of number of reported symptoms.
18. Same spiel as the fit: print out to maximum bootstrap-specific value
19. And save the maxmimum bootstrap-specific value
20. saveout all the data as .rds files

## Upon completion, we should have gpBoots.rds, gpDerivBoots.rds, and gpFitBoots.rds for use in the subsequent plotting script

## Plotting steps (assuming above .rds files are derived and locally available)

1. load needed libraries.
2. Establish some functions that we'll use throughout the markdown. One is to plot bootstraps with a fair amount of ggplot specifications, and the other is to find the furthest extent of symptoms reported in each bootstrap. For context: the significance testing we've employed uses the total number of bootstraps to derive p-values. If certain levels of symptoms are not included in some bootstraps, this changes how we can consider statistical significance. To maximize statistical certainty, only symptom ranges that were included in all bootstraps are considered. Therefore, we need this function `find_furthest_nonzero` to determine what the highest symptom count that was included across all all 10k iterations was. Finally, we set a color palette that will be used throughout.
3. Caclulate AIC in linear vs. nonlinear models as further confirmation of nonlinear relationship between g~symptoms across subscales.
4. Plot our reference linear model. This is equivalent to the reported correlation. The linear fit is also extracted to plot as a reference on future plots. `basic` is the name of the first plot, which is a hexplot.
5. The next chunk loads in the bootstrapped fits (`gpFitBoots.rds`) and plots the p-factor. We'll extract subfactors from the big all-subfactor-bootstrap dataframe. The maximum extent of each subfactor is the maximum value for that subfactor in masterdf. They are extracted accordingly: see the previous bootstrapping script (above) for reference. Next, we derive some plot elements, such as the threshold for borderline clinical and clinical thresholds from t-scores ([see this paper for reference](https://www.nature.com/articles/s41380-022-01522-w)).
6. The next chunk does the same, but extends this processing to internalizing and externalizing subscales (while extracting the information for subscales for future plotting).
7. Now we derive separate models for the clinical and subclinical portions of the sample. Effect sizes are evaluated via r^2 and compared. Both models are plotted onto the same background as that generated in the previous step.
8. Note the is an animation generated (chunk starting with ggplot2, gganimate, and dplyr library calls). No need to install these libraries if you don't need the .gif.
9. Now clinical vs. subclinical bootstrapped linear fits are loaded in for the p-factor, internalizing, and externalizing. This was threaded over 5 jobs (iterations 1-2000, 2001-4000, etc.), so 5 files were generated and need to be loaded in and merged.
10. Betas for linear estimates of g~p in subclinical and clinical boxplots are generated. Then for g~Internalizing, and g~Externalizing. Significance values are visualized as red dots (observed) vs. 10,000 null permutation values. Note externalizing does not exceed conventional significance value.
11. Now we load in bootstrapped fits (`gpFitBoots.rds`). The structure of which values belong to which subscale is specified manually. Now we'll plot the fits (using functions specified at the top), and do our manual significance testing of derivatives at 95/99% confidence (note both are depicted in the manuscript). Specifically, we'll find where derivatives (slopes) are positive or negative across all 10,000 iterations. This is a way of testing if the 99% confidence interval for the true slope excludes 0. We'll then plot significant locations (spans of symptoms where slope is significant) using geom_raster. This is equivalent to the approach used in [this paper](https://www.sciencedirect.com/science/article/pii/S1878929320300360), but manual bootstrapping was needed in this situation (hence sherlock/slurm). All credit for this idea to the esteemed Dr. Bart "Benevolent" Larsen. We'll plot out the significant slopes for all subscales. 95% confidence is predicated on the >9500 line: change to >9900 for 99% confidence intervals.
