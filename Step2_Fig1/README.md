# Figure 1: Non-linear fits to cbcl subscales

Here, the goal is to submit 10k iterations of non-linear fits of g to cbcl subscales on sherlock (our computing cluser, slurm-based). So we'll take masterdf from the sample construction step and run 10,000 [bootstraps](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) to gauge model consistency over different permutations of participants. Because we have repeated measures (2x), we need to take care to remove entire participants for each bootstrap, not just observations.

These analyses take place over three scripts: sample construction (done already), [Fig2_Boots_xSectional.R](https://github.com/WilliamsPanLab/gp/blob/master/Slurm/Fig2_Boots_xSectional.R), and then plotting the values derived from the bootstrapping ([Fig1.Rmd](https://github.com/WilliamsPanLab/gp/blob/master/Figures/code/Fig1.Rmd))

## The steps (assuming sample construction is done)

1. Port masterdf over to Sherlock with scp. Subsequent steps in Fig2_Boots_xSectional.R.
2. module load R/4.1 on sherlock and open R (terminal)
3. Load needed master dataframe
4. Glean number of subjects for bootstrapping purposes later
5. Convert all Child Behavioral Checklist scores of interest to numeric
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
18. And save the maxmimum bootstrap-specific value
19. saveout all the data as .rds files

## Upon completion, we should have gpBoots.rds, gpDerivBoots.rds, and gpFitBoots.rds for use in the subsequent plotting script
