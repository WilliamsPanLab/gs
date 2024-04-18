# Sample Construction

Sample construction is arguably the most important step. Here, we will be parsing participants from timepoints 1 and 2 who have complete data on our variables of interest. 

It's important to note that these measurements, even from past timepoints, can change over time in NDAR. Different releases have different kinds of and extents of missing data. This is uncommon, but can induce minor month-to-month changes in exactly which subjects have missing data when it occurs. Further, some measures are available from multiple instruments, with variable missingness and variable missingness over time.

The easiest way to follow this walkthrough might be to download the .Rmd to follow the code on your local machine, and follow the .md file online. 

## The steps - organized by "chunks" in the Rmd

1. After loading the libraries we'll need, load in our primary measures of interest: child mental health. Specifically, the [Child Behavioral Checklist](https://nda.nih.gov/data_structure.html?short_name=abcd_cbcls01). Some objects are initialized here, and children are checked for having missing data at either timepoint. For maximal interpretability of analyses, children without data at one timepoint are excluded.

2. Chunk 2 takes in the family ID of each participant. We also need this for maximal interpretability, to avoid overrepresentation of certain families in our results. The [acspw03](https://nda.nih.gov/data_structure.html?short_name=acspsw03) instrument is used here. 

3. Chunk 3 just converts data to numeric format.

4. Chunk 4 loads in cognitive data. Note this and subsequent chunks are modeled after a [nice paper](https://pubmed.ncbi.nlm.nih.gov/30595399/), for which we explicitly demonstrate functional equivalence as a supplementary result. The basis of the derived cognitive score for each child (at each timepoint) is largely the NIH toolbox, but two additional tasks are included. In chunk 6 the cognitive data is merged in ([source 1](https://nda.nih.gov/data_structure.html?short_name=abcd_tbss01), [2](https://nda.nih.gov/data_structure.html?short_name=abcd_ps01), and [3](https://nda.nih.gov/data_structure.html?short_name=lmtp201)). Note not all tasks gathered at timepoint1 were repeated at timepoint2: we used data/tasks available at both timepoints.

5. Chunk 5 preps cognitive data for factorization: converts data to numeric, calculates correct trials, and notes + excludes participants with missing data (unfortunately we do lose a lot of kids here)

6. Chunk 6 further preps for factorization by randomly removing one member from each family with more than one child represented in the study. Love that full-on families are completing ABCD, but between later results and deriving *g*, we don't want these dedicated participant-families to drive results more than single-participant families.

7. Chunk 7 runs PCA to derive *g* (first factor). As noted, a supplementary analysis demonstrates this produces an equivalent outcome to more complex approaches.

8. Chunk 8 just checks for complete data at both timepoints.

9. Chunk 9 handles a participants tsv, which is [a useful file with centralized info](https://collection3165.readthedocs.io/en/stable/recommendations/#2-the-bids-participants-files-and-matched-groups). We use this one to derive sex and parent income, as child-endorsed gender identification is missing for every single observation for every single child (as of now, from KSADS). Same drill as prior: kids with missing data are excluded. We also save out the primary dataframe here (gp_masterdf.rds)

10. Chunk 13 imports higher-resolution demographic data for robustness analyses conducted in R1. Note in R1 there was one kid with extreme outlier data (> 7 S.D. from mean in otherwise normally distributed g score) was omitted.
   
11. Chunk 11, then saves out the master dataframe to be bootstrapped in subsequent steps.

12. Chunk 12 calculates and saves out tertile thresholds from participant t-scores on the cbcl.

13. Chunk 13 loads in puberty data (where available) for robustness R1 analyses

14. Chunk 14 loads in BPM data (where available) for robustness R1 analyses

15. Chunk 15 evaluates number of observations retained at each site. Allows us to predict the most populated site for least artificially induced variability in subsequent bootstraps

16. Chunks 16 and 17 plots missing data as alluvial and pie charts (before and after)

## gp_masterdf.rds from step 11 is the starting point for subsequent steps
