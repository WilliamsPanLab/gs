# Sample Construction

Sample construction is the most important step! Here, we will be parsing participants from timepoints 1 and 2 who have complete data on our variables of interest. 

It's important to note that these measurements, even from past timepoints, can change over time in NDAR. Different releases have different kinds of and extents of missing data. This is uncommon, but can induce minor month-to-month changes in exactly which subjects have missing data when it occurs. Further, some measures are available from multiple instruments, with variable missingness and variable missingness over time.

You can download the .Rmd to follow the code on your local machine, and follow the .md file online. 

## The steps

1. First, we'll load the libraries we need (chunk 1).

2. Next, we'll load in our primary measures of interest: child mental health (chunk 2). Specifically, the [Child Behavioral Checklist](https://nda.nih.gov/data_structure.html?short_name=abcd_cbcls01). Some objects are initialized here, and children are checked for having missing data at either timepoint. For maximal interpretability of analyses, children without data at one timepoint are excluded.

3. Chunk 3 processes adult mental health data, that is, the parents. Specifically, the [Parent Adult Self Report Raw Scores Aseba](https://nda.nih.gov/data_structure.html?short_name=pasr01) instrument is used. Similarly, participants without data at both timepoints are exlcuded. Very few kids are lost here.

4. Chunk 4 takes in the family ID of each participant. We also need this for maximal interpretability, to avoid overrepresentation of certain families in our results. The [acspw03](https://nda.nih.gov/data_structure.html?short_name=acspsw03) instrument is used here.

5. Chunk 5 just converts data to numeric format

6. Chunk 6 loads in cognitive data. Note this and subsequent chunks are modeled after a [nice paper](https://pubmed.ncbi.nlm.nih.gov/30595399/), for which we explicitly demonstrate functional equivalence as a supplementary result. The basis of the derived cognitive score for each child (at each timepoint) is largely the NIH toolbox, but two additional tasks are included. In chunk 6 the cognitive data is merged in ([source 1](https://nda.nih.gov/data_structure.html?short_name=abcd_tbss01), [2](https://nda.nih.gov/data_structure.html?short_name=abcd_ps01), and [3](https://nda.nih.gov/data_structure.html?short_name=lmtp201))

7. Chunk 7 preps cognitive data for factorization: converts data to numeric, calculates correct trials, and notes + excludes participants with missing data (unfortunately we do lose a lot of kids here)

8. Chunk 8 further preps for factorization by randomly removing one member from each family with more than one child represented in the study. Love that full-on families are completing ABCD, but between later results and deriving *g*, we don't want these dedicated participant-families to drive results more than single-participant families.

9. Chunk 9 runs PCA to derive *g* (first factor). As noted, a supplementary analysis demonstrates this produces an equivalent outcome to more complex approaches.

10. Chunk 10 just checks for complete data at both timepoints.

11. Chunk 11 calculates the adult p factor. For multiple reasons, discussed in-text, we are using the sum of endorsed symptoms. Several items are reverse scored, which presents a potentially challenging interpretative framework. Specifically, it's not clear that saying you are happy is equivalent to negative psychiatric symptoms. To make the fewest assumptions, these items are omitted rather than reversed in adult totals (subtracted from the total sum).

12. Chunk 12 handles a participants tsv, which is [a useful file with centralized info](https://collection3165.readthedocs.io/en/stable/recommendations/#2-the-bids-participants-files-and-matched-groups). We use this one to derive sex and parent income, as child-endorsed gender identification is missing for every single observation for every single child (as of now, from KSADS). Same drill as prior: kids with missing data are excluded. We also save out the primary dataframe here (gp_masterdf.rds)

13. Chunk 13 plots missing data. It's time to see where we lost participants. This makes an alluvial plot.

14. Chunk 14 plots missing data as pie charts (before and after)

15. Loads in [ksads data](https://nda.nih.gov/data_structure.html?short_name=abcd_ksad501) so we have a child-report version of their p factor score. Note the ksads is messy: branching logic, variably-expressed missing data, and big (~1,000 variables). These data are also saved out after merging with masterdf, and further missingness is plotted.

16. This chunk prepares data for temporal precedence analysis: largely a replication of Romer & Pizzagalli [2021](https://pubmed.ncbi.nlm.nih.gov/34332330/). This uses a sep. dataframe because we are explicitly interested in temporal precedence from timepoint 1 to timepoint 2. This does not establish causality, but does add credence to the ways our models are constructed through the rest of the paper.

## Upon completion, we should have masterdf, masterdf2, and OutDFTmpPrec to bootstrap on our institution's computer cluster