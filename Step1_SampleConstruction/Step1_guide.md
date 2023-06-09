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

