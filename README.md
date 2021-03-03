# qPCRTools
R package to help analyze qPCR data

Installation

```{r}
install.packages("devtools")
devtools::install_github("kevincjnixon/qPCRTools")
```

Usage - ddCt method

```{r}
library(qPCRTools)
easyRT() #Run interactively
easyRT(showEB=F, showStat=F) #Run interactively, but don't plot error bars (showEB=F) or stats (showStat=F)
```

When running interactively:
1. File browser will pop up, choose text-delimited file to analyze
2. You will be prompted if this is 'bioRad' input, enter Y/N (see below for details)
3. Select a standard deviation threshold to filter Ct values (if SD of triplicates is greater than this threshold, an outlier will be selected and removed from analysis)
4. Levels of sample indicators will appear. You should use format 'Condition*delim*Replicate' in your setup. Indicate what *delim* is. Leave black if it is a space.
5. Select your reference condition
6. Indicate how many reference genes you used (at least 1). More than 1 reference gene will use geoMean of Ct values of those genes for reference.
7. Indicate your reference gene(s)
8. Indicate how you want to average Ct values (geoMean = geometric mean or mean)
9. Enter a title for the plot
10. Choose if you want to write out results to .csv (Y/N)
  a. if Y, enter filename (ending in .csv)
    i. if file exists, overwrite? (Y/N) - 'Y' will overwrite, 'N' will offer to select new filename
  b. if N, you're done

Plot will be produced with showing conditions as colours on bars, plotting genes on x-axis and relative gene expression (delta-delta Ct method) on y-axis. Error bars are SEM (if showEB=T; default). Stats are from Benjamini-Hochberg-corrected pairwise-t-tests (if showStat=T; default).

Setup for data input:
Tab-delimited files output from qPCR machine.

bioRad CFX machine: (Select "Y" at first prompt)

| Well | Fluor | Target | Content | Sample | Biological Set Name |  Cq  | Cq Mean | Cq Std. Dev |
|------|-------|--------|---------|--------|---------------------|------|---------|-------------|
| A01  | SYBR  | Gene   | Unkn    | WT-1   |                     | 24.37| 24.37   |      0      |

There can be more columns, but it is necessary to have the first 8 in this order. There is no header, and no blank first column.

Other format: (Select "N" at first prompt)

*Header 10 lines long*

| Position |  Flag | Sample | Detector | Task |  Ct  |
|----------|-------|--------|----------|------|------|
|    A1    | Passed| WT-1   | Gene     | Unkn | 24.37|

There can be more columns, but it is necessary to have these columns in this order with these names

Note that 'Sample' in both cases is in the format 'Condition*delimiter*replicate' (i.e. WT is the condition, '-' is the delimiter, and 1 is the replicate)

There should be technical triplicates for each gene and each condition/replicate, and at least two biological replicates per condition.
