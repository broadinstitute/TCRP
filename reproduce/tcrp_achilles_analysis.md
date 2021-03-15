TCRP Achilles Analysis
================
William Colgan
15 March 2021

# Test correlation by K

Plots includes 464 genes with at least one 6 sigma outlier cell line.
Points represent average across 5 trials for TCRP and 10 trials for few
shot random forest for each gene and K
value.

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# HNF1B

The paper reports HNF1B prediction performance of .6 for TCRP and .19
for random
forest

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Performance relative to standard random forest

Standard random forest uses the same features but uses standard
five-fold cross validation. K value is fixed at 5 for both TCRP and few
shot random
forest.

## Few shot vs standard RF

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## TCRP vs standard RF

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
