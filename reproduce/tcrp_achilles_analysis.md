TCRP Achilles Analysis
================
William Colgan
26 February 2021

# Purpose

To analyze the results from our attempt to reproduce Fig. 2a of
*Learning predictive models of drug response that translate across
biological contexts. Nature Cancer*

# Test correlation by K

Plots includes 464 genes that kill at least one cell line. Points
represent average across 5 trials for TCRP and 10 trials for few shot
random forest for each gene and K
value

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# HNF1B

The paper claims an HNF1B prediction performance of .6 for TCRP and .19
for random
forest.

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Performance relative to standard random forest

Standard random forest uses the same features but uses standard cross
validation. K value is fixed at 5 for both TCRP and few shot random
forest.

## Few shot vs standard RF

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## TCRP vs standard RF

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

![](tcrp_achilles_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
