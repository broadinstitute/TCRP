# Introduction

This repository is an attempt to reproduce the results presented in

Jianzhu Ma, Samson H. Fong, Christopher J. Bakkenist, John Paul Shen, Soufiane Mourragui, Lodewyk F. A. Wessels, Marc Hafner, Roded Sharan, Jian Peng, Trey Ideker.  *Learning predictive models of drug response that translate across biological contexts. Nature Cancer*.

Specifically, we attempt to reproduce Fig. 2a from the paper which presents the correlation between predicted and actual CRISPR dependency as a function of the number of samples from the target tissue used for few-shot learning.

# Data

TODO: Add info about data processing

Data to run reproduce_achilles_results.py can be downloaded
# TCRP Model

## Modifications

* To get the TCRP model running locally we removed the CUDA dependency by deleting all .cuda() calls
* We noticed an issue in data_loading.py which made the size of the training set used for meta-learning dependent on K. Basically, only K cell lines from each lineage were used which meant the results for different K values were not comparable. We modified data_loading.py and meta_learner_cv.py to add a fix_train_set_issue flag which allows us to fix this issue.
* We noticed an issue in tcrp_cv.py which made the number of lineages selected for cross validation dependent on K. Again, this issue meant the results for different K value were not comparable. We added a fix_lineage_selection_issue flag to our code which allows us to fix this issue.

## Reproducing results

The few-shot learning results presented in Fig. 2a can be reproducing with the reproduce_achilles_results.py script.

```python reproduce_achilles_results.py --fix_lineage_selection_issue True --fix_train_set_issue True --trials_for_each_K 5 --genes HNF1B ESR1```
