# Introduction

This repository is an attempt to reproduce the results presented in

Jianzhu Ma, Samson H. Fong, Christopher J. Bakkenist, John Paul Shen, Soufiane Mourragui, Lodewyk F. A. Wessels, Marc Hafner, Roded Sharan, Jian Peng, Trey Ideker.  *Learning predictive models of drug response that translate across biological contexts. Nature Cancer*.

Specifically, we attempt to reproduce Fig. 2a from the paper which presents the correlation between predicted and actual CRISPR dependency as a function of the number of samples from the target tissue used for few-shot learning.

# Data

There are four into datafiles used to prepare the dataset. 
All downloads can be obtained from the depmap portal at the downloads page, which is here: https://depmap.org/portal/download/
They are enumerated as follows:

1. ceresgeneeffects.csv which is the CEREs data matrix, which can be downloaded at the depmap portal under the dataset name "Achilles Avana Public 17Q4 v2" on the left navigation bar. It is also available with this link: https://ndownloader.figshare.com/files/9562258
2. sample_info.csv which contains cell line information about the CCLE 19Q4 release. This can be downloaded from the depmap portal under the "DepMap Public 19Q4" dataset on the left navigation bar. It is also available with this link: https://ndownloader.figshare.com/files/20274744
3. CCLE_expression.csv which contains cell line gene expression information from the CCLE 19Q4 release. This can be downloaded from the depmap portal under the "DepMap Public 19Q4" dataset on the left navigation bar. It is also available with this link: https://ndownloader.figshare.com/files/20234346
4. CCLE_mutations.csv which contains cell line mutation information from the CCLE 19Q4 release. This can be downloaded from the depmap portal under the "DepMap Public 19Q4" dataset on the left navigation bar. It is also available with this link: https://ndownloader.figshare.com/files/20274747

To prepare the data matricies as input to the Random Forest and TCRP pipelines, please run 

```python preprocess_data.py```

You may need to tweak the filepaths of the four input matrices to your situation. Also make sure you have the required libraries installed. You might notice that this code takes some time to run, as I (David) didn't optimize it. 

Then, you will need create a pkl file mapping gene KOs to their set of PPI partners, to include those features when predicting over. We have a script to do this that also includes a gene KO's paralogs, so its not exactly what is described in the TCRP paper, but it should be fairly close as a reproduction.

Our script requires use of our internal data management tool, Taiga, so you may not be able to run it. However, if you would like to try, run 

```python get_ppi_features_dict.py```

Then, make the final dictionary mapping gene KOs to features to use for that gene KO by running

```python generate_feature_dict_pkl.py```

The code contains comments that should help you understand what is happening at each stage


Data to run reproduce_achilles_results.py can be downloaded
# TCRP Model

## Modifications

* To get the TCRP model running locally we removed the CUDA dependency by deleting all .cuda() calls
* We noticed an issue in data_loading.py which made the size of the training set used for meta-learning dependent on K. Basically, only K cell lines from each lineage were used which meant the results for different K values were not comparable. We modified data_loading.py and meta_learner_cv.py to add a fix_train_set_issue flag which allows us to fix this issue.
* We noticed an issue in tcrp_cv.py which made the number of lineages selected for cross validation dependent on K. Again, this issue meant the results for different K value were not comparable. We added a fix_lineage_selection_issue flag to our code which allows us to fix this issue.

## Reproducing results

The random forest model with the fewshot paradigm can be reproduced with the run_random_forest_with_feature_selection.py script. This runs one gene at a time, allowing parallization.

```python run_random_forest_with_feature_selection.py ***GENE OF INTEREST***```

The random forest model with a traditional train/test paradigm can be reproduced with the run_random_forest_with_feature_selection_regular_full_run.py script. This runs over all genes.

```python run_random_forest_with_feature_selection_regular_full_run.py```

The few-shot learning results presented in Fig. 2a can be reproducing with the reproduce_achilles_results.py script.

```python reproduce_achilles_results.py --fix_lineage_selection_issue True --fix_train_set_issue True --trials_for_each_K 5 --genes HNF1B ESR1```
