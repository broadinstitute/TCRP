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

You may need to tweak the filepaths of the four input matrices to your situation. Also make sure you have the required libraries installed. You might notice that this code takes some time to run, as I (David) didn't optimize it. Briefly, this script defines all lineages to fewshot over, cell line to lineage mappings, gets the gene KOs to evaluate over (defined as gene KOs with at least one cell line that has CERES score of at least +- 6 std from the mean of all cell lines for that gene), removes expression features with bottom 10th percentile stdev, defines all genes with expression correlation > 0.4 to the expression of each gene KO (co-expression features as defined in the paper), and processes the CCLE mutations file to exclude silent mutations, one hot encode the remaining mutations into a feature matrix (as described in paper, 1 for a gene with >= 1 mutation, 0 otherwise), and exclude any gene with less than 10 cell lines with >= 1 mutation for that gene (as described in the paper). 

Then, you will need create a pkl file mapping gene KOs to their set of PPI partners (since TCRP paper included these as well), to include those features when predicting over. We have a script to do this that also includes a gene KO's paralogs and the gene KO itself, so its not exactly what is described in the TCRP paper, but it should be fairly close as a reproduction.

Our script requires use of our internal data management tool, Taiga, so you may not be able to run it. However, if you would like to try, run

```python get_ppi_features_dict.py```

Then, make the final dictionary mapping gene KOs to features to use for that gene KO. You can do this by running

```python generate_feature_dict_pkl.py```

This script simply collates all the data generated from the previous steps into a single dictionary mapping gene KOs to the list of features to use to simplify things for downstream scripts.

The code contains comments that should help you understand what is happening at each stage

The prepared to run matrices to run reproduce_achilles_results.py which reproduces the TCRP results can be downloaded [here]("https://drive.google.com/drive/folders/1Hyn65w7UyxCEsTUE2U1yk4JhcavAygR6?usp=sharing").

# Random Forest

## Reproducing results

The random forest model with the fewshot paradigm can be reproduced with the run_random_forest_with_feature_selection.py script. This runs one gene at a time, allowing parallization.

```python run_random_forest_with_feature_selection.py ***GENE OF INTEREST***```

The random forest model with a traditional train/test paradigm can be reproduced with the run_random_forest_with_feature_selection_regular_full_run.py script. This runs over all genes.

```python run_random_forest_with_feature_selection_regular_full_run.py```

These scripts run random forest models on the features filtered from the previous section (see "Data" section). Briefly, features are filtered in as close a way as possible to the Ideker paper.

# TCRP Model

## Modifications

* To get the TCRP model running locally we removed the CUDA dependency by deleting all .cuda() calls.
* We noticed an issue in data_loading.py which made the size of the training set used for meta-learning dependent on K. Basically, only K cell lines from each lineage were used which meant the results for different K values were not comparable. We modified data_loading.py and meta_learner_cv.py to add a fix_train_set_issue flag which allows us to fix this issue.
* In tcrp_cv.py the number of lineages selected for cross validation dependent on K. We added a fix_lineage_selection_issue flag which allows us select lineages based on min_lines_per_lineage instead. To reproduce the results we set min_lines_per_lineage to 15 similar to the paper.

## Issues

* When K = 1 the first meta-learner training epoch is always selected. This is becasue train_corr is always -1 since correlation is undefined for the single example in unseen_train_loader. We didn't correct this issue but it likely decreases the performance of 1-shot learning.
* In data_loading.py the number of cell lines from each lineage selected for training is dependent on K. This issue makes K values less comparable so we added a fix_train_set_issue flag to fix it. However, testing revealed that the impact was minimal so we didn't correct this issue when we reproduced the results.


## Reproducing results

The few-shot learning results presented in Fig. 2a can be reproducing with the reproduce_achilles_results.py script. This runs one or more genes at a time, allowing parallization.

```python reproduce_achilles_results.py --fix_lineage_selection_issue True --trials_for_each_K 5 --genes ***GENES OF INTEREST***```

# Results

CSV outputs from these scripts can be found in **/reproduce/results/** and **/reproduce/tcrp_achilles_analysis.md** includes some basic plots generated by tcrp_achilles_analysis.Rmd
