# This code is all CDS code
# Author: David Wu (dwu@broadinstitute.org)

import collections
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy.stats import pearsonr
from tqdm import *


K_FOLD_NUMBER = 5

# Function to evaluate the performance via mse, rmse, mae, and r2 of a given set of predicted output variables vs the actual ground truth values
def evaluate(x, y):
    predictions_list = list()
    ground_truth_list = list()
    
    # Intialize k-fold cv
    kfold = KFold(K_FOLD_NUMBER, shuffle=True, random_state=0)
    for train, test in kfold.split(x):
        # Split dataset into train-test splits
        train_dataset = x.iloc[train]
        train_labels = y.iloc[train]
        
        test_dataset = x.iloc[test]
        test_labels = y.iloc[test]
        
        # Train test model
        model = RandomForestRegressor(n_estimators=100, max_depth=8, min_samples_leaf=5, n_jobs=-1, random_state=0)
        model.fit(train_dataset, train_labels)
        
        predictions = model.predict(test_dataset)
        predictions_list += list(predictions)

        ground_truth_list += test_labels.tolist()
        
    # Return overall performance
    return pearsonr(predictions_list, ground_truth_list)[0]


def main():
    # Read input data
    feature_matrix = pd.read_csv("achilles_data/feature_matrix.csv", header=0, index_col=0)
    labels_matrix = pd.read_csv("achilles_data/labels_matrix.csv", header=0, index_col=0)

    # Read metadata
    with open('achilles_data/gene_to_features.pkl', 'rb') as f:
        gene_to_features = pickle.load(f)

    all_genes_df = pd.read_csv('achilles_data/to_run_genes_feature_selection.csv', header=0)
    all_genes = all_genes_df['Gene'].tolist()

    # Run for each gene, 5 fold cv and get correlation across all 5 folds
    gene_to_performance = dict()
    for gene in tqdm(all_genes):
        final_features_to_use = gene_to_features[gene]
        
        if len(final_features_to_use) == 0:
            continue

        final_features_to_use = list(final_features_to_use)
        final_feature_matrix = feature_matrix[final_features_to_use]

        correlation = evaluate(final_feature_matrix, labels_matrix[gene])
        gene_to_performance[gene] = correlation


    with open('results_dict.pkl', 'wb') as handle:
        pickle.dump(gene_to_performance, handle)


main()
