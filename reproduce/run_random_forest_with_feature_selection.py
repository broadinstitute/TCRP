import collections
import pandas as pd
import pickle
import random
import sys
from scipy.special import comb
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor

# Generate 10 random splits of cell lines for a tissue to fewshot over
def get_10_splits(cell_lines_in_tissue, num_cell_lines_in_training_set):
    set_of_combs = set()
    while len(set_of_combs) < min(comb(len(cell_lines_in_tissue), num_cell_lines_in_training_set, exact=True), 10):
        random_training_set = tuple(random.sample(cell_lines_in_tissue, num_cell_lines_in_training_set))
        set_of_combs.add(random_training_set)
    return set_of_combs

# Convert defaultdict to regular for pickling
def default_to_regular(d):
    if isinstance(d, collections.defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def main():
    logfile = open("log.txt", "w+")

    feature_matrix = pd.read_csv("achilles_data/feature_matrix.csv", header=0, index_col=0)
    labels_matrix = pd.read_csv("achilles_data/labels_matrix.csv", header=0, index_col=0)

    # K shot learning to perform
    FEWSHOT_START = 0
    FEWSHOT_END = 10

    # get gene to predict over
    argv = sys.argv[1:]
    if len(argv) != 1:
        raise ValueError('A very specific bad thing happened.')

    gene = argv[0]

    # Read in necessary files
    with open('achilles_data/lineages_to_fewshot_over.pkl', 'rb') as f:
        lineages_to_fewshot_over = pickle.load(f)

    with open('achilles_data/gene_to_features.pkl', 'rb') as f:
        gene_to_features = pickle.load(f)

    # Get feature set
    final_features_to_use = gene_to_features[gene]
    logfile.write(f"Number of features used: {len(final_features_to_use)}\n")

    # Quit if gene has no features
    print(f"Gene: {gene} has {len(final_features_to_use)} features!")
    if len(final_features_to_use) == 0:
        quit()

    # prepare feature matrix
    final_features_to_use = list(final_features_to_use)
    feature_matrix = feature_matrix[final_features_to_use]

    all_cell_lines = feature_matrix.index.tolist()

    # Dictionary mapping tissue to gene to list of k shot performances
    tissue_to_gene_to_corrlation_list = collections.defaultdict(lambda : collections.defaultdict(lambda : collections.defaultdict(list)))

    # Replicate TCRP fewshot learning process on random forest
    for tissue in lineages_to_fewshot_over:
        sys.stderr.write(f"Tissue: {tissue}")
        cell_lines_in_tissue = lineages_to_fewshot_over[tissue]
        for num_cell_lines_in_training_set in range(FEWSHOT_START, FEWSHOT_END + 1):
            cell_lines_to_include_in_training = get_10_splits(cell_lines_in_tissue, num_cell_lines_in_training_set)
            for cell_line_set in cell_lines_to_include_in_training:
                test_set_cell_lines = [element for element in cell_lines_in_tissue if element not in cell_line_set]
                training_set_cell_lines = [element for element in all_cell_lines if element not in test_set_cell_lines]

                if len(set(test_set_cell_lines) & set(training_set_cell_lines)) != 0:
                    raise ValueError('A very specific bad thing happened.')

                train_features = feature_matrix.loc[training_set_cell_lines]
                test_features = feature_matrix.loc[test_set_cell_lines]

                train_labels = labels_matrix.loc[training_set_cell_lines][gene]
                test_labels = labels_matrix.loc[test_set_cell_lines][gene]

                model = RandomForestRegressor(n_estimators=100, max_depth=8, min_samples_leaf=5, n_jobs=-1, random_state=0)
                model.fit(train_features, train_labels)
                predictions = model.predict(test_features)

                correlation = pearsonr(test_labels.tolist(), list(predictions))[0]
                tissue_to_gene_to_corrlation_list[tissue][gene][num_cell_lines_in_training_set].append(correlation)

    tissue_to_gene_to_corrlation_list = default_to_regular(tissue_to_gene_to_corrlation_list)

    # ouput results
    with open(f'{gene}_results_dict.pkl', 'wb') as handle:
        pickle.dump(tissue_to_gene_to_corrlation_list, handle)


main()
