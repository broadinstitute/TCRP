# This code is all CDS code
# Author: David Wu (dwu@broadinstitute.org)

import pandas as pd
import pickle


def main():
    # Read in feature matrix
    feature_matrix = pd.read_csv('achilles_data/feature_matrix.csv', header=0, index_col=0)

    # Read in coexpression dict
    with open('achilles_data/gene_to_coexpression_features.pkl', 'rb') as f:
        gene_to_coexpression_features = pickle.load(f)

    # Read in ppi dictionary
    with open('achilles_data/ppi_features.pkl', 'rb') as handle:
        ppi_features = pickle.load(handle)

    # Get all genes to predict over
    all_genes_df = pd.read_csv('achilles_data/to_run_genes_feature_selection.csv', header=0)
    all_genes = all_genes_df['Gene'].tolist()

    # Get features to use for this gene
    ppi_features_final = dict()
    for depmap_id in ppi_features:
        hgnc_symbol = depmap_id.split(' ')[0]
        ppi_features_final[hgnc_symbol] = ppi_features[depmap_id]

    # All initial features
    all_initial_features = set(feature_matrix.columns)
    
    # Population final dictionary mapping gene KO to features to use
    gene_to_features = dict()
    for gene in all_genes:
        final_features_to_use = set()
        if gene in ppi_features_final:
            for partner in ppi_features_final[gene]:
                expression_feature_name = partner + " EXPRESSION"
                if expression_feature_name in all_initial_features:
                    final_features_to_use.add(expression_feature_name)

                mutation_feature_name = partner + " MUTATION"
                if mutation_feature_name in all_initial_features:
                    final_features_to_use.add(mutation_feature_name)

        if gene in gene_to_coexpression_features:
            for partner in gene_to_coexpression_features[gene]:
                expression_feature_name = partner + " EXPRESSION"
                if expression_feature_name in all_initial_features:
                    final_features_to_use.add(expression_feature_name)

                mutation_feature_name = partner + " MUTATION"
                if mutation_feature_name in all_initial_features:
                    final_features_to_use.add(mutation_feature_name)

        gene_to_features[gene] = sorted(list(final_features_to_use))
    # END get features to use for this gene

    # Dump the final dict
    with open('achilles_data/gene_to_features.pkl', 'wb') as handle:
        pickle.dump(gene_to_features, handle)


main()
