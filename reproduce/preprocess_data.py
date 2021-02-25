import collections
import numpy as np
import pandas as pd
import pickle
from scipy.stats import zscore, pearsonr

# YOU MAY HAVE TO CHANGE THE BELOW LINES
achilles_dataseta_file = 'data/ceresgeneeffects.csv'
cell_line_info_file = 'data/sample_info.csv'
ccle_expression_file = 'data/CCLE_expression.csv'
ccle_mutation_file = 'data/CCLE_mutations.csv'
# YOU MAY HAVE TO CHANGE THE ABOVE LINES

# Read in CERES data
GENE_EFFECT_MATRIX = pd.read_csv(achilles_dataseta_file, header=0, index_col=0)
# Get all cell lines in CERES data
all_achilles_cell_lines = GENE_EFFECT_MATRIX.columns.tolist()

# Read in the sample info file and subset it to all CEREs cell lines
cell_line_info_df = pd.read_csv(cell_line_info_file, header=0)
cell_line_info_df = cell_line_info_df.loc[cell_line_info_df['CCLE_Name'].isin(all_achilles_cell_lines)]

# Remove cell lines where the lineage is not defined
cell_line_info_df = cell_line_info_df[cell_line_info_df['lineage'].notna()]
# Get all unique lineages
all_lineages = set(cell_line_info_df['lineage'].tolist())
print(f"Number of lineages included: {len(all_lineages)}")

# Dictionaries mapping between a cell line's DepMap_ID (used in CCLE files) and CCLE_Name (used in CEREs files)
ccle_id_to_depmap_id = dict()
depmap_id_to_ccle_id = dict()
# Read through sample info file to populate above dictionaries
for index, row in cell_line_info_df.iterrows():
    ccle_id = row["CCLE_Name"]
    depmap_id = row['DepMap_ID']
    if ccle_id not in ccle_id_to_depmap_id:
        ccle_id_to_depmap_id[ccle_id] = depmap_id
        depmap_id_to_ccle_id[depmap_id] = ccle_id
    else:
        print('ERROR!')
        quit()

# This next section outputs a pickled dictionary mapping all lineages to all their cell lines for later use
lineage_to_cell_lines = dict()
for lineage in all_lineages:
    cell_line_info_df_ = cell_line_info_df.loc[cell_line_info_df['lineage'] == lineage]
    lineage_to_cell_lines[lineage] = set(cell_line_info_df_['DepMap_ID'])

with open('achilles_data/lineage_to_cell_lines.pkl', 'wb') as handle:
    pickle.dump(lineage_to_cell_lines, handle)
# End section

# Dictionary containing the lineages of cell lines with >= 15 cell lines. These are the cell lines we will few shot learn over
lineages_to_few_shot_over = dict()
for lineage in all_lineages:
    # Subset dataframe containing sample information to each lineage, and check how many rows (cell lines) are in each lineage
    cell_line_info_df_ = cell_line_info_df.loc[cell_line_info_df['lineage'] == lineage]
    if cell_line_info_df_.shape[0] >= 15:
        lineages_to_few_shot_over[lineage] = set(cell_line_info_df_['DepMap_ID'])

# Rename CERES matrix cell line names to their DepMap_IDs
#-------------------------------- GENE EFFECT MATRIX -------------------------------------------------------------------------------------------------------------------#
GENE_EFFECT_MATRIX.columns = [ccle_id_to_depmap_id[element] for element in GENE_EFFECT_MATRIX.columns]
all_cell_lines_to_include = GENE_EFFECT_MATRIX.columns.tolist()

GENE_EFFECT_MATRIX = GENE_EFFECT_MATRIX.T

# Get a csv listing all genes to run 
genes_to_run_df_dict = collections.defaultdict(list)
GENE_EFFECT_MATRIX_Z_SCORE = GENE_EFFECT_MATRIX.apply(zscore)
all_genes_to_test_over = list()
for gene in GENE_EFFECT_MATRIX_Z_SCORE.columns:
    column_list = GENE_EFFECT_MATRIX_Z_SCORE[gene].tolist()
    cell_lines_with_high_z_score = [element for element in column_list if element >= 6 or element <= -6]
    if len(cell_lines_with_high_z_score) > 0:
        all_genes_to_test_over.append(gene)
        genes_to_run_df_dict["Gene"].append(gene)

genes_to_run_df = pd.DataFrame.from_dict(genes_to_run_df_dict)
genes_to_run_df.to_csv("achilles_data/to_run_genes_feature_selection.csv", header=True, index=False)
#-------------------------------- END GENE EFFECT MATRIX -------------------------------------------------------------------------------------------------------------------#

#--------------------------------- EXPRESSION DATA -------------------------------------------------------------------------------------------------------------------#
# Read in CCLE Expression data and subset to cell lines in CERES file. Rename expression feature columns as well.
CCLE_EXPRESSION_MATRIX = pd.read_csv(ccle_expression_file, header=0, index_col=0)
CCLE_EXPRESSION_MATRIX = CCLE_EXPRESSION_MATRIX.T
CCLE_EXPRESSION_MATRIX = CCLE_EXPRESSION_MATRIX[all_cell_lines_to_include]
CCLE_EXPRESSION_MATRIX = CCLE_EXPRESSION_MATRIX.T

# The following code excludes all genes whos stdev is in the lower 10th percentile, as describe in the TCRP paper
expression_stdevs = CCLE_EXPRESSION_MATRIX.std()
tenth_percentile_expression_stdev = np.percentile(expression_stdevs.tolist(), 10)

to_keep_gene_expression_genes = list()
for gene in expression_stdevs.index:
    if expression_stdevs[gene] > tenth_percentile_expression_stdev:
        to_keep_gene_expression_genes.append(gene)

CCLE_EXPRESSION_MATRIX = CCLE_EXPRESSION_MATRIX[to_keep_gene_expression_genes]
# Section end

# This following code maps HGNC gene symbols, which are what is used as gene names in the CEREs data,
# to depmap gene IDs, which are gene names used in the CCLE matricies, and vice versa
hgnc_id_to_depmap_id = dict()
depmap_id_to_hgnc_id = dict()
for depmap_id in CCLE_EXPRESSION_MATRIX.columns:
    hgnc_id = depmap_id.split(' ')[0]
    if hgnc_id in hgnc_id_to_depmap_id or depmap_id in depmap_id_to_hgnc_id:
        print("ERROR!")
        quit()
    else:
        hgnc_id_to_depmap_id[hgnc_id] = depmap_id
        depmap_id_to_hgnc_id[depmap_id] = hgnc_id
# Section end

# This next section finds gene expression features that correlate with each gene KO we want to predict
# Co expression is defined in the TCRP paper as > 0.4 or < 0.4 correlation with the gene KO's expression
gene_to_coexpression_features = collections.defaultdict(set)
all_genes_to_test_over = [element for element in all_genes_to_test_over if element in hgnc_id_to_depmap_id]
for gene in all_genes_to_test_over:
    depmap_id = hgnc_id_to_depmap_id[gene]
    gene_expression = CCLE_EXPRESSION_MATRIX[depmap_id]
    for column in CCLE_EXPRESSION_MATRIX.columns:
        if column != depmap_id:
            partner_expression = CCLE_EXPRESSION_MATRIX[column]
            corr = pearsonr(gene_expression.tolist(), partner_expression.tolist())[0]
            if abs(corr) > 0.4:
                gene_to_coexpression_features[gene].add(column)
# Section end

# convert nested defaultdict to regular dict for pickling
def default_to_regular(d):
    if isinstance(d, collections.defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d

gene_to_coexpression_features = default_to_regular(gene_to_coexpression_features)

# output pickled file of all gene expression features to use based on co expression
with open('achilles_data/gene_to_coexpression_features.pkl', 'wb') as handle:
    pickle.dump(gene_to_coexpression_features, handle)

# Rename columns in ccle expression matrix so differentiate them from mutation features
CCLE_EXPRESSION_MATRIX.columns = [element + " EXPRESSION" for element in CCLE_EXPRESSION_MATRIX.columns] 
#--------------------------------- EXPRESSION END -------------------------------------------------------------------------------------------------------------------#


#--------------------------------- MUTATAION DATA -------------------------------------------------------------------------------------------------------------------#
# Read in CCLE Mutation data and subset to rows that contain info for cell lines in CERES and also exclude Silent mutations
CCLE_MUTATION_DATA = pd.read_csv(ccle_mutation_file, header=0, index_col=0, dtype=str)
CCLE_MUTATION_DATA = CCLE_MUTATION_DATA.loc[CCLE_MUTATION_DATA['DepMap_ID'].isin(all_cell_lines_to_include)]
CCLE_MUTATION_DATA = CCLE_MUTATION_DATA.loc[CCLE_MUTATION_DATA['Variant_Classification'] != "Silent"]

# One-hot encode mutation data as described in the TCRP paper
CCLE_MUTATION_MATRIX_DICT = collections.defaultdict(dict)
all_genes_with_mutations = set()
for index, row in CCLE_MUTATION_DATA.iterrows():
    depmap_id = row['DepMap_ID']
    gene_symbol = row['Hugo_Symbol']
    entrez_id = row['Entrez_Gene_Id']
    depmap_gene_string = f"{gene_symbol} ({entrez_id})"
    all_genes_with_mutations.add(depmap_gene_string)
    if depmap_gene_string not in CCLE_MUTATION_MATRIX_DICT[depmap_id]:
        CCLE_MUTATION_MATRIX_DICT[depmap_id][depmap_gene_string] = 1.0

# Process data to be fed into pandas to create a dataframe
CCLE_MUTATION_MATRIX_DF_DICT = collections.defaultdict(list)
all_genes_with_mutations = sorted(list(all_genes_with_mutations))
for cell_line in all_cell_lines_to_include:
    CCLE_MUTATION_MATRIX_DF_DICT['Cell Line'].append(cell_line)
    for gene in all_genes_with_mutations:
        if gene in CCLE_MUTATION_MATRIX_DICT[cell_line]:
            CCLE_MUTATION_MATRIX_DF_DICT[gene + " MUTATAION"].append(CCLE_MUTATION_MATRIX_DICT[cell_line][gene])
        else:
            CCLE_MUTATION_MATRIX_DF_DICT[gene + " MUTATAION"].append(0.0)

# Create mutation csv feature matrix file
CCLE_MUTATION_MATRIX = pd.DataFrame.from_dict(CCLE_MUTATION_MATRIX_DF_DICT)
CCLE_MUTATION_MATRIX = CCLE_MUTATION_MATRIX.set_index('Cell Line')

# This next section includes gene mutation features that have more than 10 counts for all cell lines,
# as described in the TCRP paper
mutation_sums = CCLE_MUTATION_MATRIX.sum()

to_keep_gene_mutation_genes = list()
for gene in mutation_sums.index:
    if mutation_sums[gene] >= 10:
        to_keep_gene_mutation_genes.append(gene)

CCLE_MUTATION_MATRIX = CCLE_MUTATION_MATRIX[to_keep_gene_mutation_genes]
# Section end
#--------------------------------- MUTATAION DATA END -------------------------------------------------------------------------------------------------------------------#

# Merge mutation and expression data into one feature matrix
FEATURE_MATRIX = CCLE_EXPRESSION_MATRIX.merge(CCLE_MUTATION_MATRIX, left_index=True, right_index=True)

# Output data, namely, feature matrix, labels matrix, and the dictionary mapping 9 lineages to few shot over and their cell lines
FEATURE_MATRIX.to_csv("achilles_data/feature_matrix.csv", header=True, index=True)
GENE_EFFECT_MATRIX.to_csv("achilles_data/labels_matrix.csv", header=True, index=True)

with open('achilles_data/lineages_to_fewshot_over.pkl', 'wb') as handle:
    pickle.dump(lineages_to_few_shot_over, handle)
    