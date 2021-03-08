# This code is all CDS code
# Author: David Wu (dwu@broadinstitute.org)

import collections
import pickle
from taigapy import TaigaClient
from tqdm import *

# Read in CDS team's custom data file mapping genes to their ppi partners and paralogs
tc = TaigaClient()
related_features = tc.get(name='input-data-0181', version=5, file='MatchRelated')

# Populate dictionary
gene_to_related_features = collections.defaultdict(set)
for index, row in tqdm(related_features.iterrows(), total=related_features.shape[0]):
    target = row['target']
    partner = row['partner']
    gene_to_related_features[target].add(partner)

# Convert nested defaultdict to regular dict to facilitate pickling
def default_to_regular(d):
    if isinstance(d, collections.defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d

gene_to_related_features = default_to_regular(gene_to_related_features)

# Dump dictionary to pickle file
with open('achilles_data/ppi_features.pkl', 'wb') as handle:
    pickle.dump(gene_to_related_features, handle)
    