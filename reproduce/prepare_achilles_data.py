import pickle
import numpy as np
import h5py
import os

gene_to_features = pickle.load(open("./achilles_data/gene-to-features.pkl",'rb'))
lineages = pickle.load(open("./achilles_data/lineage-to-cell-lines.pkl",'rb'))
labels = h5py.File("./achilles_data/labels-matrix.hdf5", 'r')
features = h5py.File("./achilles_data/feature-matrix.hdf5", 'r')
genes = ["HNF1B","ESR1","KLB","FAM183A","JAK3","CHST8","KCNK13","SPI1","SALL4","WNT16"]

for gene in genes:
    os.mkdir("./data/"+gene)
    for lineage in lineages.keys():
        np.save(open("./data/"+gene+"/" + lineage + "_feature_description.npy",'wb'),gene_to_features[gene])
        lineage_labels = labels["data"][np.isin(labels["dim_0"], list(lineages[lineage])),:]
        lineage_labels = lineage_labels[:,np.isin(labels["dim_1"],[gene])]
        lineage_labels = np.reshape(lineage_labels,[np.size(lineage_labels)])
        np.save(open("./data/"+gene+"/" + lineage + "_"+gene+"_label.npy",'wb'),lineage_labels)
        lineage_features = features["data"][np.isin(features["dim_0"], list(lineages[lineage])),:]
        lineage_features = lineage_features[:,np.isin(features["dim_1"],list(gene_to_features[gene]))]
        np.save(open("./data/"+gene+"/" + lineage + "_"+gene+"_feature.npy",'wb'),lineage_features)
