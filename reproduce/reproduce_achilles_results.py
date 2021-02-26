import time
import argparse
import numpy as np
import random
import torch
import torch.nn.functional as F
import torch.optim as optim
import os
import glob
from torch.autograd import Variable
import sys
import torch.nn as nn
import pickle
import copy
import pandas as pd
import h5py

sys.path.append('../code/')
from data_loading import *
from utils import *
from score import *
from inner_loop import InnerLoop
from mlp import mlp
from meta_learner_cv import *

# Training settings
parser = argparse.ArgumentParser()
parser.add_argument('--fix_lineage_selection_issue', type=bool, default=False, help='Whether to fix the lineage selection issue or not')
parser.add_argument('--min_lines_per_lineage', type=int, default=15, help='Minimun number of cell lines per tissue')
parser.add_argument('--fix_train_set_issue', type=bool, default=False, help='Whether to fix the train issue or not')
parser.add_argument('--trials_for_each_K', type=int, default=1, help='number of trial for each K value')
parser.add_argument('--prepare_data', type=bool, default=True, help='Whether prepare feature and label data or not')
parser.add_argument('--save_results', type=str, default="achilles_results.csv", help='file to save results to')
parser.add_argument('--genes', nargs='+', default=[], help='list of genes to test')
args = parser.parse_args()

#Use default arguments from tcrp_cv.py
#parser.add_argument('--meta_batch_size', type=int, default=32, help='Meta-learning batch size, i.e. how many different tasks we need sample')
meta_batch_size = 32
#parser.add_argument('--inner_batch_size', type=int, default=10, help='Batch size for each individual learning job')
inner_batch_size = 32
#parser.add_argument('--num_updates', type=int, default=30, help='Number of training epochs')
num_updates = 15 #set back to 30
#parser.add_argument('--num_inner_updates', type=int, default=1, help='Initial learning rate')
num_inner_updates = 1
#parser.add_argument('--num_trials', type=int, default=50, help='Number of trials for unseen tissue')
num_trials = 50
#parser.add_argument('--hidden', type=int, default=60, help='Number of hidden units of NN for single task')
hidden = 60
#parser.add_argument('--patience', type=int, default=3, help='Patience')
patience = 3
#parser.add_argument('--meta_lr', type=float, default=0.001, help='Learning rate for meta-learning update')
meta_lr = 0.001
#parser.add_argument('--inner_lr', type=float, default=0.001, help='Learning rate for ')
inner_lr = 0.001
#parser.add_argument('--tissue_num', type=int, default=12, help='Tissue number evolved in the inner update')
tissue_num = 12
#parser.add_argument('--layer', type=int, default=1, help='Number of layers of NN for single task')
layer = 1

# Prepare data
data_path = './data/'
tissue_map = pickle.load(open('./achilles_data/lineage-to-cell-lines.pkl', 'rb'))

if args.prepare_data:
	gene_to_features = pickle.load(open("./achilles_data/gene-to-features.pkl",'rb'))
	labels = h5py.File("./achilles_data/labels-matrix.hdf5", 'r')
	features = h5py.File("./achilles_data/feature-matrix.hdf5", 'r')

	for gene in args.genes:
	    os.mkdir(data_path+gene)
	    for lineage in tissue_map.keys():
	        np.save(open(data_path+gene+"/" + lineage + "_feature_description.npy",'wb'),gene_to_features[gene])
	        lineage_labels = labels["data"][np.isin(labels["dim_0"], list(tissue_map[lineage])),:]
	        lineage_labels = lineage_labels[:,np.isin(labels["dim_1"],[gene])]
	        lineage_labels = np.reshape(lineage_labels,[np.size(lineage_labels)])
	        np.save(open(data_path+gene+"/" + lineage + "_"+gene+"_label.npy",'wb'),lineage_labels)
	        lineage_features = features["data"][np.isin(features["dim_0"], list(tissue_map[lineage])),:]
	        lineage_features = lineage_features[:,np.isin(features["dim_1"],list(gene_to_features[gene]))]
	        np.save(open(data_path+gene+"/" + lineage + "_"+gene+"_feature.npy",'wb'),lineage_features)

results_gene_list = []
results_trial_list = []
results_K_list = []
results_train_cor_mean_list = []
results_test_cor_mean_list = []

for drug in args.genes:
	print(drug)
	for K in range(1,10,2):
		for t in range(args.trials_for_each_K):

			if(args.fix_lineage_selection_issue):
				# Hold out lineages with at least min_lines_per_lineage cell lines in them
				cv_feature_list, cv_label_list, meta_tissue_index_list, test_feature_list, test_label_list, test_tissue_list  = load_data_cell_line(tissue_map, drug, args.min_lines_per_lineage, data_path)
			else:
				# Hold out lineages with at least K cell lines in them. This is invalid because it means different K values are not comparable.
				cv_feature_list, cv_label_list, meta_tissue_index_list, test_feature_list, test_label_list, test_tissue_list  = load_data_cell_line(tissue_map, drug, K, data_path)

			# Following Code from tcrp_cv.py
			best_fewshot_train_corr_list = []
			best_fewshot_test_corr_list = []

			for i, test_tissue in enumerate(test_tissue_list):

				#print 'Validate tissue', test_tissue, cv_feature_list[i].shape

				meta_dataset = dataset(cv_feature_list[i], cv_label_list[i])
				test_dataset = dataset(test_feature_list[i], test_label_list[i])

				meta_learner = MetaLearner( meta_dataset, test_dataset, K, meta_lr, inner_lr, layer, hidden, tissue_num, meta_batch_size, inner_batch_size, num_updates, num_inner_updates, meta_tissue_index_list[i], num_trials )

				best_fewshot_train_loss, best_fewshot_train_corr, best_fewshot_test_loss, best_fewshot_test_corr, best_train_corr_model = meta_learner.train(args.fix_train_set_issue)

				print test_tissue, 'best few-shot train loss', best_fewshot_train_loss, 'best few-shot train corr', best_fewshot_train_corr
				print 'best few-shot test loss', best_fewshot_test_loss, 'best few-shot test corr', best_fewshot_test_corr

				if best_fewshot_train_corr != -1:
					best_fewshot_train_corr_list.append(best_fewshot_train_corr)

				if best_fewshot_test_corr != -1:
					best_fewshot_test_corr_list.append( best_fewshot_test_corr )

			print 'Avg_test corr', np.asarray(best_fewshot_train_corr_list).mean(), np.asarray(best_fewshot_test_corr_list).mean()

			results_gene_list.append(drug)
			results_trial_list.append(t)
			results_K_list.append(K)
			results_train_cor_mean_list.append(np.nanmean(np.asarray(best_fewshot_train_corr_list)))
			results_test_cor_mean_list.append(np.nanmean(np.asarray(best_fewshot_test_corr_list)))

df = pd.DataFrame(list(zip(results_gene_list, results_trial_list, results_K_list, results_train_cor_mean_list, results_test_cor_mean_list)),
			   columns =['gene', 'trial', 'K', 'train_cor_mean', 'test_cor_mean'])

df.to_csv(open(args.save_results,'w'))
