import os
import tqdm

import argparse
import yaml
from sklearn.metrics import confusion_matrix, roc_auc_score, auc, precision_recall_curve, accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
from typing import Dict, List
import sys
from sklearn.neighbors import kneighbors_graph
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, balanced_accuracy_score
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import kneighbors_graph
import utils 
from collections import defaultdict

from sklearn.metrics import confusion_matrix, roc_auc_score, auc, precision_recall_curve, accuracy_score




def main(
	config:Dict
	) -> None:
	
	drug:str
	tissues:List[str]
	thresh:float 
	num_iters: int
	
	min_train_size: float
	max_train_size: float
	train_size_step: float
	


	geneset_map = {}

	drug, tissues, genesets, num_iters, min_train_size, max_train_size, train_size_step, rng_seed = utils.unpack_parameters(config["EXPERIMENT_PARAMS"])
	model_name, reg_min, reg_max, reg_step, model_max_iters = utils.unpack_parameters(config["MODEL_PARAMS"])


	rng = np.random.RandomState(rng_seed)

	for g in genesets:
		geneset_map[g]=utils.fetch_gene_set(g)

	feature_sets = ["TARGETS","FEATURES"] + genesets + [x+"+FEATURES" for x in genesets] + ["FEATURES+TARGETS"]
	

	
	

	targets = utils.fetch_drug_targets(drug)
	source_dataset = utils.DRUG_DATASET_MAP[drug]
	
	for tissue in tissues:
		model, param_grid = utils.make_model_and_param_grid(model_name,reg_min, reg_max,reg_step, model_max_iters)

		features = pd.read_csv("../data/expression/{ds}/{d}/{t}/features.csv".format(d=drug, t=tissue,ds=source_dataset))	
		expr = pd.read_csv("../data/expression/{ds}/{d}/{t}/expression_full.csv".format(d=drug,t=tissue, ds=source_dataset))	
		resp = pd.read_csv("../data/expression/{ds}/{d}/{t}/response.csv".format(d=drug,t=tissue,ds = source_dataset))
		
		if features.shape[0] != resp.shape:
			print("fixing feature shape")
			obs_order = [x for x in expr['Run_ID'].values if x in features['Run_ID'].values]
			
			features.set_index("Run_ID",inplace = True)
			features = features.loc[obs_order,:]
			features.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			expr.set_index("Run_ID",inplace = True)
			expr = expr.loc[obs_order,:]
			expr.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			resp.set_index("Run_ID",inplace = True)
			resp = resp.loc[obs_order,:]
			resp.reset_index(inplace = True, drop = False, names = ['Run_ID'])

		
		results = defaultdict(list)
		res_dir = "../results/baselines/{d}/{t}/".format(d=drug, t=tissue)
		os.makedirs(res_dir, exist_ok = True)
		

		for feature_set in feature_sets:
			if  feature_set == 'TARGETS':
				X = np.log2(expr[targets].values+1)
				y = resp['Response'].values
			elif feature_set == "FEATURES":
				X = features[features.columns[1:]].values
			elif feature_set in genesets:
				X = np.log2(expr[targets + [x for x in geneset_map[feature_set] if x in expr.columns]].values+1)
			else:
				fs, mod = feature_set.split("+")

				if fs == "FEATURES" and mod =="TARGETS":
					X1 = features[features.columns[1:]].values
					X2 = expr[targets].values
					X = np.hstack((X1,X2))
				elif fs in genesets and mod == "FEATURES":
					X1 = np.log2(expr[targets + [x for x in geneset_map[fs] if x in expr.columns]].values+1)
					X2 = features[features.columns[1:]].values
					X = np.hstack((X1,X2))
			y = resp['Response'].values
			
		
			for train_pct in tqdm.tqdm(np.round(np.arange(min_train_size,max_train_size,train_size_step),2)):
				clf = GridSearchCV(model, param_grid)
				
				for i in range(num_iters):
					X_train, X_test, y_train, y_test = train_test_split(X, y,train_size = train_pct, 
						random_state = rng,shuffle = True, stratify = y)
					clf.fit(X_train,y_train)

					train_preds_bin = clf.predict(X_train)
					train_preds_prob = clf.predict_proba(X_train)

					test_preds_bin = clf.predict(X_test)
					test_preds_prob = clf.predict_proba(X_test)


					train_acc = accuracy_score(y_train, train_preds_bin)
					test_acc = accuracy_score(y_test, test_preds_bin) 

					train_roc = roc_auc_score(y_train,train_preds_prob[:,1])
					test_roc = roc_auc_score(y_test,test_preds_prob[:,1])

					tn, fp, fn, tp = confusion_matrix(y_test, test_preds_bin,labels = [0,1]).ravel()
					

					results['drug'].append(drug)
					results['tissue'].append(tissue)
					results['Feature Set'].append(feature_set)
					results['Train Percentage'].append(train_pct)
					results['Train Accuracy'].append(train_acc)
					results['Train ROC_AUC'].append(train_roc)
					results['Test Accuracy'].append(test_acc)
					results['Test ROC_AUC'].append(test_roc)
					results['Test TN'].append(tn)
					results['Test FP'].append(fp)
					results['Test FN'].append(fn)
					results['Test TP'].append(tp)

		results = pd.DataFrame(results)
		res_name = res_dir + "LR_baselines.csv"
		results.to_csv(res_name, index = False)
		



	






if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()
	
	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
	

