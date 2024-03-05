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

	feature_sets = ["TARGETS","FEATURES","FEATURES+TARGETS"]
	

	
	

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
		
		os.makedirs(res_dir, exist_ok = True)

		g = kneighbors_graph(expression[de_genes].values,n_neighbors=5)
		G = nx.from_scipy_sparse_array(g)
		print(nx.is_connected(G))
		A = g.toarray()

		for train_pct in tqdm.tqdm()
		
		
			
	
	






if __name__ == '__main__':
	main({})

