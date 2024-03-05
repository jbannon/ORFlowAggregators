from typing import Dict, List
import sys
from sklearn.neighbors import kneighbors_graph
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import seaborn as sns
import os



def main(
	config:Dict
	) -> None:
	
	drug:str = "Atezo"
	tissue:str = "BLCA"
	thresh:float = 0.01
	cutoff: int = 150
	drug_tissue_map = {'Atezo':['BLCA','KIRC'],
	'Ipi':['SKCM'],
	'Ipi+Pembro':['SKCM'],
	'Nivo':['KIRC','SKCM'],
	'Pembro':['SKCM','STAD']}

	for drug in ['Atezo', 'Nivo','Ipi','Ipi+Pembro','Pembro']:
		for tissue in drug_tissue_map[drug]:
			outdir = "../figs/exploratory/{d}/{t}/".format(d=drug,t=tissue)
			os.makedirs(outdir, exist_ok = True)
			resp = pd.read_csv("../data/{d}/{t}/response.csv".format(d=drug,t=tissue))
			feat = pd.read_csv("../data/{d}/{t}/immune_features.csv".format(d=drug,t=tissue))
			data = resp.merge(feat,on='Run_ID')
			for feature in feat.columns[1:]:
				if 'Flag' not in feature:
					print(feature)
					sns.boxplot(x='Response',y=feature,data = data)
					plt.savefig("{od}{f}".format(od=outdir,f=feature))
					plt.close()
	
	





if __name__ == '__main__':
	main({})

