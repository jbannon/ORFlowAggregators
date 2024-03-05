from typing import Dict, Tuple
import matplotlib.pyplot as plt
import sys
import numpy as np
import networkx as nx
import ot
import tqdm 
import time
import pandas as pd 



df = pd.read_csv("../data/expression/cri/Atezo/PANCAN/expression.csv")
de = pd.read_csv("../data/genesets/Atezo_PANCAN_DE.csv")

de = de[de['Thresh.Value']==0.01]
de = de[de['Count']>150]
print(de)

X = np.log2(df[df.columns[1:]].values+1)
v = np.var(X,axis=0)
print(v.shape)
print(np.argsort(-v))
print(v[np.argsort(-v)])
print(X.shape)