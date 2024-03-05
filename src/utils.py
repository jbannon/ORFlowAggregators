import numpy as np
from scipy.stats import kurtosis,mode, skew, beta
import sys
from typing import List, Dict, Union
import networkx as nx
import pandas as pd
from collections import defaultdict
import numpy as np
from scipy.stats import kurtosis,mode, skew, beta
import sys
from typing import List, Dict, Union
import networkx as nx
import pandas as pd
from collections import defaultdict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC


DRUG_TISSUE_MAP = {"Atezo":["KIRC","BLCA"],"Pembro":["SKCM","STAD"],
    "Nivo":["KIRC","SKCM"], "Ipi":["SKCM"], "Ipi+Pembro":["SKCM"],
    "erlotinib":['LUAD'],"crizotinib":['LUAD'],'sorafenib':["LUAD"],'sunitinib':["LUAD"]}


DRUG_DATASET_MAP = {
    'sorafenib':'ccle',
    'erlotinib':'ccle',
    'crizotinib':'ccle',
    'sunitinib':'ccle',
    'Nivo':'cri',
    'Ipi':'cri',
    'Pembro':'cri',
    'Atezo':'cri',
    'Ipi+Pembro':'cri'
}

DATASET_FEAT_NAME_MAP = {'ccle':"metab_features.csv",'cri':"immune_features.csv"}

DRUG_TARGET_MAP = {'Atezo':'PD-L1','Pembro':'PD1','Nivo':'PD1','Ipi':'CTLA4'}

TARGET_GENE_MAP = {'PD-L1':'CD274', 'PD1':'PDCD1', 'CTLA4':'CTLA4'}


def fetch_gene_set(
    gs_name:str
    ):
    gs_path = "../data/genesets/{gs}.txt".format(gs=gs_name.upper())
    with open(gs_path,"r") as istream:
        lines = istream.readlines()
    lines = [x.rstrip() for x in lines]
    return lines

def unpack_parameters(
    D:Dict
    ):
    if len(D.values())>1:
        return tuple(D.values())
    else:
        return tuple(D.values())[0]

def make_file_path(
    base_dir:str, 
    pathNames:List[str], 
    fname:str, 
    ext:str
    )->str:
    
    pathNames.append("")
    ext = ext if ext[0]=="." else "." + ext

    base_dir = base_dir[:-1] if base_dir[-1]=="/" else base_dir
    path = [base_dir]
    path.extend(pathNames)
    path ="/".join([x for x in path])
    
    file_path = "".join([path,fname,ext])
 
    

    return file_path


def fetch_drug_targets(
    drug:str
    ) -> List[str]:
    if drug in DRUG_TARGET_MAP.keys():
        targets = [TARGET_GENE_MAP[DRUG_TARGET_MAP[drug]]]
    elif drug == "Ipi+Pembro":
        t1 = TARGET_GENE_MAP[DRUG_TARGET_MAP["Ipi"]]
        t2 = TARGET_GENE_MAP[DRUG_TARGET_MAP["Pembro"]]
        targets = [t1,t2]
    else:
        fname = "../data/genesets/{d}_targets.txt".format(d=drug)
        with open(fname, "r") as istream:
            lines = istream.readlines()
        targets = [x.rstrip() for x in lines]

    return targets

def empirical_bayes_gene_selection(
    df:pd.DataFrame,
    probThresh:float,
    nRuns:int = 200):
    

    temp = df.copy(deep=True)
    
    temp["Hits"] = temp['Count']
    temp['Misses'] = nRuns-temp['Count']
    temp['Proportion'] = temp["Count"]/nRuns

    a, b,loc, scale = beta.fit(temp[temp["Proportion"]<1.0]["Proportion"].values,fscale = 1, floc = 0)

    temp["EB_avg"] = (temp["Hits"]+a)/(nRuns+a+b)
            
    cut_point = np.amin(temp[temp["EB_avg"]>=probThresh]["Count"].values)
            
   
    temp = temp[temp['Count'] >= cut_point]
  
    gene_list = list(temp['Gene'].values)
    return gene_list

    
def make_model_and_param_grid(
    model_name:str,
    reg_min:float,
    reg_max:float,
    reg_step:float,
    model_max_iters:int
    ):
    
    preproc = ('preproc',StandardScaler())
    
    if model_name == 'LogisticRegression':
        classifier = ('clf',LogisticRegression(class_weight = 'balanced',max_iter = model_max_iters))
    elif model_name == 'LinearSVM':
        classifier = ('clf',LinearSVC(class_weight = 'balanced',max_iter = model_max_iters))
    

    param_grid = {'clf__C':np.arange(reg_min,reg_max,reg_step)}                            
    
    model = Pipeline([preproc,classifier])

    return model, param_grid

