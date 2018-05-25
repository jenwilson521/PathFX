# written to return the lin distance using pre-calculated
# sanchez IC
# written 50-24-18 JLW

import pickle,os,csv
import cui_query as cq
import networkx as nx
import argparse
import itertools
import pandas as pd
import numpy as np

# load the network
# cq.combG = nx.read_gpickle("../rscs/combinedCuiOntoGraph.gpickle")
cq.combG = nx.read_gpickle("../rscs/cuiGraph.gpickle")
cq.all_nodes = len(cq.combG.nodes())
cq.all_cui_nodes = len([k for k in cq.combG.nodes() if cq.combG.node[k]["type"] == "CuiNode"])


#ic_dir = '../results/ic_values/'
ic_dir = '../results/ic_values_v2/'
allf = [f for f in os.listdir(ic_dir)]
ic_dic = dict([(f.split('_')[0],os.path.join(ic_dir,f)) for f in allf])

def get_ic(cui):
	if cui in ic_dic:
		ic = pickle.load(open(ic_dic[cui],'rb'))
		return ic
	elif cui not in cq.combG:
#		print('CUI not in graph: '+cui)
		return 0.
	else:
#		print('CUI IC not computed: '+cui)
		return 0.

def lin_sim(iri1, iri2, only_cui=True): #, method=sanchez_ic):
	if iri1 not in cq.combG or iri2 not in cq.combG: # skip if both are not in the graph
		sim = 0
	else: # try computing similarity
		lcs = cq.get_lcs(iri1, iri2, only_cui)
		ic_lcs = get_ic(lcs) # method(lcs, only_cui)
		ic_c1 = get_ic(iri1) # method(iri1, only_cui)
		ic_c2 = get_ic(iri2) # method(iri2, only_cui)
		if ic_c1 >0. and ic_c2 > 0.:    
			sim = 2*ic_lcs/(ic_c1+ic_c2)
		else: # don't divide by zero
			sim = 0.
	return sim


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-cui_list', action='store', dest='clist', help='the first CUI')
parser.add_argument('-a','--analysis_name',dest='aname',help='Name of analysis, will be appended to output files; experiment date is suggested')
parser.add_argument('-d','--dir',dest='res_dir',help='Results directory. If none provided, a directory will be created matching the analysis name in the ../results/ dir')
args = parser.parse_args()

if args.res_dir is None:
	rdir = os.path.join('..','results',args.aname)
else:   
	rdir = args.res_dir

print('Loading the CUI file')
cui_file = os.path.join(rdir,args.clist)
cui_list = sorted(pickle.load(open(cui_file,'rb')))
cui_ind = dict([(cui,i) for (i,cui) in enumerate(cui_list)])
num_cui = len(cui_list)

sem_sim = np.zeros((num_cui,num_cui))
for c1,c2 in itertools.combinations(cui_list,2):
	lin = lin_sim(c1, c2, only_cui=True)
	if lin > 0.:
		c1_ind = cui_ind[c1]
		c2_ind = cui_ind[c2]
		sem_sim[c1_ind,c2_ind] = lin

# fill the diagonal
print('Calculating pairwise similarity')
for c in cui_list:
	c_ind = cui_ind[c]
	sem_sim[c_ind,c_ind] = 1.

sem_res = pd.DataFrame(sem_sim,index=cui_list,columns=cui_list)
outf = os.path.join(rdir,'lin_pandas_matrix.pkl')
pickle.dump(sem_res,open(outf,'wb'))

