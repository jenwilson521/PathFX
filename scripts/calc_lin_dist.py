# written to return the lin distance using pre-calculated
# sanchez IC
# written 50-24-18 JLW

import pickle,os,csv
import cui_query as cq
import networkx as nx
import argparse

# load the network
cq.combG = nx.read_gpickle("combinedCuiOntoGraph.gpickle")
cq.all_nodes = len(cq.combG.nodes())
cq.all_cui_nodes = len([k for k in cq.combG.nodes() if cq.combG.node[k]["type"] == "CuiNode"])


ic_dir = 'ic_values/'
allf = [f for f in os.listdir(ic_dir)]
ic_dic = dict([(f.split('_')[0],os.path.join(ic_dir,f)) for f in allf])

def get_ic(cui):
	if cui in ic_dic:
		ic = pickle.load(open(ic_dic[cui],'rb'))
		return ic
	elif cui not in cq.combG:
		print('CUI not in graph: '+cui)
		return 0.
	else:
		print('CUI IC not computed: '+cui)
		return 0.

def lin_sim(iri1, iri2, only_cui=True): #, method=sanchez_ic):
	lcs = cq.get_lcs(iri1, iri2, only_cui)
	ic_lcs = get_ic(lcs) # method(lcs, only_cui)
	ic_c1 = get_ic(iri1) # method(iri1, only_cui)
	ic_c2 = get_ic(iri2) # method(iri2, only_cui)
	if ic_c1 >0. and ic_c2 > 0.:    
		sim = 2*ic_lcs/(ic_c1+ic_c2)
	else:
		sim = 0.
	return sim

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-c1', action='store', dest='c1', help='the first CUI')
parser.add_argument('-c2', action='store', dest='c2', help='the second CUI')
args = parser.parse_args()
[c1,c2] = [args.c1,args.c2]

print(lin_sim(c1, c2, only_cui=True))

