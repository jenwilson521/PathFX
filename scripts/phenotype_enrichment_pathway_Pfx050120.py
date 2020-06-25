### PathFX version created: Pfx050120

# IMPORT STATEMENTS
import pickle,csv,os
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from datetime import datetime

# DATA UPDATE-SPECIFIC VARIABLES
new_hash_sum_files = pickle.load(open('../rscs/pfx041520_0.82_spec_nbhd_hash.pkl','rb'))
new_hash = pickle.load(open('../rscs/pfx041520_0.8_node_to_hashID.pkl','rb'))
dtd = pickle.load(open('../rscs/pfxDB050620_dint.pkl','rb'))
unn = pickle.load(open('../rscs/pfx041520_unique_nodes.pkl','rb'))
db = pickle.load(open('../rscs/pfxDB050620_dbid2name.pkl','rb')) # Drugbank ID to name


# BASE METHODS
rdb = dict([(v.lower(),k) for (k,v) in db.items()])

def write_neighborhood_to_file(pth_dic,outf):
	all_paths = defaultdict(float) 
	for (pth,pscore) in pth_dic.items():
		if '@' in pth:
			a=pth.split('@')[-2]
			b=pth.split('@')[-1]
			if pscore > all_paths[(a,b)]:
				all_paths[(a,b)] = pscore
		else:
			a=pth
			all_paths[(a,'')] = pscore
	for ((a,b),pscore) in all_paths.items():
		outf.write('\t'.join([a,b,str(pscore),'\n']))
	outf.close()

def create_visualization_files(dname,tlist,phen_file,net_file,outdir):
	print('creating files for visualization')
	# initialize a new network file, with phenotypes
	fN_name = os.path.join(outdir,dname + '_merged_neighborhood__withDrugTargsAndPhens.txt')
	fN  = open(fN_name,'w') # full network file
	hed = ['node1','node2','edge_score','\n'] # add header row
	fN.write('\t'.join(hed))
	for t in tlist:
		fN.write('\t'.join([dname,t,'1.0','\n']))
	net_lines = [l.strip() for l in open(net_file,'rU').readlines()]
	net_lines = [l for l in net_lines if len(l) > 2] # remove self-loops to drug-target proteins
	n = fN.write('\n'.join(net_lines)+'\n')

	# initialize a node-type file
	nT_name = os.path.join(outdir,dname +'_network_nodeType.txt')
	nT = open(nT_name,'w') # node type file
	hed = ['node_name','node_type','\n']
	nT.write('\t'.join(hed))
	nT.write('\t'.join([dname,'drug','\n']))
	for t in tlist:
		nT.write('\t'.join([t,'drug_target','\n']))
	phen_read = csv.DictReader(open(phen_file,'rU'),delimiter='\t')
	for row in phen_read:
		ph = row['phenotype']
		nT.write('\t'.join([ph,'phenotype','\n']))

		gene_list = row['genes']
		if ',' in gene_list:
			for g in gene_list.split(','):
				fN.write('\t'.join([ph,g,'1.0','\n'])) # multiple genes associated
		else:
			fN.write('\t'.join([ph,gene_list,'1.0','\n'])) # single genes associated
	nT.close()
	fN.close()


def do_network(tlist,aname,outdir,dname,doCluster):
	all_dics = []
	for t in tlist:
		# get node hash number
		nfname = new_hash[t]
		spn_name = "spn" + nfname
		tfile = new_hash_sum_files[spn_name]
			
		if os.path.exists(tfile):
			spec_dic = pickle.load(open(tfile,'rb')) 
			outf = open(os.path.join(outdir,'_'.join([t,'specific','neighborhood','.txt'])),'w')
			write_neighborhood_to_file(spec_dic,outf) # save the specific network
			all_dics+= spec_dic.items() 
		else:
			print('target not analyzed, not in interactome: '+t)
	print('merging')
	mergf = os.path.join(outdir,'_'.join([dname,'merged','neighborhood','.txt']))
	outf = open(mergf,'w')
	write_neighborhood_to_file(dict(all_dics),outf) # save the merged neighborhood
	froot = os.path.split(mergf)[-1]

	print('starting phenotypic enrichment')
	cmd = 'python ../scripts/get_network_associations_Pfx050120.py -f %s -a %s -d %s -n %s' % (froot,dname,outdir,len(tlist))
	os.system(cmd)
	
	phen_file = [f for f in os.listdir(outdir) if 'merged_neighborhood__assoc_table_.txt' in f][0] # find the association file
	net_file = [f for f in os.listdir(outdir) if 'merged_neighborhood_.txt' in f][0] # find the protein-protein network
	ph_fpath = os.path.join(outdir,phen_file)
	net_fpath = os.path.join(outdir,net_file)
	create_visualization_files(dname,tlist,ph_fpath,net_fpath,outdir)

	if doCluster:
		print('starting semantic similarity')
		fname = [f for f in os.listdir(outdir) if 'cui_list_.txt' in f][0] # find the cui list
		# cmd = 'python ../scripts/calc_lin_matrix.py -cui_list %s -a %s -d %s' % (fname,dname,outdir)	
		cmd = 'python ../scripts/calc_lin_matrix_umls.py -cui_list %s -a %s -d %s' % (fname,dname,outdir)	
		os.system(cmd)

		fname = [f for f in os.listdir(outdir) if 'lin_pandas_matrix' in f][0] # find the matrix object
		cmd = 'python ../scripts/plot_and_cluster_phenotypes.py -f %s -a %s -d %s' % (fname,dname,outdir)
		os.system(cmd)
		print('plotted phenotype clustering')

	
def main():
	parser=OptionParser()	
	parser.add_option('-d','--drug_name',dest='dname',help='Drugname, either marketed or DrugBankID')
	parser.add_option('-a','--analysis_name',dest='aname',help='Name for the Analysisi, this will also be the dirname for the results')
	parser.add_option('-i','--interactome_name',dest='iname',default="",help='The name of the interactome, default is the tissue non-specific interactome')
	parser.add_option('-t','--targets',dest='dtargs',default="",help='A comma-separated list of drug targets')
	parser.add_option('-c','--cluster',dest='docluster',default=False,help='Optional parameter to do phenotype clustering. Default = False')

	(options,args) = parser.parse_args()
	if options.dtargs:
		print(options.dtargs)
		targets = [t for t in options.dtargs.split(',')]
		print('user provided '+str(len(targets))+' targets')
	else:
		targets = []

	print('\n\nRUNNING PATHFX\n\n')
	dname = options.dname
	if dname in dtd:
		targets+= dtd[dname]
		print('PathFX added '+str(len(targets))+' targets from Drugbank')
		# check for spaces
		dname = dname.replace(' ','')
	elif dname.capitalize() in dtd:
		targets+= dtd[dname.capitalize()]
		print('PathFX added '+str(len(targets))+' targets from Drugbank')
		# check for spaces
		dname = dname.capitalize().replace(' ','')
	elif dname.lower() in rdb:
		dname = rdb[dname.lower()]
		targets+= dtd[dname]	
		print('PathFX added '+str(len(targets))+' targets from Drugbank')
	elif dname.replace(' ','') in rdb:
		dname = rdb[dname.replace(' ','')]
		targets+= dtd[dname]
		print('PathFX added '+str(len(targets))+' targets from Drugbank')
	else:
		print('No drug bank targets found')
	targets = list(set(targets)) # remove redundancies

	global fpath
	global scr
	if options.iname:
		(fpath,scr) = intome_data[options.iname]
	else:
		fpath = '../results/new_intome_rands_0.77/'
	#	scr = '0.85' # use the default interatome
		scr = '0.77' # use the default interatome

	aname = options.aname.replace(' ','_').lower()
	# dtime = datetime.now().isoformat().replace(':','-')
	# aname = '_'.join([aname,dtime])
	outdir = os.path.join('../results/',aname,dname)
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	print('\n'+dname)
	print('results in: '+outdir)
	print('full target list: '+' '.join(targets))

	# check which targets are in the network
	t_in_n = list(set(targets).intersection(set(unn)))
	if len(t_in_n) > 0:
		print('for analysis (targets within the interactome): '+' '.join(t_in_n))
	else:
		print('no targets provided or none were found in the interaction network')

	doCluster = options.docluster

	do_network(t_in_n,aname,outdir,dname,doCluster)

if __name__ == "__main__":
	main()

