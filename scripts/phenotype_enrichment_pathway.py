# revision of the analysis pipeline
# to draw on the completed randomizations instead of creating
# all neighborhoods individualls
# written 11-8-17 JLW

import pickle,csv,os
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from datetime import datetime

# drugs to targets
dtf = os.path.join('..','rscs','mapped_drugs_to_targets.pkl')
dtf = os.path.join('..','rscs','drug_intome_targets.pkl')
dtd = pickle.load(open(dtf,'rb'))
dbf = os.path.join('..','rscs','drugbankid_to_name.pkl')
db = pickle.load(open(dbf,'rb'))
rdb = dict([(v.lower(),k) for (k,v) in db.items()])
file_swap = pickle.load(open(os.path.join('..','rscs','file_swap.pkl'),'rb'))
sum_file_swap = pickle.load(open(os.path.join('..','rscs','summary_file_swap.pkl'),'rb'))

# default results
intome_data = {'cardiac':('../rsc/FIXTHIS','0.96')} # tuples of file paths and score values 
long_name = 'CYP2C19star1,CYP2C19star10,CYP2C19star11,CYP2C19star13,CYP2C19star14,CYP2C19star15,CYP2C19star16,CYP2C19star18,CYP2C19star19,CYP2C19star22,CYP2C19star23,CYP2C19star24,CYP2C19star25,CYP2C19star26,CYP2C19star5,CYP2C19star6,CYP2C19star8,CYP2C19star9'
short_name = 'CYP2C19star1,5-10,11,13-16,18-19,22-26'

def clean_node_name(sn):
	if '+' in sn:
		sn = sn.replace('+','-')
	if ':' in sn:
		sn = sn.replace(':','-')
	if '/' in sn:
		sn = sn.replace('/','-')
	if ' ' in sn:
		sn = sn.replace(' ','')
	if '*' in sn:
		sn = sn.replace('*','star')
	if long_name in sn:
		sn = sn.replace(long_name,short_name)

	return sn

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

def create_specific_dic(pth_dic):
	spec_dic = []
	for (pth,s) in pth_dic.items():
		if '@' in pth:
			pg = pth.split('@')[-1]
		else:
			spec_dic.append((pth,1)) # the target stays in the specific network
			continue
		pg = clean_node_name(pg)
		grf = os.path.join(fpath,'summary','_'.join([pg,str(scr),'randPathScores','.pkl']))
		if grf in sum_file_swap:
			grf = sum_file_swap[grf] # swap in the case of long file names
		pg_scores = pickle.load(open(grf,'rb'))
		avg = np.mean(pg_scores)
		if (s-avg) >0:
			spec_dic.append((pth,s))

	spec_dic = dict(spec_dic)
	return spec_dic

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
			for g in gene_list:
				fN.write('\t'.join([ph,g,'1.0','\n'])) # multiple genes associated
		else:
			fN.write('\t'.join([ph,gene_list,'1.0','\n'])) # single genes associated
	nT.close()
	fN.close()


def do_network(tlist,aname,outdir,dname,doCluster):
	all_dics = []
	for t in tlist:
		t = clean_node_name(t)
		tfile = os.path.join(fpath,'pth_dic_'+t+'_scr'+scr+'.pkl')
		if tfile in file_swap:
			tfile = file_swap[tfile]
		if os.path.exists(tfile):
			tdic = pickle.load(open(tfile,'rb'))
			outf = open(os.path.join(outdir,'_'.join([t,'neighborhood','.txt'])),'w')
			write_neighborhood_to_file(tdic,outf) # save individual network
			spec_dic = create_specific_dic(tdic)
			outf = open(os.path.join(outdir,'_'.join([t,'specific','neighborhood','.txt'])),'w')
			write_neighborhood_to_file(spec_dic,outf) # save the specific network
			all_dics+= spec_dic.items() 
		else:
			print('target not analzyed, not in interactome: '+t)
	print('merging')
	mergf = os.path.join(outdir,'_'.join([dname,'merged','neighborhood','.txt']))
	outf = open(mergf,'w')
	write_neighborhood_to_file(dict(all_dics),outf) # save the merged neighborhood
	froot = os.path.split(mergf)[-1]

	print('starting phenotypic enrichment')
#	cmd = 'python get_network_associations_v2.py -f %s -a %s -d %s' % (froot,dname,outdir)
	cmd = 'python ../scripts/get_network_associations_v3.py -f %s -a %s -d %s' % (froot,dname,outdir)
	os.system(cmd)
#	print(cmd)
	
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
		targets = [t for t in options.dtargs.split(',')]
	else:
		targets = []

	print('\n\nRUNNING PATHFX\n\n')
	dname = options.dname
	if dname in dtd:
		targets+= dtd[dname]
	elif dname.lower() in rdb:
		dname = rdb[dname.lower()]
		targets+= dtd[dname]	
	elif dname.replace(' ','') in rdb:
		dname = rdb[dname.replace(' ','')]
		targets+= dtd[dname]
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
	dtime = datetime.now().isoformat().replace(':','-')
	aname = '_'.join([aname,dtime])
	outdir = os.path.join('../results/',aname,dname)
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	print('\n'+dname)
	print('results in: '+outdir)
	print('targets: '+' '.join(targets))

	doCluster = options.docluster

	do_network(targets,aname,outdir,dname,doCluster)

if __name__ == "__main__":
	main()

