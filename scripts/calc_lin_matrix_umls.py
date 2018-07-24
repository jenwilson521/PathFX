# Written to get the lin matrix from
# the UMLS database instead of the CUI Graph
# Written 7-24-18, JLW

import pickle,os,csv
import argparse
import itertools
import pandas as pd
import numpy as np


### Might need to call some commands here to start the database server? ###
# cmd = "perlbrew use 5.26.0"
# os.system(cmd)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-cui_list', action='store', dest='clist', help='the first CUI')
parser.add_argument('-a','--analysis_name',dest='aname',help='Name of analysis, will be appended to output files; experiment date is suggested')
parser.add_argument('-d','--dir',dest='res_dir',help='Results directory. If none provided, a directory will be created matching the analysis name in the ../results/ dir')
args = parser.parse_args()

if args.res_dir is None:
	rdir = os.path.join('..','results',args.aname)
else:   
	rdir = args.res_dir

cui_file = os.path.join(rdir,args.clist)
out_file = os.path.join(rdir,'lin_pandas_matrix.txt')
cmd = "umls-similarity.pl --config itfc_confg.txt --infile %s --measure lin --matrix --realtime > %s"%(cui_file,out_file) 
os.system(cmd)



