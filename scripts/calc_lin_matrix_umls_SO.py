# Written to get the lin matrix from
# the UMLS database instead of the CUI Graph
# Written 7-24-18, JLW

import pickle,os,csv
import argparse
import itertools
import pandas as pd
import numpy as np


### Might need to call some commands here to start the database server ###
# cmd = "perlbrew use 5.26.0"
# os.system(cmd)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-cui_list', action='store', dest='clist', help='the first CUI')
parser.add_argument('-outf', action='store', dest='outf', help='the name for the output file')
args = parser.parse_args()

cui_file = args.clist
out_file = args.outf

cmd = "umls-similarity.pl --config itfc_confg.txt --infile %s --measure lin --matrix --realtime > %s"%(cui_file,out_file) 
print(cmd)
os.system(cmd)



