import pickle,os,csv
import numpy as np

### 1. Enter the name of the analysis
### 2. Enter the name of the drug. This can be a drug bank ID or the name of drug that isn't in drug bank
### 3. (Optional) Enter the targets. This only works if you have entered a drug from DrugBank. For investigational drugs, you will need to enter the targets associated with that drug.


analysis_name = 'PathFX_vis_test'
dbid2namef = '../rscs/drugbankid_to_name.pkl'
dbid2n = pickle.load(open(dbid2namef,'rb'))
dbids = [k for k in dbid2n.keys()]

for drug_name in ['DB03173','DB00592','DB03835']:
	cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
	os.system(cmd)
