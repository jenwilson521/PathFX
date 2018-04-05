import pickle,os,csv

### 1. Enter the name of the analysis
### 2. Enter the name of the drug. This can be a drug bank ID or the name of drug that isn't in drug bank
### 3. (Optional) Enter the targets. This only works if you have entered a drug from DrugBank. For investigational drugs, you will need to enter the targets associated with that drug.


analysis_name = 'PathFX_demo'
drug_name = 'Jenimab'
drug_targets = 'EGFR,ABL1,ABCB1'

### call the algorithm
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s -t %s'%(drug_name,analysis_name,drug_targets)
os.system(cmd)
