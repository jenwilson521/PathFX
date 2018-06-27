import pickle,os,csv
from datetime import datetime

startTime = datetime.now()
### 1. Enter the name of the analysis
### 2. Enter the name of the drug. This can be a drug bank ID or the name of drug that isn't in drug bank
### 3. (Optional) Enter the targets. This only works if you have entered a drug from DrugBank. For investigational drugs, you will need to enter the targets associated with that drug.

analysis_name = 'PathFX_mike_demo'

drug_name = 'metformin'
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
os.system(cmd)
print(datetime.now() - startTime)

startTime = datetime.now()
drug_name = 'Montelukast'
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
os.system(cmd)
print(datetime.now() - startTime)

#drug_name = 'pravastatin'
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
#os.system(cmd)
#
#
#drug_name = 'propofol'
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
#os.system(cmd)



