import pickle,os,csv

### 1. Enter the name of the analysis
### 2. Enter the name of the drug. This can be a drug bank ID or the name of drug that isn't in drug bank
### 3. (Optional) Enter the targets. This is only optional if you have entered a drug from DrugBank. For investigational drugs, you will need to enter the targets associated with that drug.


analysis_name = 'PathFX_testV2'
drug_name = 'Metformin'
drug_targets = 'GENE1'

### call the algorithm without phenotype clustering
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s -t %s'%(drug_name,analysis_name,drug_targets)
print(cmd)
os.system(cmd)

#drug_name = 'FakeDrug'
#drug_targets = 'CDK4,CDK7'
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s -t %s'%(drug_name,analysis_name,drug_targets)
#os.system(cmd)

### call the algorithm with phenotype clustering
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s -c %s'%(drug_name,analysis_name,'True')
#os.system(cmd)
