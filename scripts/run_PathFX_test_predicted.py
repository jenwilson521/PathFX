import pickle,os,csv


analysis_name = 'predicted_test'
drug_name = 'DB00199'

### call the algorithm without phenotype clustering
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
os.system(cmd)

