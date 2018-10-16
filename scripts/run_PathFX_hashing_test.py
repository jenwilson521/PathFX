# written 10-16-18 JLW
# testing implementation of hashed files
import pickle,os,csv
from datetime import datetime

analysis_name = 'PathFX_hashing_test'

for drug_name in ['metformin','pravastatin','montelukast','propofol']:
	startTime = datetime.now()
	# drug_name = 'metformin'
	cmd = 'python phenotype_enrichment_pathway_beta.py -d %s -a %s'%(drug_name,analysis_name)
	print(cmd)
	os.system(cmd)
	print(datetime.now() - startTime)

#startTime = datetime.now()
#drug_name = 'Montelukast'
#cmd = 'python phenotype_enrichment_pathway_beta.py -d %s -a %s'%(drug_name,analysis_name)
#print(cmd)
#os.system(cmd)
#print(datetime.now() - startTime)

#drug_name = 'pravastatin'
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
#os.system(cmd)
#
#
#drug_name = 'propofol'
#cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
#os.system(cmd)



