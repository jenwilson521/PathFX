# written to create gene-cui files
# for JAMIA examples
# re-written 3-10-21 JLW

import pickle,os,csv

g2c = pickle.load(open('../rscs/merged_genes_to_cuis.pkl','rb'))

for gname in ['IL1R2','FCGR2B']:
	cui_terms = g2c[gname]
	outf = open(gname+'_rel_cuis.txt','w')
	n = outf.write('\n'.join(cui_terms))
	outf.close()
