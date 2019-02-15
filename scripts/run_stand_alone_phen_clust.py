# Example to call phenotype clustering stand-alone
# 14 Feb 2019 JLW

import os

analysis_name = 'testing_phen_clust'
rdir = os.path.join('..','results',analysis_name)

if not os.path.exists(rdir):
	os.makedirs(rdir)

cui_list = '../data/CUI_terms_for_phen_clustering_test_clean.txt'  # a new-line delimited file  of CUI terms
lin_file = os.path.join(rdir,analysis_name+'_lin_matrix_out.txt')

# Call the UMLS perl script with configuration file
cmd = "python calc_lin_matrix_umls_SO.py -cui_list %s -outf %s"%(cui_list,lin_file)
print(cmd)
os.system(cmd)

# Call the plotting script once UMLS script is complete
cmd = "python plot_and_cluster_phenotypes_SO.py -f %s -a %s -d %s"%(lin_file,analysis_name,rdir)
print(cmd)
os.system(cmd)
