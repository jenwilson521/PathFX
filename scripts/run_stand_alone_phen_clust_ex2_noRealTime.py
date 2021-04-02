# Example to call phenotype clustering stand-alone
# new example developed June-2020 JLW
# rewritten Feb 2021 to rerun with UMLS 2020AA

import os, datetime

# track time to get a sense 
start = datetime.datetime.now()

### FIRST FCGR2B
analysis_name = 'FCGR2B_phenotype_analysis_noRealTime'
rdir = os.path.join('..','results',analysis_name)
if not os.path.exists(rdir):
	os.makedirs(rdir)

cui_list_f = 'FCGR2B_rel_cuis.txt'
cui_list = os.path.join(rdir,cui_list_f) # a new-line delimited file  of CUI terms
gname = cui_list_f.replace('_rel_cuis.txt','')
asub_name = 'jamia_ex_'+gname
lin_file = os.path.join(rdir,asub_name+'_lin_matrix_out.txt')

# Call the UMLS perl script with configuration file
cmd = "python calc_lin_matrix_umls_SO_noRealtime.py -cui_list %s -outf %s"%(cui_list,lin_file)
print(cmd)
os.system(cmd)

# Call the plotting script once UMLS script is complete
cmd = "python plot_and_cluster_phenotypes_SO.py -f %s -a %s -d %s"%(lin_file,analysis_name,rdir)
print(cmd)
os.system(cmd)
print('Total taken for '+analysis_name)
print(datetime.datetime.now() - start)


### SECOND IL1R2
start = datetime.datetime.now()
analysis_name = 'IL1R2_phenotype_analysis_noRealTime'
rdir = os.path.join('..','results',analysis_name)
if not os.path.exists(rdir):
	os.makedirs(rdir)

cui_list_f  = 'IL1R2_rel_cuis.txt'
cui_list = os.path.join(rdir,cui_list_f) # a new-line delimited file  of CUI terms
gname = cui_list_f.replace('_rel_cuis.txt','')
asub_name = 'jamia_ex_'+gname
lin_file = os.path.join(rdir,asub_name+'_lin_matrix_out.txt')

# Call the UMLS perl script with configuration file
cmd = "python calc_lin_matrix_umls_SO_noRealtime.py -cui_list %s -outf %s"%(cui_list,lin_file)
print(cmd)
os.system(cmd)

# Call the plotting script once UMLS script is complete
cmd = "python plot_and_cluster_phenotypes_SO.py -f %s -a %s -d %s"%(lin_file,analysis_name,rdir)
print(cmd)
os.system(cmd)

print('Total taken for '+analysis_name)
print(datetime.datetime.now() - start)


