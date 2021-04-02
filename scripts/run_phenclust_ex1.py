# written to demonstrate an example analysis
# using Metform-network-phenotypes
# re-written 3-11-21 JLW

import pickle,os,csv,datetime 
# track time to get a sense
start = datetime.datetime.now()

cui_list_f = "Metformin_merged_neighborhood__cui_list_.txt"

analysis_name = "Metformin"
rdir = os.path.join('..','results',analysis_name)
if not os.path.exists(rdir):
        os.makedirs(rdir)

cui_list = os.path.join(rdir,cui_list_f)
lin_file = os.path.join(rdir,analysis_name+'_lin_matrix_out.txt')

# Call the UMLS perl script with configuration file
cmd = "python calc_lin_matrix_umls_SO.py -cui_list %s -outf %s"%(cui_list,lin_file)
print(cmd)
os.system(cmd)

# Call the plotting script once UMLS script is complete
cmd = "python plot_and_cluster_phenotypes_SO.py -f %s -a %s -d %s"%(lin_file,analysis_name,rdir)
print(cmd)
os.system(cmd)
print('Total taken for '+analysis_name)
print(datetime.datetime.now() - start)
