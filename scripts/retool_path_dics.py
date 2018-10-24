# written to implement hashing for accessing node dics
# written 10-12-18 JLW

import os,pickle
from collections import defaultdict

# make new dir for storing results
ndir = os.path.join('..','results','node_results')
if not os.path.exists(ndir):
	os.makedirs(ndir)
nsdir = os.path.join('..','results','node_summary_results')
if not os.path.exists(nsdir):
	os.makedirs(nsdir)

new_hash_map = defaultdict(list)
sum_hash_map = {}
unique_nodes = set()

rdir = '../results/new_intome_rands_0.77/'
allf = [(i+1,f) for (i,f) in enumerate(sorted([f for f in os.listdir(rdir)]))] #23853 files

# also look at summary files
srdir = os.path.join(rdir,'summary')
sumf_dic = dict([(f.replace('_0.77_randPathScores_.pkl','').upper(),f) for f in os.listdir(srdir)]) # force all to uppercase

for (counter,f) in allf:
	if 'pth_dic' in f: # and '(' in f:

		print(counter)
		# create new file name and subdirectory
		file_name = 'n'+"{0:0=6d}".format(counter)+'.pkl'
		sdir = os.path.join(ndir,file_name[0:4])
		if not os.path.exists(sdir):
			os.makedirs(sdir)

		# get the node name, assign new file path
		nname = f.replace('pth_dic_','').replace('_scr0.77.pkl','')
		new_f_path = os.path.join(sdir,file_name)
		# print(nname)

		# map all forms of the gene name back to the new file
#		new_hash_map[nname] = new_f_path
		new_hash_map[nname.upper()] = new_f_path
#		new_hash_map[nname.lower()] = new_f_path 

		# error checking for parentheses in file name
		if '(' in f or ')' in f or "'" in f:
			f = f.replace('(','\(').replace(')','\)').replace("'",r"\'")
		# copy the new file to the new location
		cmd = 'cp %s %s' % (os.path.join(rdir,f),new_f_path) # original command
		# print(cmd)
		# os.system(cmd)

		unique_nodes.add(nname.upper())

		# repeat process if the node has a summary file
		if nname.upper() in sumf_dic: 
			sumf = sumf_dic[nname.upper()]
			if '(' in sumf or ')' in sumf or "'" in sumf:
				sumf = sumf.replace('(','\(').replace(')','\)').replace("'",r"\'")
			sumf_path = os.path.join(srdir,sumf)

			# follow the originally assigned node number for naming the summary files
			ssdir = os.path.join(nsdir,'sumn'+file_name[0:4])
			if not os.path.exists(ssdir):
				os.makedirs(ssdir)
			sum_file_name = 'sumn'+"{0:0=6d}".format(counter)+'.pkl'
			new_sf_path = os.path.join(ssdir,sum_file_name)

			cmd = 'cp %s %s' % (sumf_path,new_sf_path)
			# print(cmd)
			# os.system(cmd)
			
			sum_hash_map[nname.upper()] = new_sf_path

	else: # these are duplicates of files with long names
		pass

# manually map a lost file? I don't know that the original node dic was kept, maybe only the version of the file with the abbreviated name
nname = 'CYP2C19star1,CYP2C19star10,CYP2C19star11,CYP2C19star13,CYP2C19star14,CYP2C19star15,CYP2C19star16,CYP2C19star18,CYP2C19star19,CYP2C19star22,CYP2C19star23,CYP2C19star24,CYP2C19star25,CYP2C19star26,CYP2C19star5,CYP2C19star6,CYP2C19star8,CYP2C19star9'
sumf_path = os.path.join(srdir,'CYP2C19star1,5-10,11,13-16,18-19,22-26_0.77_randPathScores_.pkl')
counter+=1
file_name = 'n'+"{0:0=6d}".format(counter)+'.pkl' # but I can't find the actual node pth dic file?
sum_file_name = 'sumn'+"{0:0=6d}".format(counter)+'.pkl'
ssdir = os.path.join(nsdir,'sumn'+file_name[0:4])
new_sf_path = os.path.join(ssdir,sum_file_name)

sum_hash_map[nname.upper()] = new_sf_path
cmd = 'cp %s %s' % (sumf_path,new_sf_path)
print(cmd)
os.system(cmd)

pickle.dump(unique_nodes,open('../rscs/unique_network_nodes.pkl','wb'))
pickle.dump(new_hash_map,open('../rscs/gene_to_hash_map.pkl','wb'))
pickle.dump(sum_hash_map,open('../rscs/gene_to_sum_hash_map.pkl','wb'))	

