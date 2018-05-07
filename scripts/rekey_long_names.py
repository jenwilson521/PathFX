# written to create shorter names
# for nodes (mostly gene variants) with long file names
# written 5-3-18 JLW

import pickle,os,csv,sys

fpath = '../results/new_intome_rands_0.77/'
allf = os.listdir(fpath)

len_thr = 100

long_names = []
for f in allf:
	if len(f) > len_thr:
		long_names.append(f)

file_swap = {}
for (i,lf) in enumerate(sorted(long_names,reverse=True)):
	old_path = os.path.join(fpath,lf)
	
	new_name = '_'.join(['node',str(i),'.pkl'])
	new_path = os.path.join(fpath,new_name)
	fdic = pickle.load(open(old_path,'rb'))
	pickle.dump(fdic,open(new_path,'wb'))
	#cmd = "cp %s %s"%(old_path,new_path)
	#n = os.system(cmd)

	file_swap[old_path] = new_path

pickle.dump(file_swap,open(os.path.join('..','rscs','file_swap.pkl'),'wb'))

### Repeat similar process for summary files
fpath = '../results/new_intome_rands_0.77/summary/'
allf = os.listdir(fpath)
long_names = []
for f in allf:
	if len(f) > len_thr:
		long_names.append(f)
file_swap = {}
for (i,lf) in enumerate(sorted(long_names,reverse=True)):
	old_path = os.path.join(fpath,lf)
	new_name = '_'.join(['node',str(i),'.pkl'])
	new_path = os.path.join(fpath,new_name)
	fdic = pickle.load(open(old_path,'rb'))
	pickle.dump(fdic,open(new_path,'wb'))

	file_swap[old_path]=new_path		

pickle.dump(file_swap,open(os.path.join('..','rscs','summary_file_swap.pkl'),'wb'))

#def clean_node_name(sn):
#        if '+' in sn:
#                sn = sn.replace('+','-')
#        if ':' in sn:
#                sn = sn.replace(':','-')
#        if '/' in sn:
#                sn = sn.replace('/','-')
#        if ' ' in sn:
#                sn = sn.replace(' ','')
#        if '*' in sn:
#                sn = sn.replace('*','star')
#        return sn

	
