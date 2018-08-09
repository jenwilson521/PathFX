# written to create a clustergram of all
# phenotypes discovered by PathFX
# re-written 5-25-18 JLW

import csv, pickle, os
import matplotlib
matplotlib.use("AGG")
from scipy.cluster.hierarchy import dendrogram,set_link_color_palette
from fastcluster import linkage
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex, colorConverter
import seaborn as sns
import pandas as pd
import argparse

def separate_links(lk,link_dist,nl):
	node_list = nl
	node_dic = defaultdict(list)
	while len(node_list) > 0:
		for (i,[n1,n2,ndist,num]) in enumerate(lk):
			if ndist <= link_dist:
				new_node = i + len(df) # new nodes are N...2N-2
				if n1 < len(df): # the node is still a single row from data-frame
					node1 = [df.index[int(n1)]]
					node_list.remove(df.index[int(n1)])
				else: # the node is a list of merged, single rows
					node1 = node_dic[n1]
					node_dic.pop(n1,'')
				if n2 < len(df):
					node2 = [df.index[int(n2)]]
					node_list.remove(df.index[int(n2)])
				else:
					node2 = node_dic[n2]
					node_dic.pop(n2,'')
				node_dic[new_node] =  node1+node2
			elif n1 < len(df) or n2 < len(df): # these are singletons that don't get merged
				if n1 < len(df):
					new_node = n1
					node1 = [df.index[int(n1)]]
					node_dic[new_node] = node1
					node_list.remove(df.index[int(n1)])
				if n2 < len(df):
					new_node = n2
					node2 = [df.index[int(n2)]]
					node_dic[new_node] = node2
					node_list.remove(df.index[int(n2)])
	return node_dic

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-f', action='store', dest='linf', help='the file containing a pandas data matrix of similarity values')
parser.add_argument('-a','--analysis_name',dest='aname',help='Name of analysis, will be appended to output files; experiment date is suggested')
parser.add_argument('-d','--dir',dest='res_dir',help='Results directory. If none provided, a directory will be created matching the analysis name in the ../results/ dir')
args = parser.parse_args()

if args.res_dir is None:
        rdir = os.path.join('..','results',args.aname)
else:
        rdir = args.res_dir

# Name the analysis
tname = args.aname
print('loading semantic similarity results for plotting')
# Load the data frame, process
full_name = os.path.join(rdir,args.linf)
#df = pickle.load(open(full_name,'rb'))
df = pd.read_table(full_name,delimiter='\t',header=0)#,na_values="NULL")
df = df.dropna(axis=1,how='all')
df = df.fillna(-1)

# caluclate linkage
link = linkage(df, metric='euclidean', method='ward')
# link_distances = [(i,ndist) for (i,[n1,n2,ndist,num]) in enumerate(link)]
# for lc in [2.2,2.7,3.2,3.7,4.2]:
link_cutoff = 1.7 # set for now, automate later?
print('Separating the tree with: '+str(link_cutoff))

# separate the clusters
separated = separate_links(link,link_cutoff,list(df.index))
print('Saved cluster membership')
print('Number of new clusters: '+str(len(separated)))
outfname = os.path.join(rdir,'disease_clusters_lin_'+str(link_cutoff)+'.pkl')
pickle.dump(separated,open(outfname,'wb'))
p_cut_off = len(separated)

# Mapping to get phenotypes in a readable format
c2phen = pickle.load(open('../rscs/cuis_to_all_phens.pkl','rb'))
rc2ph = pickle.load(open('../rscs/remaining_cui_to_phen.pkl','rb'))

def get_ph_str(cui_list):
	ph_all = []
	for c in cui_list:
		if c in c2phen:
			ph_list = c2phen[c]
		elif c in rc2ph:
			ph_list = [rc2ph[c]]
		else:
			ph_list = [c]		
		for ph in ph_list:
			ph_all.append(ph)

	ph_all = list(set(ph_all))
	return '|'.join(ph_all)

def get_summary_words(cui_list):
	word_counts = defaultdict(int)
	for c in cui_list:
		# find the appropriate source of phenotypes
		if c in c2phen:
			ph_list = c2phen[c]
		elif c in rc2ph:
			ph_list = [rc2ph[c]]
		else:
			ph_list = [c]

		# clean up each word, remove numbers or sinle letters
		for ph in ph_list:
			for w in ph.split(' '):
				if '(' in w and ')' in w: # these are OMIM associations with single letters
					pass
				else: # these are real world and might need to be cleaned up
					if '{' in w or '}' in w or ',' in w or '?' in w:
						w = w.replace('{','').replace('}','').replace(',','').replace('?','')
					w = ''.join([i for i in w if not i.isdigit()]) # remove digits
					if len(w) >2:
						word_counts[w.lower()]+=1
	if len(word_counts) > 0: # if words are found, clean up and count
		(top_words,counts) = zip(*sorted(word_counts.items(),key=lambda x:x[1],reverse=True))
		return ','.join(top_words[0:5])
	else:
		return '|'.join(cui_list)
	
# write the cluster to file to see if they make sense
outfname = os.path.join(rdir,'cluster_membership_'+str(link_cutoff)+'_'+str(p_cut_off)+'.txt')
outf = open(outfname,'w')
for (clid,cui_list) in separated.items():
	ph_string = get_ph_str(cui_list)
	outf.write('\t'.join([str(clid),ph_string,'\n']))
outf.close()

# plot subsets of the dendrogram
fig,ax = plt.subplots(figsize=(10,10))
den = dendrogram(link,labels=df.index,orientation='right',distance_sort=True,truncate_mode='lastp',p=p_cut_off)
ax.yaxis.set_label_position("right")
ax.set_title(tname+'\nUnlabeled dendrogram')
no_spine = {'left': True, 'bottom': True, 'right': True, 'top': True}
sns.despine(**no_spine)
ax.yaxis.set_label_position("right")
plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9, left=0.5)
outfname = os.path.join(rdir,tname+'_unlabeled_dendogram_full_'+str(link_cutoff)+'.png')
plt.savefig(outfname,format='png',dpi=300)
plt.close()

# for merged, take the intersection of their top words
leaf_label_dic = {}
for (leaf_key,cui_list) in separated.items():
	top_words = get_summary_words(cui_list)
	leaf_label_dic[leaf_key]=top_words

# mapping function to pass to dendrogram
def llf(xx):
	if xx in leaf_label_dic:
		return leaf_label_dic[xx]
	else:
		print(xx)
		return 'CHECK THIS'

# plot with leaf labels
fig,ax = plt.subplots(figsize=(10,10))
den = dendrogram(link,leaf_label_func=llf,orientation='right',distance_sort=True,truncate_mode='lastp',p=p_cut_off)
ax.yaxis.set_label_position("right")
ax.set_title(tname+'\nMost common words labeling')
no_spine = {'left': True, 'bottom': True, 'right': True, 'top': True}
sns.despine(**no_spine)
ax.yaxis.set_label_position("right")
plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9, left=0.5)
plt.tick_params(axis='both', which='major', labelsize=10)
outfname = os.path.join(rdir,tname+'_labeledClusters_dendogram_full_'+str(link_cutoff)+'.png')
plt.savefig(outfname,format='png',dpi=300)
plt.close()

print('Finished plotting')

