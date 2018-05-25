# written to compute similarity between CUI terms
# compiled using Maulik's code
# written 5-23-18 JLW


import pickle,os,csv,re,sys
import networkx as nx
import numpy as np

sys.setrecursionlimit(10000)

stopWords = set([
    "a", "also", "although", "am", "an", "and", "are",
    "as", "at", "back", "be", "became", "because", "become",
    "becomes", "becoming", "been", "being", "bill", "both",
    "bottom", "but", "by", "call", "can", "con",
    "could", "de", "do", "done", "eg", "etc", "even", "ever", 
    "find", "for", "found", "from", "get", "give", "go",
    "had", "has", "have", "he", "her", "here", "hers", "herself", "him", "himself", "his",
    "how", "however", "if", "in", "inc", 
    "into", "is", "it", "its", "itself", "keep", "may", "me", "mine", "my", "myself", "name", "namely", "of", "onto", "our",
    "ours", "ourselves", "please", "put", "should", "show", "such", "take", "that", "the", "their", "them",
    "themselves", "these", "they", "this", "those", "though",
    "thru", "to", "us", "via", "was", "we", "were", "what", "whatever", "when",
    "whence", "whenever", "where", "whereafter", "whereas", "whereby",
    "wherein", "whereupon", "wherever", "whether", "which", "whither",
    "who", "whoever", "whom", "whose", "why", "will", "would", "yet", "you", "your", "yours", "yourself", "yourselves"])

def get_ancestors_bfs(term, only_cui=True):
    ancestors = set(nx.bfs_predecessors(combG, term).keys())
    if only_cui:
        cui_ancestors = set([])
        for k in ancestors:
            if combG.node[k]["type"] == "CuiNode": cui_ancestors.add(k)
        return cui_ancestors
    else: return ancestors
    
def get_descendants_bfs(term, only_cui=True):
    _bfs = nx.bfs_edges(combG, term, reverse=True)
    descendants = set([k[1] for k in _bfs])
    if only_cui:
        cui_descendants = set([])
        for k in descendants:
            if combG.node[k]["type"] == "CuiNode": cui_descendants.add(k)
        return cui_descendants
    else: return descendants

def get_leaves_bfs(term, only_cui=True):
    descendants = get_descendants_bfs(term, only_cui)
    leaves = set([])
    for k in descendants:
        in_cuis = [m[0] for m in combG.in_edges(k) if combG.node[m[0]]["type"] == "CuiNode"]
        if len(in_cuis) == 0: leaves.add(k)
    return leaves

def get_cui_mapper(iri):
    if combG.node[iri]["type"] == "CuiNode": return set([iri])
    cui_terms = set([])
    for k in combG[iri]:
        if combG.node[k]["type"] == "CuiNode": cui_terms.add(k)
    return cui_terms

def get_lcs(term1, term2, only_cui=True):
    ans1 = get_ancestors_bfs(term1, only_cui)
    ans2 = get_ancestors_bfs(term2, only_cui)
    cmn_ans = ans1.intersection(ans2)
    lcs = min(cmn_ans, key=lambda x: nx.bfs_predecessors(combG, term1)[x])
    return lcs

def get_depth(term, only_cui=True):
    top_nodes = set([])
    ancestors = get_ancestors_bfs(term, only_cui)
    for k in ancestors:
        if len(combG[k]) == 0: top_nodes.add(k)
    if len(top_nodes) == 0: return 0
    shortest_depth = min([nx.shortest_path_length(combG, term, k) for k in top_nodes])
    return shortest_depth

def sanchez_ic(iri, only_cui=True):
    # term = iri
    cui_terms = get_cui_mapper(iri)
    ic_set = []
    for term in cui_terms:
        # print("Computing Sanchez IC for " + term)
        leaves = get_leaves_bfs(term, only_cui)
        max_nodes = all_cui_nodes if only_cui else all_nodes
        ancestors = get_ancestors_bfs(term, only_cui)
        #print len(leaves), len(ancestors), max_nodes
        ic = - np.log((len(leaves)/float(len(ancestors)) + 1)/float(max_nodes+1))
        ic_set.append(ic)
        # print(ic)
    ic = np.mean(ic_set)
    return ic
    
def lin_sim(iri1, iri2, only_cui=True, method=sanchez_ic):
    lcs = get_lcs(iri1, iri2, only_cui)
    ic_lcs = method(lcs, only_cui)
    ic_c1 = method(iri1, only_cui)
    ic_c2 = method(iri2, only_cui)
    sim = 2*ic_lcs/(ic_c1+ic_c2)
    return sim

#combG = nx.read_gpickle("combinedCuiOntoGraph.gpickle")
#lin_sim(iri1, iri2)
