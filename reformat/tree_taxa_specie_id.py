#!/usr/bin/env python
import re
import sys
import argparse

def readtree(tree_file):
    with open(tree_file) as tf:
        for line in tf:
            tree=line.rstrip()
            tree=tree[0:-1]
    return tree

def read_seq2taxa(seq2taxa_file,change):
    if change=="to_specie":
        a,b=1,0
    else:
        a,b=0,1
    seq2taxa={}
    with open(seq2taxa_file) as stf:
        for line in stf:
            line=line.rstrip().split('\t')
            seq2taxa[line[a]]=line[b]
    return seq2taxa

def name_mapping(tree, seq2taxa):
    pr=re.sub(':-?[0-9]+\.[0-9]+','',tree)
    pr=re.sub('\(','',pr)
    pr=re.sub('\)[0-9]*','',pr)
    leaves=pr.split(',')
    for i in leaves:
        tree=re.sub(i, seq2taxa[i], tree)
    return tree

if __name__ == '__main__':
    app = argparse.ArgumentParser(description="Change the name of the tree leaves. From specie name to taxa ID, or vice versa.")
    app.add_argument("tree_file",type=str,help="Input tree file (.nwk)")
    app.add_argument("seq2taxa_file",type=str,help="Table with two columns: specie name, taxa ID")
    app.add_argument("out_file",type=str,help="Output tree file (.nwk)")
    app.add_argument("change",type=str,choices=["to_specie","to_taxa"]) 
    args=app.parse_args()
    tree=readtree(args.tree_file)
    seq2taxa=read_seq2taxa(args.seq2taxa_file,args.change)
    tree_TID=name_mapping(tree, seq2taxa)
    out=open(args.out_file, 'w')
    out.write(str(tree_TID)+';')
    out.close()
