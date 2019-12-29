#!/usr/bin/env python
import re
import sys
#import math
from scipy.special import comb


def readtree(tree_file):
    with open(tree_file) as tf:
        for line in tf:
            tree=line.rstrip()
            tree=tree[0:-1]  # Exclude ;
    if tree[-3:]=='100':
        tree=tree[0:-3]   # Exclude 100. Note that the TOL doesn't have bootstrap support. 
    return tree

def read_msa(msa_file):
    msa={}
    with open(msa_file) as mf:
        #reads line by line the file without running it all in one shot on the RAM
        for line in mf:
            line=line.rstrip()
            if line[0]=='>':
                upid=line[1:]
                continue
            msa[upid]=msa.get(upid,'')+line  # Get multi-line fasta sequences
    return msa


class Node:
    def __init__(self):
        self.name=""
        self.label=''
        self.distance=float (0)
        self.bootstrap=float (-1)
        self.identity=float(0)
        self.left=0
        self.right=0
        self.parent=0
        self.children=[]
        self.extra=[]
        self.profile=None
        self.level=0
        self.shallowness=0
        self.below=[]
        self.above=[]
        self.in_partition=None
        self.tree_compliance=''
        self.shared='0'


def partition(node, nodes):
    part_in=[]
    for i in node.children:
        i=int(i)-1
        if nodes[i].level==1:
            ##print i
            node=nodes[i]
            for j in node.children:
                part_in.append(nodes[j].name)
                #below.append(nodes[j].label)
            ##print part_in
            ##print below
        else:
            ##print i
            node=nodes[i]
            if nodes[i].children!=None:
                ##print nodes[i].children
                part_in= part_in + partition(node, nodes)
                #below.append(nodes[i].label)
                ##print part_in
                ##print below
            else:
                part_in.append(nodes[i].name)
                #below.append(nodes[i].label)
                ##print part_in
                ##print below
    return part_in

def below(node, nodes):
        low=[]
        for i in node.children:
            i=int(i)-1
            if nodes[i].level==1:
                #print i
                node=nodes[i]
                for j in node.children:
                    low.append(nodes[j].label)
                    ##print low
            else:
                ##print i
                node=nodes[i]
                if nodes[i].children!=None:
                    low.append(nodes[i].label)
                    low=low+below(node, nodes)
                    ##print low
                else:
                    low.append(nodes[i].label)
                    ##print low
        return low


def newick2nodes(tree, msa):
    nodes={}
    all_nodes=[]
    
    # Get leaves
    pr=re.sub(':-?[0-9]+\.[0-9]+','',tree)  # Remove edge distance ,   Remove 0.10617, 0.07338, ... from (HUMAN:0.10617,(MOUSE:0.07338,BOVIN:0.12290)46:0.03801)
    pr=re.sub('\(','',pr)  
    pr=re.sub('\)[0-9]*','',pr) 
    leaves=pr.split(',')

    # Get subtrees (internal nodes + leaves) name, label, distance & bootstrap (if internal node) 
    # Simplify tree to something like:  ((((((1,(2,3)9)10,4)11,5)12,6)13,7)14,8)
    n_nodes=len(leaves)+len(re.findall('\(',tree))  # Number of leaves + number of internal nodes
    for i in range(n_nodes):
        nodes[i]=Node()
        if i<len(leaves):  # If leaves
            nodes[i].name=leaves[i]  # Leaf name: specie / taxa
            nodes[i].label=str(i+1)
            nodes[i].bootstrap=None
            nodes[i].left=None
            nodes[i].right=None
            nodes[i].children=None
            nodes[i].level=0
            nodes[i].shallowness=0
        else:   # If internal nodes
            nodes[i].name=str(i+1)
            nodes[i].label=str(i+1)
    clades=re.findall('[A-Z:0-9.,-]+',tree)
    acc=len(leaves)
    tmp=0
    t=tree
    for j in range(len(clades)):
        c=clades[j].split(',')
        ###print c
        for p in c:
            if p!='':
                ###print p
                p=p.split(':')
                m = re.search('[A-Z]+', p[0])  # The leaves are represented as something like HUMAN:0.29, while the nodes are like ):0.53
                if m==None:
                    ###print acc
                    acc=acc+1
                    nodes[acc-1].distance=p[1]
                    nodes[acc-1].bootstrap=p[0]
                    t=re.sub(str(p[0]+':'+p[1]),str(acc),t)
                else:
                    tmp=tmp+1
                    nodes[tmp-1].distance=p[1]
                    t=re.sub(str(p[0]+':'+p[1]),str(tmp),t)
    print t
    
    # Reformat simplified tree t to something like: 012345nx6nxn6n5nnxn4nnxn3nnxn2nnxn1nnxn0 
    count=0
    par=''
    for i in range(len(t)):  # For character in simplified tree t
        if t[i]=='(':
            par=par+str(count)
            count=count+1
        elif t[i]==')':
            count=count-1
            par=par+str(count)
        elif t[i]==',':  # Comma
            par=par+'x'
        else:
            par=par+'n'  # Leaf
    print par
    
    # For each internal node, get parent, left & right children 
    for i in range(len(par)):
       for j in range(i+1, len(par)):
           if par[i].isdigit() and par[i]==par[j]:  # If the same clade
               ##print i, j
               if j!=len(t)-1:
                   P=t[j+1]  # Get leaf (next to "(", so +1)
                   if j+2<len(t):
                       if par[j+2]=='n':
                           P=t[j+1:j+3]
                           ##print P
                   no=int(par[i])  # Number node
                   cl=t[i:j+1]     # Clade
                   p_cl=par[i:j+1]
                   ind=int(P)-1
                   ##print ind
                   par=par[0:i]+'k'+par[i+1:j]+'k'+par[j+1:]  # Replace node number by k
                   #print cl,'parent', P
                   if re.findall('[0-9]',p_cl[1:-1]):
                       # Find next internal node
                       s_i_c=p_cl.find(str(no+1))
                       if s_i_c!=0:
                          e_i_c=p_cl[s_i_c+1:].find(str(no+1))
                          e_i_c=s_i_c+e_i_c
                          n_cl=cl[:s_i_c]+cl[e_i_c+2:]
                          p_cl=p_cl[:s_i_c]+p_cl[e_i_c+2:]
                       elif s_i_c==0:
                          cl=cl[1:]
                          e_i_c=p_cl[1:].find(str(no+1))
                          n_cl=cl[e_i_c+1:]
                          p_cl=p_cl[e_i_c+1:]
                   else:
                       n_cl=cl
                   #print(s_i_c,e_i_c,n_cl,p_cl,cl)
                   #print n_cl
                   #print p_cl
                   if n_cl[2].isdigit():
                       l=n_cl[1:3]
                   else:
                       l=n_cl[1]
                   if n_cl[-3].isdigit():
                       r=n_cl[-3:-1]
                   else:
                       r=n_cl[-2]
                   #print ind,l,r,P
                   nodes[ind].left=l
                   nodes[ind].right=r
                   nodes[int(l)-1].parent=P
                   nodes[int(r)-1].parent=P
                   nodes[ind].children=[]
    for i in nodes.iterkeys():
        if nodes[i].right or nodes[i].left:  # If internal node, so it has left/right children
            nodes[i].children=[nodes[i].left, nodes[i].right]
    r_ind=max(nodes.keys())  # Max node index
    r=nodes[r_ind].label     # Max node label
    for i in nodes.keys():
        if i!=r_ind:
            if nodes[i].parent==0:
                nodes[i].parent=r

    # Organize the root
    to_root=[]
    for i in nodes.iterkeys():
        if nodes[i].parent==r:
            to_root.append(i)  # The nodes/leaves directly connected to root
    nodes[r_ind].children=[]
    nodes[r_ind].name='root'
    missing=[]
    for i in to_root:
        nodes[r_ind].children.append(i+1)  # Get the children of root
        if nodes[i].children==None:
            nodes[r_ind].extra.append(i+1) # Extra children of root
        else:
            missing.append(i)

    if len(missing)>=2:
        l_c=min(nodes[i].children for i in missing)
        r_c=max(nodes[i].children for i in missing)
        for i in missing:
            if nodes[i].children==l_c:
                nodes[r_ind].left=str(i+1)
            elif nodes[i].children==r_c:
                nodes[r_ind].right=str(i+1)
            else:
                nodes[r_ind].extra.append(i+1)
    else:
        l=min(missing)+1
        r=min(nodes[r_ind].extra)
        nodes[r_ind].left=l
        nodes[r_ind].right=r
        nodes[r_ind].extra.remove(r)
    nodes[r_ind].bootstrap=None
    for i in nodes.iterkeys():
        all_nodes.append(nodes[i].label)
    all_nodes=set(all_nodes)
    for i in nodes.iterkeys():
        if int(nodes[i].label)>len(leaves):
            nodes[i].in_partition=partition(nodes[i], nodes)
            nodes[i].below=below(nodes[i], nodes)
            bl=set(nodes[i].below)
            ab=list(all_nodes-bl)
            ab.remove(nodes[i].label)
            nodes[i].above=ab
    for i in nodes.iterkeys():
        if nodes[i].bootstrap=='':
            nodes[i].bootstrap=0
    return nodes

def seq_to_profile(seq):
    aa='ACDEFGHIKLMNPQRSTVWY-'
    profile=[[0.0 for j in range(len(seq))]for i in range(len(aa))]
    for j in range(len(seq)):
        for i in range(len(aa)):
            if seq[j]==aa[i]:
                profile[i][j]=profile[i][j]+1
    return profile

def compute_perc_identity(prof1, prof2):
    if len(prof1)==len(prof2):
        match=0
        for j in range(len(prof1[0])):
            m=0
            for i in range(len(prof1)):
                if prof1[i][j]==prof2[i][j]:
                    m=m+1
            if m==len(prof1):
                match=match+1
        match=(float(match)/len(prof1[0]))*100
    return match

def profiles2profile(p1, p2):
    if len(p1)==len(p2):
        profile=[[0.0 for j in range(len(p1[0]))] for i in range(len(p1))]
        for j in range(len(p1[0])):
            for i in range(len(p1)):
                ##print i,j
                profile[i][j]=(p1[i][j]+p2[i][j])/2
    return profile

def msa_on_tree(msa,nodes):
    ##print msa.keys()
    for n in nodes.iterkeys():
        #nodes[n].profile=[[0.0 for i in len(msa[n])] for j in range(21)]
        if nodes[n].name in msa.keys():  # If leaf
            nodes[n].profile=seq_to_profile(msa[nodes[n].name])
            nodes[n].identity=None
        else:  # If internal node
            ##print n, nodes[n].name
            l=int(nodes[n].left)-1
            r=int(nodes[n].right)-1
            nodes[n].profile=profiles2profile(nodes[l].profile, nodes[r].profile)
            nodes[n].identity=compute_perc_identity(nodes[l].profile, nodes[r].profile)
    return nodes

def level_of_nodes(nodes):
    for n in nodes.iterkeys():
        if nodes[n].left and nodes[n].right:
            nodes[n].level=int(min(nodes[int(nodes[n].left)-1].level, nodes[int(nodes[n].right)-1].level))+1
        elif nodes[n].left!=None and nodes[n].right==None:
            nodes[n].level=int(nodes[int(nodes[n].left)-1].level)+1
        elif nodes[n].right==None and nodes[n].left!=None:
            nodes[n].level=int(nodes[int(nodes[n].right)-1].level)+1
        else:
            nodes[n].level=0

    ml=max(nodes[n].level for n in nodes.keys())
    for n in nodes.iterkeys():
        nodes[n].level=float(nodes[n].level/ml)  # Normalize by the maximum level ?? So level between [0,1]?
                                                 # Question: it will be 0 or 1, since level was set as integer
    return nodes

# def shallowness_of_nodes(nodes, no_l):
#     r_ind=max(nodes.keys())
#     #print r_ind
#     for n in nodes.iterkeys():
#         if n<no_l:
#             nodes[n].shallowness=0
#         elif n!=r_ind:
#             nodes[n].shallowness=int(min(min(nodes[int(x)-1].shallowness for x in nodes[n].above), min(nodes[int(y)-1].shallowness for y in nodes[n].below)))+1
#         else:
#             #print 'computing root'
#             nodes[n].shallowness=int(min(nodes[int(y)-1].shallowness for y in nodes[n].below))+1
#     return nodes

# def shallowness_of_nodes(nodes, no_l):
#     r_ind=max(nodes.keys())
#     ##print r_ind
#     for n in nodes.iterkeys():
#         if n<no_l:
#             nodes[n].shallowness=0
#         elif n!=r_ind:
#             P=int(nodes[n].parent)-1
#             pc=nodes[P].children
#             for c in pc:
#                 if c==n+1:
#                     #print pc, c, n
#                     pc.remove(c)
#             nodes[n].shallowness=int(min(min(nodes[int(x)-1].shallowness for x in pc), min(nodes[int(y)-1].shallowness for y in nodes[n].children), nodes[P].shallowness))+1
#         else:
#             nodes[n].shallowness=int(min(nodes[int(y)-1].shallowness for y in nodes[n].children))+1
#     return nodes

def shallowness_of_nodes(nodes,no_l):
    for n in nodes.iterkeys():
        if n<no_l:
            nodes[n].shallowness=0
        else:
            X= min((no_l-len(nodes[n].in_partition)),len(nodes[n].in_partition))
            if X!=0:
                nodes[n].shallowness=(1.0/comb(no_l,X))
                ##print n, X,(1.0/comb(no_l,X),10)
            #print(nodes[n].name,nodes[n].left,nodes[n].right,no_l-len(nodes[n].in_partition),len(nodes[n].in_partition),X)
    return nodes

def taxa2nodes(tree):
    nodes={}
    all_nodes=[]
    le=tree.replace('(','')
    le=le.replace(')','')
    leaves=le.split(',')
    n_nodes=len(leaves)+len(re.findall('\(',tree))
    for i in range(n_nodes):
        nodes[i]=Node()
        if i<len(leaves):
            nodes[i].name=leaves[i]
            nodes[i].label=str(i+1)
            nodes[i].left=None
            nodes[i].right=None
            nodes[i].children=None
        else:
            nodes[i].name=str(i+1)
            nodes[i].label=str(i+1)
    t=tree
    for i in range(len(leaves)):
        t=t.replace(leaves[i], nodes[i].label)
    st=len(leaves)
    tr=''
    for i in t:
        tr=tr+i
        if i==')':
            st=st+1
            tr=tr+str(st)
    t=tr
    #print tree
    print t
    count=0
    par=''
    for i in range(len(t)):
        if t[i]=='(':
            count=count+1
            par=par+str(count)
        elif t[i]==')':
            par=par+str(count)
            count=count-1
        elif t[i]==',':
            par=par+'x'
        else:
            par=par+'n'
    print(par)
    for i in range(len(par)):
       for j in range(i+1, len(par)):
           if par[i].isdigit() and par[i]==par[j]:
               #print i, j
               if j!=len(t)-1:
                   P=t[j+1]
                   if j+2<len(t):
                       if par[j+2]=='n':
                           P=t[j+1:j+3]
                           ##print P
                   no=int(par[i])
                   cl=t[i:j+1]
                   p_cl=par[i:j+1]
                   ind=int(P)-1
                   #print ind
                   par=par[0:i]+'k'+par[i+1:j]+'k'+par[j+1:]
                   ##print cl,'parent', P
                   if re.findall('[0-9]',p_cl[1:-1]):
                       s_i_c=p_cl.find(str(no+1))
                       if s_i_c!=0:
                          e_i_c=p_cl[s_i_c+1:].find(str(no+1))
                          e_i_c=s_i_c+e_i_c
                          ##print e_i_c
                          n_cl=cl[:s_i_c]+cl[e_i_c+2:]
                          p_cl=p_cl[:s_i_c]+p_cl[e_i_c+2:]
                       elif s_i_c==0:
                          cl=cl[1:]
                          e_i_c=p_cl[1:].find(str(no+1))
                          ##print e_i_c
                          n_cl=cl[e_i_c+1:]
                          p_cl=p_cl[e_i_c+1:]
                   else:
                       n_cl=cl
                   #print n_cl
                   #print p_cl
                   if n_cl[2].isdigit():
                       l=n_cl[1:3]
                   else:
                       l=n_cl[1]
                   if n_cl[-3].isdigit():
                       r=n_cl[-3:-1]
                   else:
                       r=n_cl[-2]
                   #print l,r,P
                   nodes[ind].left=l
                   nodes[ind].right=r
                   nodes[int(l)-1].parent=P
                   nodes[int(r)-1].parent=P
                   nodes[ind].children=[]
                   next_p=str(int(p_cl[0])+1)
                   if p_cl.find(str('nxnnxn'+p_cl[0]))!=-1 or p_cl.find(str('nxnxn'+p_cl[0]))!=-1 or p_cl.find(str('nxnnxnn'+p_cl[0]))!=-1 or p_cl.find(str('nxnxnn'+p_cl[0]))!=-1:
                       #print 'not binary split, no extra p'
                       #print p_cl
                       s=p_cl.find('nx')
                       if p_cl[s-1]=='n':
                           s=s-1
                       #print s
                       p_cl=p_cl[s:]
                       #print p_cl
                       ex_cl=n_cl[s:-1]
                       #print ex_cl
                       ex_cl=ex_cl.split(',')
                       #print ex_cl
                       for n in ex_cl:
                           if n!=l and n!=r:
                               nodes[ind].extra.append(n)
                               nodes[ind].children.append(n)
                               nodes[int(n)-1].parent=P
                   if p_cl.find(str('nxnx'+next_p))!=-1 or p_cl.find(str('nxnnx'+next_p))!=-1:
                       #print 'not binary split'
                       s=p_cl.find(next_p)
                       e=p_cl[s+1:].find(next_p)
                       e=s+e+2
                       ex_cl=n_cl[1:s]+n_cl[e:-1]
                       ex_cl=ex_cl.split(',')
                       for n in ex_cl:
                           if n!=l and n!=r:
                               nodes[ind].extra.append(n)
                               nodes[ind].children.append(n)
                               nodes[int(n)-1].parent=P
    for i in nodes.iterkeys():
        if nodes[i].right or nodes[i].left:
            nodes[i].children=[nodes[i].left,nodes[i].right]
            nodes[i].children=nodes[i].children+nodes[i].extra
    for i in nodes.iterkeys():
        all_nodes.append(nodes[i].label)
    all_nodes=set(all_nodes)
    for i in nodes.iterkeys():
        if int(nodes[i].label)>len(leaves):
            nodes[i].in_partition=partition(nodes[i], nodes)
            nodes[i].below=below(nodes[i], nodes)
            bl=set(nodes[i].below)
            ab=list(all_nodes-bl)
            ab.remove(nodes[i].label)
            nodes[i].above=ab
    for i in nodes.iterkeys():
        if nodes[i].in_partition!=None and len(nodes[i].in_partition)==len(leaves):
            nodes[i].name='root'
    return nodes

def read_seq2taxa(seq2taxa_file):
    seq2taxa={}
    with open(seq2taxa_file) as stf:
        for line in stf:
            line=line.rstrip().split('\t')
            seq2taxa[line[0]]=line[1]
    return seq2taxa

def converting_names(nodes,no_l,seq2taxa):
    for i in nodes.iterkeys():
        if int(nodes[i].label)<=no_l:
            nodes[i].name=seq2taxa[nodes[i].name]
        else:
            for j in range(len(nodes[i].in_partition)):
                nodes[i].in_partition[j]=seq2taxa[nodes[i].in_partition[j]]
    return nodes

def classes2bootstrap(nodes1d, nodes3d, no_l, t_nodes, tree_file):
    multisplit=[]
    binary=[]
    f1d=[]
    nf1d=[]
    ms1d=[]
    f1d3d=[]
    f3d1d=[]
    f3d=[]
    nf3d=[]
    ms3d=[]
    for k in t_nodes.iterkeys():
        if t_nodes[k].children!=None:
            if len(t_nodes[k].children)>2 and t_nodes[k].name!='root':
                multisplit.append(sorted(t_nodes[k].in_partition))
            else:
                binary.append(sorted(t_nodes[k].in_partition))
    for i in nodes1d.iterkeys():
        for j in nodes3d.iterkeys():
            if int(nodes1d[i].label)>int(no_l)and nodes1d[i].name!='root':  # If internal node
                if nodes3d[j].in_partition!=None and sorted(nodes1d[i].in_partition)== sorted(nodes3d[j].in_partition):  # If the same branch
                    f1d3d.append(nodes1d[i].bootstrap)
                    f3d1d.append(nodes3d[j].bootstrap)
                    nodes1d[i].shared='S'
                    nodes3d[j].shared='S'
    for i in nodes1d.iterkeys():
        if int(nodes1d[i].label)>int(no_l) and nodes1d[i].name!='root':
            #print nodes1d[i].label
            if sorted(nodes1d[i].in_partition) in binary:
                f1d.append(nodes1d[i].bootstrap)
                #print str(tree_file[0:7])
                nodes1d[i].tree_compliance='F'
            elif sorted(nodes1d[i].in_partition) in multisplit:
                ms1d.append(nodes1d[i].bootstrap)
                #print str(tree_file[0:7])
                nodes1d[i].tree_compliance='M'
            else:
                nf1d.append(nodes1d[i].bootstrap)
                nodes1d[i].tree_compliance='N'    
    for i in nodes3d.iterkeys():
        if int(nodes3d[i].label)>int(no_l)and nodes3d[i].name!='root':
            #print nodes3d[i].label
            if sorted(nodes3d[i].in_partition) in binary:
                f3d.append(nodes3d[i].bootstrap)
                #print str(tree_file[0:7])
                nodes3d[i].tree_compliance='F'
            elif sorted(nodes3d[i].in_partition) in multisplit:
                ms3d.append(nodes3d[i].bootstrap)
                #print str(tree_file[0:7])
                nodes3d[i].tree_compliance='M'
            else:
                nf3d.append(nodes3d[i].bootstrap)
                nodes3d[i].tree_compliance='N'
    return f1d, nf1d, ms1d,f3d, nf3d, ms3d,f1d3d,f3d1d



if __name__ == '__main__':
    tree_file1d=sys.argv[1]
    tree_file3d=sys.argv[2]
    msa_file=sys.argv[3]
    taxa_t_file=sys.argv[4]   # TOL tree
    seq2taxa_file=sys.argv[5]   # Seq to taxa list
    out_file_b_id_sh=sys.argv[6]
    out_file=sys.argv[7]
    tree1d=readtree(tree_file1d)
    tree3d=readtree(tree_file3d)
    taxa=readtree(taxa_t_file)
    msa=read_msa(msa_file)
    no_l=len(msa.keys())   # Number sequences
    seq2taxa=read_seq2taxa(seq2taxa_file)
    nodes1d=newick2nodes(tree1d, msa)
    nodes1d=msa_on_tree(msa,nodes1d)
    nodes1d=level_of_nodes(nodes1d)
    nodes1d=shallowness_of_nodes(nodes1d, no_l)
    nodes1d=converting_names(nodes1d,no_l,seq2taxa)
    nodes3d=newick2nodes(tree3d, msa)
    nodes3d=msa_on_tree(msa,nodes3d)
    nodes3d=level_of_nodes(nodes3d)
    nodes3d=shallowness_of_nodes(nodes3d, no_l)
    nodes3d=converting_names(nodes3d,no_l,seq2taxa)
    t_nodes=taxa2nodes(taxa)
    t_nodes=level_of_nodes(t_nodes)
    t_nodes=shallowness_of_nodes(t_nodes, no_l)
    f1d, nf1d, ms1d,f3d, nf3d, ms3d,f1d3d,f3d1d=classes2bootstrap(nodes1d, nodes3d,no_l, t_nodes,tree_file1d)
    out_file_F1d=str(out_file)+'_1d_found'
    out_file_NF1d=str(out_file)+'_1d_not_found'
    out_file_MS1d=str(out_file)+'_1d_multi_split'
    out_f=open(out_file_F1d, 'a+')
    for i in f1d:
        out_f.write(str(i)+'\n')
    out_f.close()
    out_nf=open(out_file_NF1d,'a+')
    for i in nf1d:
        out_nf.write(str(i)+'\n')
    out_nf.close()
    out_ms=open(out_file_MS1d, 'a+')
    for i in ms1d:
        out_ms.write(str(i)+'\n')
    out_ms.close()
    out_file_F3d=str(out_file)+'_3d_found'
    out_file_NF3d=str(out_file)+'_3d_not_found'
    out_file_MS3d=str(out_file)+'_3d_multi_split'
    out_f3d=open(out_file_F3d, 'a+')
    for i in f3d:
        out_f3d.write(str(i)+'\n')
    out_f3d.close()
    out_nf3d=open(out_file_NF3d,'a+')
    for i in nf3d:
        out_nf3d.write(str(i)+'\n')
    out_nf3d.close()
    out_ms3d=open(out_file_MS3d, 'a+')
    for i in ms3d:
        out_ms3d.write(str(i)+'\n')
    out_ms3d.close()
    out_file_F1d3d=str(out_file)+'_found_1d_3d'
    out_f1d3d=open(out_file_F1d3d, 'a+')
    for i in range(len(f1d3d)):
        out_f1d3d.write(str(f1d3d[i])+' '+str(f3d1d[i])+'\n')
    out_f1d3d.close()
    out_file_b_id_sh1d='1d_'+str(out_file_b_id_sh)
    out=open(out_file_b_id_sh1d, 'a+')
    for i in nodes1d.iterkeys():
        if nodes1d[i].bootstrap!=None:
            out.write(str(tree_file1d[0:7])+' '+str(tree_file1d[8:14])+' '+str(nodes1d[i].bootstrap)+' '+str(nodes1d[i].identity)+' '+ str(nodes1d[i].level)+' '+str(nodes1d[i].shallowness)+' '+str(nodes1d[i].tree_compliance)+' '+str(nodes1d[i].shared)+'\n')
    out.close()
    out_file_b_id_sh3d='3d_'+str(out_file_b_id_sh)
    out=open(out_file_b_id_sh3d, 'a+')
    for i in nodes3d.iterkeys():
      if nodes3d[i].bootstrap!=None:
          out.write(str(tree_file3d[0:7])+' '+str(tree_file3d[8:14])+' '+str(nodes3d[i].bootstrap)+' '+str(nodes3d[i].identity)+' '+ str(nodes3d[i].level)+' '+str(nodes3d[i].shallowness)+' '+str(nodes3d[i].tree_compliance)+' '+str(nodes3d[i].shared)+'\n')
    out.close()

