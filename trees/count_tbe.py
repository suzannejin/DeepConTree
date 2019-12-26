
__description__ = '''

Given a tree with the TBE support of another tree (computed using Booster),

it calculates the average TBE and the TBE of the deepest node. 

'''



# ===========
# READ FILES
# ===========

def read_tree(tree_file):
    with open(tree_file) as tf:
        for line in tf:
            tree=line.rstrip()
            tree=tree[0:-1]  # Exclude ;
    if tree[-3:]=='100':
        tree=tree[0:-3]   # Exclude 100. Note that the TOL doesn't have bootstrap support. 
    return tree


# ===========
# TREE INFOS
# ===========    
    
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

def newick2nodes(tree):
    nodes={}
    all_nodes=[]
    
    # Get leaves
    pr=re.sub(':-?[0-9]+\.[0-9]+','',tree)  # Remove edge distance ,   Remove 0.10617, 0.07338, ... from (HUMAN:0.10617,(MOUSE:0.07338,BOVIN:0.12290)46:0.03801)
    pr=re.sub('\(','',pr)  
    pr=re.sub( '\)[0-9]*[\.0-9]*','',pr)  #pr=re.sub('\)[0-9]*','',pr)
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
    clades=re.findall('[A-Z\/:0-9.,-_]+',tree) #clades=re.findall('[A-Z:0-9.,-]+',tree)
    acc=len(leaves)
    tmp=0
    t=tree
    for j in range(len(clades)):
        c=clades[j].split(',')
        for p in c:
            if p!='':
                p=p.split(':')
                m = re.search('[A-Z]+', p[0])  # The leaves are represented as something like HUMAN:0.29, while the nodes are like ):0.53
                if m==None:
                    acc=acc+1
                    nodes[acc-1].distance=p[1]
                    nodes[acc-1].bootstrap=p[0]
                    t=re.sub(str(p[0]+':'+p[1]),str(acc),t)
                else:
                    tmp=tmp+1
                    nodes[tmp-1].distance=p[1]
                    t=re.sub(str(p[0]+':'+p[1]),str(tmp),t)
    #print(t)
    
    # Reformat simplified tree t to something like: 012345nx6nxn6n5nnxn4nnxn3nnxn2nnxn1nnxn0 
    count=65  # ASCII starting from "A"
    par=''
    for i in range(len(t)):  # For character in simplified tree t
        if t[i]=='(':
            par=par+str(unichr(count))
            count=count+1
        elif t[i]==')':
            count=count-1
            par=par+str(unichr(count))
        elif t[i]==',':  # Comma
            par=par+'x'
        else:
            par=par+'n'  # Leaf
    #print(par)

    # For each internal node, get parent, left & right children 
    for i in range(len(par)):
       for j in range(i+1, len(par)):
           if re.match('[A-Z]',par[i]) and par[i]==par[j]:  # If the same clade
               ##print i, j
               if j!=len(t)-1:
                   P=t[j+1]  # Get leaf (next to "(", so +1)
                   if j+2<len(t):
                       if par[j+2]=='n':
                           P=t[j+1:j+3]
                           ##print P
                   no=ord(par[i])  # Node number
                   cl=t[i:j+1]     # Clade
                   p_cl=par[i:j+1]
                   ind=int(P)-1
                   ##print ind
                   par=par[0:i]+'k'+par[i+1:j]+'k'+par[j+1:]  # Replace node number by k
                   #print cl,'parent', P
                   if re.findall('[A-Z]',p_cl[1:-1]):
                       # Find next internal node
                       s_i_c=p_cl.find(str(unichr(no+1)))
                       if s_i_c!=0:
                          e_i_c=p_cl[s_i_c+1:].find(str(unichr(no+1)))
                          e_i_c=s_i_c+e_i_c
                          n_cl=cl[:s_i_c]+cl[e_i_c+2:]
                          p_cl=p_cl[:s_i_c]+p_cl[e_i_c+2:]
                       elif s_i_c==0:
                          cl=cl[1:]
                          e_i_c=p_cl[1:].find(str(unichr(no+1)))
                          n_cl=cl[e_i_c+1:]
                          p_cl=p_cl[e_i_c+1:]
                   else:
                       n_cl=cl
                   #print(s_i_c,e_i_c,n_cl,p_cl,cl)
                   if re.match('[0-9]',n_cl[2]):
                       l=n_cl[1:3]
                   else:
                       l=n_cl[1]
                   if re.match('[0-9]',n_cl[-3]):
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
    return(nodes,len(leaves))
    
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

def shallowness_of_nodes(nodes,no_l):
    for n in nodes.iterkeys():
        if n<no_l:
            nodes[n].shallowness=0
        else:
            X= min((no_l-len(nodes[n].in_partition)),len(nodes[n].in_partition))
            if X!=0:
                nodes[n].shallowness=X
                #nodes[n].shallowness=(1.0/comb(no_l,X))
    return(nodes)



# =====
# MAIN
# =====
    
def get_avg_tbe(tree):
    import re
    tbe=re.findall("\)([0-9]+\.[0-9]+)",tree)
    total=0.0
    for t in tbe:
        total+=float(t)
    avgtbe=total/len(tbe)
    return(avgtbe)
    
def get_deepest_node_tbe(nodes):
    maxshal = 0
    maxn = None
    for n in nodes.iterkeys():
        if nodes[n].shallowness>maxshal:
            maxshal=nodes[n].shallowness
            maxn=n
    return(nodes[maxn].bootstrap)
    
    
if __name__=='__main__':
    
    import re
    import sys
    import argparse
    #import math
    from scipy.special import comb
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("tree",type=str,help="The tree with the TBE support of another tree (computed using Booster).")
    args=app.parse_args()
    
    tree=read_tree(args.tree)
    nodes,no_l=newick2nodes(tree)
    nodes=shallowness_of_nodes(nodes,no_l)
    
    avgtbe=get_avg_tbe(tree)
    deeptbe=get_deepest_node_tbe(nodes)
    
    #print(tree)
    #for n in nodes.iterkeys():
    #    if n>=no_l:
    #        print(nodes[n].name,nodes[n].left,nodes[n].right,(no_l-len(nodes[n].in_partition)),len(nodes[n].in_partition),nodes[n].shallowness)

    print(avgtbe,deeptbe)
    
    