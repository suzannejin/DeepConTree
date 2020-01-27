
__description__ = '''

Given a tree with the TBE support of another tree (computed using Booster),

it calculates the average TBE and the TBE of the deepest node. 

'''

    
def get_avg_tbe(tree):
    import re
    tbe=re.findall("\)([0-9]+\.[0-9]+)",tree)
    total=0.0
    for t in tbe:
        total+=float(t)
    avgtbe=total/len(tbe)
    return(avgtbe)
    
def get_deepest_node_tbe(nodes):
    shal={}
    boot=[]
    for n in nodes.iterkeys():
        shal[n]=nodes[n].shallowness
    maxshal=max(shal.values())
    for n in shal:
        if shal[n]==maxshal:
            boot.append(nodes[n].bootstrap)
    return(min(boot))
    
def count_tbe(tbetree,typ,act,p):
    #import math
    from scipy.special import comb
    from rtrees import Node
    from rtrees import read_tree, newick2nodes, partition, below, shallowness_of_nodes, taxa2nodes
    # Read tree
    tree=read_tree(tbetree)
    if typ=="tree":
        nodes,no_l=newick2nodes(tree,p)
    else:
        nodes,no_l=taxa2nodes(tree,p)
    nodes=shallowness_of_nodes(nodes,no_l) 
    # Print
    if p:
        print(tree)
        for n in nodes.iterkeys():
            if n>=no_l:
                print(nodes[n].name,nodes[n].left,nodes[n].right,(no_l-len(nodes[n].in_partition)),len(nodes[n].in_partition),nodes[n].shallowness)
    # Calculate TBE
    if args.act=="all":
        avgtbe=get_avg_tbe(tree)
        deeptbe=get_deepest_node_tbe(nodes)
        return(avgtbe,deeptbe)
    elif args.act=="avg":
        avgtbe=get_avg_tbe(tree)
        return(avgtbe)
    elif args.act=="deep":
        deeptbe=get_deepest_node_tbe(nodes)
        return(deeptbe)
    

    
    
if __name__=='__main__':
    
    import argparse
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("tree",type=str,help="The tree with the TBE support of another tree (computed using Booster).")
    app.add_argument("typ",type=str,choices=['tree','taxa'],help="Specify if it is a normal tree (in NH format) or a TOL.")
    app.add_argument("-act",type=str,choices=['avg','deep','all'],help="Action: calculate average TBE, deepest node's TBE or both")
    app.add_argument("-p",action='store_true',help="Print process info.")
    args=app.parse_args()
    
    if args.act=="all":
        avgtbe,deeptbe=count_tbe(args.tree,args.typ,args.act,args.p)
        print(str(round(avgtbe,6)) + " " + str(round(float(deeptbe),6)))
    elif args.act=="avg":
        tbe=count_tbe(args.tree,args.typ,args.act,args.p)
        print(str(round(tbe,6)))

    
    
