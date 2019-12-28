

__description__='''

Given a tree in newick format, this script removes the edge lenths and bootstrap supports while keeping the leaf names.

'''

def read_tree(tree_file):
    with open(tree_file) as tf:
        for line in tf:
            tree=line.rstrip()
    return(tree)
    
def trim_tree(tree):
    import re
    tree=re.sub(':-?[0-9]+\.[0-9]+','',tree)  # Remove edge length
    for m in re.finditer('\)([0-9]+)',tree):
        tree=tree[:m.start(1)]+'$'*(m.end(1)-m.start(1))+tree[m.end(1):]  # Remove bootstrap support
    tree=tree.replace('$','')
    return(tree)

if __name__=='__main__':

    import argparse
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("intree",type=str,help="Input tree (NH format).")
    app.add_argument("outtree",type=str,help="Output tree (NH format).")
    args=app.parse_args()
    
    tree1=read_tree(args.intree)
    tree2=trim_tree(tree1)
    
    out=open(args.outtree,"w")
    out.write(tree2+"\n")
    out.close
    