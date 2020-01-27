


__description__='''

    Given a set of trees, it computes the RF, the average TBE and the TBE value of the deepest node between each tree.
    
    As output, it creates a list (.txt) and a matrix (.mat) with the RF/TBE values. 

'''


def read_list(filename):
    ''' Read tree list'''
    
    import collections
    
    treelist=collections.OrderedDict()
    with open(filename) as l:
        for line in l:
            line=line.strip('\n')
            fields=line.split(" ")
            path=fields[0]
            name=fields[1]
            treelist[name]=path
    return(treelist)
    
    
def booster_tbe(treelist,outf):
    ''' Run booster '''
    
    import os
    from other.run_command import __run_command
    
    for name1,path1 in treelist.items():
        for name2,path2 in treelist.items():
            out="{}/{}_{}.nwk".format(outf,name1,name2)
            if not os.path.isfile(out): 
                command="booster -i {} -b {} -o {}".format(path1,path2,out)
                __run_command(" ",command," ")


def calc_tbe(treelist,outf,p):
    ''' Determine avg TBE and deepest node's TBE '''
    
    import numpy as np
    from count_tbe import count_tbe, get_avg_tbe, get_deepest_node_tbe
    
    # Initialize matrix
    n=len(treelist)
    shape=(n,n)
    avgmat=np.zeros(shape)
    deepmat=np.zeros(shape)

    # Determine TBE
    trees=treelist.keys()
    for i in range(n):
        for j in range(n):
            name1,name2=trees[i],trees[j]
            if name1=="tol":   ### TOL or tree
                typ="taxa"
            else:
                typ="tree"
            tbetree="{}/{}_{}.nwk".format(outf,name1,name2)
            avgtbe,deeptbe=count_tbe(tbetree,typ,p)
            avgmat[i,j]=avgtbe
            deepmat[i,j]=deeptbe
    # Determine symmeric TBE
    symavg,symdeep=sym_tbe(avgmat,deepmat,treelist)
    
    return(avgmat,deepmat,symavg,symdeep)   
    

def sym_tbe(avgmat,deepmat,treelist):
    ''' Calculate symmetric TBE '''
    
    import numpy as np
    
    # Initialize matrix
    n=len(treelist)
    shape=(n,n)
    symavg=np.zeros(shape)
    symdeep=np.zeros(shape)
    
    # Determine symmetric TBE
    for i in range(n):
        for j in range(n):
            symavg[i,j]=(avgmat[i,j]+avgmat[j,i])/2
            symdeep[i,j]=(deepmat[i,j]+deepmat[j,i])/2   # Be careful with this parameter. The deepest nodes might not be the same in both trees
            symavg[j,i],symdeep[i,j]=symavg[i,j],symdeep[i,j]
            
    return(symavg,symdeep)
    
    
def write_tbe(avgmat,deepmat,symavg,symdeep,treelist,outf):

    import numpy as np
    
    n=len(treelist)
    trees=treelist.keys()

    # List file
    lf=open(outf+"/tbe.txt","w")
    lf.write("str1 str2 avgtbe deeptbe sym_avgtbe sym_deeptbe\n")
    for i in range(n):
        for j in range(n):
            name1,name2=trees[i],trees[j]
            lf.write("{}_{} {} {} {} {}\n".format(name1,name2,round(avgmat[i,j],4),round(deepmat[i,j],4),round(symavg[i,j],4),round(symdeep[i,j],4)))
    lf.close()    
    
    # Write matrix to file
    tmp=np.around(avgmat,decimals=4)
    np.savetxt(outf+"/avgtbe.mat",tmp,delimiter=" ",fmt='%.4f')
    tmp=np.around(deepmat,decimals=4)
    np.savetxt(outf+"/deeptbe.mat",tmp,delimiter=" ",fmt='%.4f') 
    tmp=np.around(symavg,decimals=4)
    np.savetxt(outf+"/symavgtbe.mat",tmp,delimiter=" ",fmt='%.4f') 
    tmp=np.around(symdeep,decimals=4)
    np.savetxt(outf+"/symdeeptbe.mat",tmp,delimiter=" ",fmt='%.4f') 
    
          
def calc_rf(treelist,outf):
    ''' Determine RF score '''
    
    import os
    from other.run_command import __run_command
    
    # Write trees.nwk
    trees=open(outf+"/trees.nwk","w")
    for path in treelist.values():
        with open(path) as f:
            for line in f:
                line=line.strip("\n")
                trees.write(line+"\n")
    
    # Calculate RF --> rf.txt
    script=os.path.dirname(os.path.abspath(__file__))+"/rf.R" 
    command="Rscript {} {}/treelist {}/trees.nwk {}/rf.txt {}/rf.xlsx".format(script,outf,outf,outf,outf)
    __run_command(" ",command," ")
    


if __name__=='__main__':
    
    import argparse
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument('-i',type=str,required=True,help='A list with the absolute path and name of the trees (in .nwk format).')  # Including exp1d, exp3d, deep1d, deep3d, nat1d, nat3d
    app.add_argument('-o',type=str,required=True,help='Output folder.')
    args=app.parse_args()
    
    # Read tree
    treelist=read_list(args.i)
    # TBE
    booster_tbe(treelist,args.o)
    avgmat,deepmat,symavg,symdeep=calc_tbe(treelist,args.o,None)
    write_tbe(avgmat,deepmat,symavg,symdeep,treelist,args.o)
    # RF
    calc_rf(treelist,args.o)
    
