

__description__='''

Reorganize the output from TREE-PUZZLE and Scaleboot into one matrix.

'''

class Puzzle:
    def __init__(self,tree,logl,diff,kh1s,sh,kh2s):
        self.tree=tree
        self.logl=logl
        self.diff=diff
        self.kh1s=kh1s
        self.sh=sh
        self.kh2s=kh2s

class Scaleboot:
    def __init__(self,tree,bp,au,si,hypothesis):
        self.tree=tree
        self.bp=bp
        self.au=au
        self.si=si
        self.hypothesis=hypothesis

def read_puzzle(filename,ntrees):
    
    # Get lines
    lines=[]
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            if line:
                lines.append(line)
                
    # Get block
    n=lines.index("COMPARISON OF USER TREES (NO CLOCK)")
    puzzle=[]
    for i in range(ntrees):
        nline=n+3+i
        line=lines[nline]
        tree=line[1]
        logl=line[5:13]
        diff=line[16:22]
        kh1s=line[38:46]
        sh=line[49:57]
        kh2s=line[72:76]
        puzzle.append(Puzzle(tree,logl,diff,kh1s,sh,kh2s))
        
    return(puzzle)


def read_scaleboot(filename):
    
    f=open(filename,"r")
    lines=f.readlines()
    scaleboot=[]
    for i in range(4,len(lines)):
        line=lines[i].strip("\n")
        fields=line.split(" ")
        fields=[x for x in fields if x]
        tree=fields[0][1]
        bp=fields[1]+" "+fields[2]
        au=fields[3]+" "+fields[4]
        si=fields[5]+" "+fields[6]
        hypothesis=fields[11]
        scaleboot.append(Scaleboot(tree,bp,au,si,hypothesis))
    return(scaleboot)



if __name__ == '__main__':
    
    import argparse
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("-puzzle",type=str,help="Report from TREE-PUZZLE.")
    app.add_argument("-scaleboot",type=str,help="Output from Scaleboot")
    app.add_argument("-out",type=str,help="Output file.")
    app.add_argument("-ntrees",type=int,help="Number of trees being compared.")
    args=app.parse_args()
    
    puzzle=read_puzzle(args.puzzle,args.ntrees)
    scaleboot=read_scaleboot(args.scaleboot)
    
    out=open(args.out,"w")
    out.write("log L;difference;1sKH;SH;2sKH;BP;AU;SI;hypothesis\n")
    for tree in range(len(puzzle)):
        puzz=puzzle[tree]
        scal=scaleboot[tree]
        out.write(puzz.logl+";"+puzz.diff+";"+puzz.kh1s+";"+puzz.sh+";"+puzz.kh2s+";"+scal.bp+";"+scal.au+";"+scal.si+";"+scal.hypothesis+"\n")
    out.close()
    
