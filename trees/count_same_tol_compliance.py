
__description__='''

It counts the number of common finds (F) of two given TOL compliance output files.

'''

class Compl:
    def __init__(self):
        self.f=[]
        self.s=[]
        self.commonF=0
        self.commonM=0


def read_compliances(filename):
    compl=Compl()
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            fields=line.split(" ")
            compl.f.append(fields[6])
            compl.s.append(fields[7])
    return(compl)

def check_compliances(compl):
    ''' Remove the last branch from 1D-trees '''
    nbranches=min(len(x.f) for x in compl)  
    for i in range(len(compl)):
        if len(compl[i].f)==nbranches+1:
            compl[i].f=compl[i].f[:-1]
        if len(compl[i].s)==nbranches+1:
            compl[i].s=compl[i].s[:-1]
    return(compl)
            
def compare_compliances(compl):
    ''' Common TOL-compliance between tree1 and tree2 '''
    ffs=[x.f.count("F") for x in compl]
    mms=[x.f.count("M") for x in compl]
    maxf=min(ffs)
    maxm=min(mms)
    fs,ms=[],[]
    for c in compl:
        f,m=0,0
        for i in range(len(c.f)):
            if c.s[i]=="S":
                if c.f[i]=="F":
                    f+=1
                elif c.f[i]=="M":
                    m+=1
        fs.append(f)
        ms.append(m)
    f=min(fs)
    m=min(ms)
    return(f,m,maxf,maxm,ffs,mms)


if __name__=='__main__':
    
    import argparse
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-i",type=str,nargs='+',help="Input files.")  # Give two input files
    args=app.parse_args()
    
    # Read files
    compliances=[]
    for f in args.i:
        compliances.append(read_compliances(f))
    compliances=check_compliances(compliances)
    
    # Determine common TOL-compliances
    f,m,maxf,maxm,ffs,mms=compare_compliances(compliances)
    print("{}/{}/{}/{} {}/{}/{}/{}".format(f,maxf,ffs[0],ffs[1],m,maxm,mms[0],mms[1]))
    
    