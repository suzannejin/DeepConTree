
__description__='''

It counts the number of common finds (F) of two given TOL compliance output files.

'''

def read_compliances(filename):
    compl=[]
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            fields=line.split(" ")
            compl.append(fields[6])
    return(compl)

def compare_compliances(compl):
    maxf=min(x.count("F") for x in compl)
    maxm=min(x.count("M") for x in compl)
    f,m=0,0
    for i in range(len(compl[0])):
        f2,m2=0,0
        for c in range(1,len(compl)):
            if compl[0][i]==compl[c][i]:
                if compl[0][i]=="F":
                    f2+=1
                elif compl[0][i]=="M":
                    m2+=1
        if f2==len(compl)-1:
            f+=1
        if m2==len(compl)-1:
            m+=1
    return(f,m,maxf,maxm)


if __name__=='__main__':
    
    import argparse
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-i",type=str,nargs='+',help="Input files.")
    args=app.parse_args()
    
    compliances=[]
    for f in args.i:
        compliances.append(read_compliances(f))
    f,m,maxf,maxm=compare_compliances(compliances)
    print("{}/{} {}/{}".format(f,maxf,m,maxm))
    
    