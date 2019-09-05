
#
# Prepare restraint .tbl files
# 

def __prepare_restraints_contacts(rr_file,fasta_file,out,restraint):

    ''' Prepare contact restraints from RR file and b-strand O-O restraints from SSE file. 
    
    rr_file:    Contact file in RR format
    fasta_file: Sequence file in fasta format.
    out:        Output filename.
    restraint:  List. "o-o" for beta-strand O-O constraints.
    
    '''


    import random
    import re
    from Bio import SeqIO


    myfas = SeqIO.read(fasta_file, "fasta")
    f = open(out,"wt")
    f.write("set echo=false end\nset wrnlev=0 end\n")



    # Calculate mean weight
    x = []
    N = 0
    sum = 0
    with open (rr_file) as file:
        for l in file:
            x.append(l.strip())
            tmp_list = l.split()
            sum += float(tmp_list[4])
            N += 1
            
    mean = sum / N


    # Write contact restraints
    f.write("!residue contacts\n")
    for con in x:
        tmp_list = con.split()
        atom1 = 'CB'
        atom2 = 'CB'
        first_res = int(tmp_list[0])
        sec_res = int(tmp_list[1])
        prob = float(tmp_list[4])
        if ( myfas.seq[first_res-1] == 'G' ):
            atom1 = 'CA'
        if ( myfas.seq[sec_res-1] == 'G' ):
            atom2 = 'CA' 
        f.write("WEIGht {}\n".format(round(prob/mean,4)))
        f.write("ASSIgn (resid\t {}  and name  {} \t)(resid\t {}  and name  {} \t) 8.0 8.0 0.0\n".format(tmp_list[0],atom1,tmp_list[1],atom2))


    # Write b-strand O-O restraints
    if "o-o" in restraint:
        f.write("!beta strand O-O\n")
        sse=""
        for line in sselines[1:]:
            fields=line.split()
            ss=fields[2]
            sse+=ss
        ssebeta=[]
        for i in re.finditer("E+",sse):
            if "E" in i.group(0):
                ssebeta.append(i.span())
        for element in ssebeta:
            start,end=element[0]+1,element[1]+1
            if end-start < 3:
                continue
            for n in range(start,end-1):
                f.write("WEIGht 1.0000\nASSIgn (resid\t {}  and name  O \t)(resid\t {}  and name  O \t) 4.6 0.1 0.1\n".format(n,n+1))


    f.write("set echo=true end\nset wrnlev=5 end\n")    
    f.close()




def __prepare_restraints_dihedral(sse_file, out):
  
    ''' Prepare dihedral angles restraints from SSE file. 
    
    sse_file:   Prediction of secondary structure elements in Spider3-local format.
    out:        Output filename.
    
    '''
    
    import re
    

    ssefile=open(sse_file,'rt')
    outputfile=open(out,'wt')

    sselines=ssefile.readlines()
    ssefile.close()


    # Secondary structure elements (SSE)
    seq=""
    sse=""
    for line in sselines[1:]:
        line=line.strip()
        fields=line.split()
      
        res,ss=fields[1],fields[2]
        seq+=res
        sse+=ss
        
    # Only alpha-helix and beta-strands/bridges
    sseAlphaBeta={}
    for i in re.finditer("(H+)|(E+)",sse):
        start,end=i.span()[0],i.span()[1]
        if "H" in i.group(0):
            if end-start >= 4:
              sseAlphaBeta[i.span()]="H"
        else:
            if end-start >= 3:
              sseAlphaBeta[i.span()]="E"


    outputfile.write("set echo=false end\nset wrnlev=0 end\n")
    

    # Phi angles
    outputfile.write("!phi\n")
    for element in sorted(sseAlphaBeta):
        ss=sseAlphaBeta[element]
        start,end=element[0],element[1]
        if ss=="H":
            constant,angle,nrange,ed=1.0,-63.5,4.5,2
        else:
            constant,angle,nrange,ed=1.0,-118.0,10.7,2
        for n in range(start,end-2):
            pos=n+1
            outputfile.write("ASSIgn (resid {} and name c) (resid {} and name n) (resid {} and name ca) (resid {} and name c) {} {} {} {}\n".format(pos,pos+1,pos+1,pos+1,constant,angle,nrange,ed))
        
    # Psi angles
    outputfile.write("!psi\n")
    for element in sorted(sseAlphaBeta):
        ss=sseAlphaBeta[element]
        start,end=element[0],element[1]
        if ss=="H":
            constant,angle,nrange,ed=1.0,-41.5,5,2
        else:
            constant,angle,nrange,ed=1.0,134,8.6,2
        for n in range(start+1,end-1):
            pos=n+1
            outputfile.write("ASSIgn (resid {} and name n) (resid {} and name ca) (resid {} and name c) (resid {} and name n) {} {} {} {}\n".format(pos,pos,pos,pos+1,constant,angle,nrange,ed))

    outputfile.close()


