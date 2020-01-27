
__description__='''
    
    Given the sequences, PDB files (experimental and predicted), etc, merge the information to be used latter in PCA.
    
    Note that the data is organized by Pfam family.

'''


class Prot:
    def __init__(self,name):
        self.name=name
        self.seq="NA"
        self.seqlen="NA"
        self.nnatcon="NA"


def calc_stat(values,funct):
    
    import statistics as st
    if funct=="average":
        value=st.mean(values)
    elif funct=="median":
        value=st.median(values)
    return(value)
    

        
if __name__=='__main__':

    import argparse
    import collections
    import os
    from read_files import Seqs
    from read_files import read_multifasta_name_seq_len
    from extrac_info_functions import get_ncontacts, get_alnlen_rmsd_pid
    
    
    app = argparse.ArgumentParser(description=__description__)
    #app.add_argument('-dataset',type=str,choices=['OMA','BenchFam'],help="Dataset")
    app.add_argument('-pfam',type=str,help="Pfam family ID")
    app.add_argument('-taxa',type=str)
    app.add_argument('-seq',type=str,help="Multi-sequence in fasta format (.fa)")
    
    # Contacts
    app.add_argument('-natf',type=str,help="Folder with the native contact & contactbench files")
    app.add_argument('-threshold',type=float,default=0.35,help="Threshold for probability of contact")
    app.add_argument('-seqsep',type=int,default=8,help="Minimum separation in sequence")
    
    # PDB files
    app.add_argument('-exppdbf',type=str,help="Folder with the experimental structures. Note: the residues/atoms should start from 1.")
    
    # Function
    app.add_argument('-funct',type=str,choices=['average','median'],default='average',help="Function to be applied")
    app.add_argument('-out',type=str,help="Output file")
    
    args=app.parse_args()
    
    
    # Sequences
    seqs=read_multifasta_name_seq_len(args.seq)  # class Seqs ( name, seq, seqlen )
    
    # class Prot (name, seq, seqlen, nnatcon)
    prots=[]
    for seq in seqs:
        prot=Prot(seq.name)
        prot.seq=seq.seq
        prot.seqlen=seq.seqlen
        prot.nnatcon=get_ncontacts(args.natf+"/"+prot.name+".con",args.threshold,args.seqsep)
        prots.append(prot)
    
    # nseq & average/median seqlen,nnatcon
    seqlen=calc_stat([x.seqlen for x in prots],args.funct)
    nseq=len(prots)
    nnatcon=calc_stat([x.nnatcon for x in prots],args.funct)

    # RMSD, alnlen, seqid between pdb from the same Pfam family
    inalnlenl,inrmsdl,inseqidl=[],[],[]
    for prot1 in prots:
        for prot2 in prots:
            if prot1.name != prot2.name:
                pdb1=args.exppdbf+"/"+prot1.name+".pdb"
                pdb2=args.exppdbf+"/"+prot2.name+".pdb"
                inalnlen,inrmsd,inseqid=get_alnlen_rmsd_pid(pdb1,pdb2)
                inalnlenl.append(inalnlen)
                inrmsdl.append(inrmsd)
                inseqidl.append(inseqid)
    # Average/median
    inalnlen=calc_stat(inalnlenl,args.funct)
    inrmsd=calc_stat(inrmsdl,args.funct)
    inseqid=calc_stat(inseqidl,args.funct)
    
    # Write output
    if not os.path.isfile(args.out):
        o=open(args.out,'w')
        o.write("pfam taxa nseq seqlen nnatcon inrmsd inseqid inalnlen\n")
    else:
        o=open(args.out,'a')
    o.write("{} {} {} {} {} {} {} {}\n".format(args.pfam,args.taxa,nseq,int(round(seqlen,0)),int(round(nnatcon,0)),round(inrmsd,2),round(inseqid,2),int(round(inalnlen,0))))
    o.close()



