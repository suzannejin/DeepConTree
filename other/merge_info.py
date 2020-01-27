
__description__='''
    
    Given the sequences, PDB files (experimental and predicted), etc, merge the information to be used latter in PCA.
    
    Note that the data is organized by Pfam family.

'''


class Info:
    
    def __init__(self,pfam,struct):
        
        self.pfam=pfam
        self.struct=struct     # Type of structure : exp, deep, nat
        
        self.nseq="NA"         # Number sequences
        self.seqlen="NA"       # Average/median sequence length
        
        self.ncon="NA"         # Average/median number contacts
        self.accuracy="NA"     # Average/median accuracy contacts
        self.precision="NA"    # Average/median precision contacts
        self.recall="NA"       # Average/median recall contacts
        self.f1="NA"           # Average/median F1-score contacts
        
        self.rmsd="NA"         # Average/median RMSD structures compared to exp
        self.tm="NA"           # Average/median TM-score structures compared to exp
        
        self.inrmsd="NA"       # Average/median pairwise RMSD between structures between proteins of the same Pfam
        self.inseqid="NA"      # Average/median pairwise sequence id between proteins of the same Pfam
        self.inalnlen="NA"     # Pairwise aligned sequence length between proteins of the same Pfam
        
        self.tc="NA"           # Average/median TC score compared to exp
        self.sp="NA"           # Average/median SP score compared to exp
        
        self.nbranch="NA"      # Number branches in the tree
        self.tolcompl="NA"     # TOL compliance
        self.toltbe="NA"       # Symmetric average TBE - TOL
        self.exptbe="NA"       # Symmetric average TBE - exp
        self.nsexptbe="NA"     # Non-symmetric average TBE - exp (as ref) - tree (as replicate)
    
        
class Files:
    def __init__(self,struct,seqs,msa,conf,pdbf,tree):
        self.struct=struct
        self.seqs=seqs
        self.msa=msa
        self.conf=conf
        self.pdbf=pdbf
        self.tree=tree
        

def get_seqlen_nseq(seqs,funct):
    '''
    Input
    seqs    : dictionary { sequence name : sequence }, obtained from a given multifasta file
    
    Output 
    nseq    : number of sequences
    avglen  : average length of the sequences
    '''
    import statistics as st
    lens=[]
    for seq in seqs.values():
        lens.append(len(seq))
    n=len(lens)
    if funct=="average":
        lens=st.mean(lens)
    elif funct=="median":
        lens=st.median(lens)
    return(lens,n)
            

def get_pred_contacts_info(foldername,threshold,seqsep,funct):
    '''
    Input
    foldernane : folder with the contacts (.con) and contactbench (.contactbench) files
    
    Output
    avgacc     : average accuracy
    avgprec    : average precision
    avgrec     : average recall
    avgf1      : average F1-score
    '''
    
    import glob
    import statistics as st
    from read_files import read_contactbench
    
    contactbenchs=glob.glob(foldername+"/*.contactbench")
    accuracy,precision,recall,f1,ncons=[],[],[],[],[]
    for filename in sorted(contactbenchs):
        cons=read_contactbench(filename)
        tp,fp,tn,fn=0.0,0.0,0.0,0.0
        ncon=0
        for con in cons:
            if con.dist<seqsep:
                continue
            if con.prob>=threshold:
                ncon+=1
                if con.true==1:
                    tp+=1.0
                else:
                    fp+=1.0
            else:
                if con.true==0:
                    tn+=1.0
                else:
                    fn+=1.0
        accuracy.append((tp+tn)/(tp+tn+fp+fn))
        precision.append(tp/(tp+fp))
        recall.append(tp/(tp+fn))
        f1.append(2*tp/(2*tp+fp+fn))
        ncons.append(ncon)
    if funct=="average":
        accuracy=st.mean(accuracy)
        precision=st.mean(precision)
        recall=st.mean(recall)
        f1=st.mean(f1)
        ncons=st.mean(ncons)
    elif funct=="median":
        accuracy=st.median(accuracy)
        precision=st.median(precision)
        recall=st.median(recall)
        f1=st.median(f1)
        ncons=st.median(ncons)
    return(ncons,accuracy,precision,recall,f1)


def get_nat_contacts(foldername,threshold,seqsep,funct):
    
    import glob
    import statistics as st
    from read_files import read_contacts
    
    nnats=[]
    contacts=glob.glob(foldername+"/*.con")
    for filename in sorted(contacts):
        cons=read_contacts(filename)
        nnat=0
        for con in cons:
            if abs(con.i - con.j)>=8:
                nnat+=1
        nnats.append(nnat)
    if funct=="average":
        nnats=st.mean(nnats)
    elif funct=="median":
        nnats=st.median(nnats)
    return(nnats)
    

def get_rmsd_tm(seqs,folder1,folder2,funct):
    
    import re
    import statistics as st
    import subprocess
    
    rmsd,tm=[],[]
    for prot in seqs:
        pdb1=folder1+"/"+prot+".pdb"
        pdb2=folder2+"/"+prot+".pdb"
        out=subprocess.Popen(['TMscore',pdb1,pdb2],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stdout,stderr=out.communicate()
        for line in stdout.split("\n"):
            line=line.strip("\n")
            mrmsd=re.match('^RMSD\s+of\s+the\s+common\s+residues=\s*(\d+\.\d+)',line)
            mtm=re.match('^TM-score\s*=\s(\d\.\d+)',line)
            if mrmsd:
                rmsd.append(float(mrmsd.group(1)))
            elif mtm:
                tm.append(float(mtm.group(1)))
                break
    if funct=="average":
        rmsd=st.mean(rmsd)
        tm=st.mean(tm)
    elif funct=="median":
        rmsd=st.median(rmsd)
        tm=st.median(tm)
    return(rmsd,tm)


def get_alnlen_rmsd_pid(seqs,folder,funct):
    
    import re
    import statistics as st
    import subprocess
    
    alnlen,rmsd,seqid=[],[],[]
    for prot1 in seqs:
        for prot2 in seqs:
            if prot1==prot2:
                continue
            pdb1=folder+"/"+prot1+".pdb"
            pdb2=folder+"/"+prot2+".pdb"
            out=subprocess.Popen(['TMalign',pdb1,pdb2],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            stdout,stderr=out.communicate()
            for line in stdout.split("\n"):
                line=line.strip("\n")
                m=re.match('^Aligned\s+length=\s+(\d+),\s*RMSD=\s*(\d+\.\d+),\s*Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)',line)
                if m:
                    alnlen.append(int(m.group(1)))
                    rmsd.append(float(m.group(2)))
                    seqid.append(float(m.group(3)))
                    break
    if funct=="average":
        alnlen=st.mean(alnlen)
        rmsd=st.mean(rmsd)
        seqid=st.mean(seqid)
    elif funct=="median":
        alnlen=st.median(alnlen)
        rmsd=st.median(rmsd)
        seqid=st.median(seqid)
    return(alnlen,rmsd,seqid)
    

def get_tc_sp(msa1,msa2):
    
    import re
    import subprocess
    
    tmp=[]
    for comp in ['sp','tc']:
        out=subprocess.Popen(['t_coffee','-other_pg','aln_compare','-al1',msa1,'-al2',msa2,'-compare_mode',comp],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stdout,stderr=out.communicate()
        lines=filter(None,stdout.split("\n"))
        m=re.match('.+\s+(\d+\.\d+)\s+(\d+\.\d+)',lines[2])
        tmp.append(m.group(2))
    sp,tc=float(tmp[0]),float(tmp[1])
    return(sp,tc)

    
def get_tbe(tbefolder,struct):
    
    from trees.count_tbe import get_avg_tbe, get_deepest_node_tbe, count_tbe
    
    toltrees=[tbefolder+"/tol_"+struct+"3d.nwk",tbefolder+"/"+struct+"3d_tol.nwk"]
    exptrees=[tbefolder+"/exp3d_"+struct+"3d.nwk",tbefolder+"/"+struct+"3d_exp3d.nwk"]
    
    toltbe,exptbe=[],[]
    for tree in toltrees:
        avgtbe,deeptbe=count_tbe(tree,"taxa","")
        toltbe.append(avgtbe)
    for tree in exptrees:
        avgtbe,deeptbe=count_tbe(tree,"tree","")
        exptbe.append(avgtbe)
    nsexptbe=exptbe[0]
    toltbe=(toltbe[0]+toltbe[1])/2
    exptbe=(exptbe[0]+exptbe[1])/2
    return(nsexptbe,toltbe,exptbe)
    

def write_all_info(info,out):
    
    from os import path
    
    if not path.exists(out):
        o=open(out,'w')
        o.write("pfam struct nseq seqlen ncontacts accuracy precision recall f1 rmsd tm inrmsd inseqid inalnlen tc sp nbranch tolcompl toltbe exptbe nsexptbe\n")
    else:
        o=open(out,'a')
    o.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(info.pfam,info.struct,info.nseq,int(info.seqlen),int(info.ncon),round(info.accuracy,2),round(info.precision,2),round(info.recall,2),round(info.f1,2),round(info.rmsd,2),round(info.tm,2),round(info.inrmsd,2),round(info.inseqid,2),int(info.inalnlen),round(info.tc,2),round(info.sp,2),info.nbranch,info.tolcompl,round(info.toltbe,2),round(info.exptbe,2),round(info.nsexptbe,2)))
    o.close()


def write_basic(info,out):
    
    from os import path
    
    if not path.exists(out):
        o=open(out,'w')
        o.write("pfam nseq seqlen\n")
    else:
        o=open(out,'a')
    o.write("{} {} {}\n".format(info.pfam,info.nseq,int(info.seqlen)))
    o.close()
    exit(1)
    

        
if __name__=='__main__':

    import argparse
    import collections
    from read_files import read_multifasta
    
    app = argparse.ArgumentParser(description=__description__)
    #app.add_argument('-dataset',type=str,choices=['OMA','BenchFam'],help="Dataset")
    app.add_argument('-pfam',type=str,help="Pfam family ID")
    app.add_argument('-seq',type=str,help="Multi-sequence in fasta format (.fa)")
    
    # Contacts
    app.add_argument('-conf',type=str,default=None,help="Folder with the predicted contact & contactbench files")
    app.add_argument('-natf',type=str,default=None,help="Folder with the native contact & contactbench files")
    app.add_argument('-threshold',type=float,default=0.35,help="Threshold for probability of contact")
    app.add_argument('-seqsep',type=int,default=8,help="Minimum separation in sequence")
    
    # PDB files
    app.add_argument('-exppdbf',type=str,default=None,help="Folder with the experimental structures. Note: the residues/atoms should start from 1.")
    app.add_argument('-deeppdbf',type=str,default=None,help="Folder with the predicted structures")
    app.add_argument('-natpdbf',type=str,default=None,help="Folder with the native structures")
    
    # Alignments 
    app.add_argument('-expmsa',type=str,default=None,help="MSA - exp")
    app.add_argument('-deepmsa',type=str,default=None,help="MSA - deep")
    app.add_argument('-natmsa',type=str,default=None,help="MSA - nat")
    
    # Trees
    app.add_argument('-exptree',type=str,default=None,help="Phylo3D-exp tree")
    app.add_argument('-deeptree',type=str,default=None,help="Phylo3D-deep tree")
    app.add_argument('-nattree',type=str,default=None,help="Phylo3D-nat tree")
    app.add_argument('-tbef',type=str,default=None,help="Folder with all TBE outputs")
    
    # Function
    app.add_argument('-funct',type=str,choices=['average','median'],default='average',help="Function to be applied")
    app.add_argument('-out',type=str,help="Output file")
    app.add_argument('-basic',action='store_true',help="If basic, then just compute the seqlen, nseq, ntolbranch, nnatcon")
    
    args=app.parse_args()
    
    
    # Sequences
    seqs=read_multifasta(args.seq)  # { sequence name : sequence }
    seqlen,nseq=get_seqlen_nseq(seqs,args.funct)
            
    
    # Declare files
    files=collections.OrderedDict()
    files["exp"]=Files("exp",seqs,args.expmsa,None,args.exppdbf,args.exptree)
    files["deep"]=Files("deep",seqs,args.deepmsa,args.conf,args.deeppdbf,args.deeptree)
    files["nat"]=Files("nat",seqs,args.natmsa,args.natf,args.natpdbf,args.nattree)
    
    
    for struct,a in files.items():
        info=Info(args.pfam,struct)
        info.nseq=nseq
        info.seqlen=seqlen
        
        # Only write basic information: Pfam, sequence length, number sequences
        if args.basic:
            write_basic(info,args.out)
            
        
        if struct=="nat" and args.nattree=="NA":   # Skip when there are no nat-structures
            continue
        
        # Contacts
        if struct=="deep":
            ncon,accuracy,precision,recall,f1=get_pred_contacts_info(a.conf,args.threshold,args.seqsep,args.funct)
            info.ncon=ncon
            info.accuracy=accuracy
            info.precision=precision
            info.recall=recall
            info.f1=f1
        else:
            info.ncon=get_nat_contacts(args.natf,args.threshold,args.seqsep,args.funct)
            info.accuracy=1.00
            info.precision=1.00
            info.recall=1.00
            info.f1=1.00
        
        # Structures
        # RMSD & TM-score compared to exp
        rmsd,tm=get_rmsd_tm(seqs,args.exppdbf,a.pdbf,args.funct)
        info.rmsd=rmsd
        info.tm=tm
        
        # Pairwise RMSD & seqid within Pfam
        alnlen,rmsd,seqid=get_alnlen_rmsd_pid(seqs,a.pdbf,args.funct)
        info.inalnlen=alnlen
        info.inrmsd=rmsd
        info.inseqid=seqid
        
        # MSA compared to exp
        sp,tc=get_tc_sp(args.expmsa,a.msa)
        info.sp=sp
        info.tc=tc
        
        # Trees
        nsexptbe,toltbe,exptbe=get_tbe(args.tbef,struct)
        info.nsexptbe=nsexptbe
        info.toltbe=toltbe
        info.exptbe=exptbe
        
        # Write info
        write_all_info(info,args.out)
        
        
    
    
