

class Contact_metrics:
    def __init__(self,accuracy,precision,recall,f1,ncontacts,threshold,seqsep):
        self.accuracy=accuracy
        self.precision=precision
        self.recall=recall
        self.f1=f1
        self.ncontacts=ncontacts
        self.threshold=threshold
        self.seqsep=seqsep
        

def get_contactbench_info(contactbench,threshold,seqsep):
    '''
    Input
    contactbench : class Contactbench (i, j, dist, prob, true)
    threshold    : a contact is considered if prob >= threshold
    seqsep       : minimun separation in sequence
    
    Output
    metrics      : class Contact_metrics (accuracy, precision, recall, f1, ncontacts)
    '''
    for con in contactbench:
        tp,fp,tn,fn=0.0,0.0,0.0,0.0
        ncon=0
        if con.dist>=seqsep:
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
    accuracy=(tp+tn)/(tp+tn+fp+fn)
    precision=tp/(tp+fp)
    recall=tp/(tp+fn)
    f1=2*tp/(2*tp+fp+fn)
    metrics=Contact_metrics(accuracy,precision,recall,f1,ncon,threshold,seqsep)
    return(metrics)


def get_ncontacts(contactfile,threshold,seqsep):
    '''
    Input
    contactfile  : contact filename
    threshold    : a contact is considered if prob >= threshold
    seqsep       : minimun separation in sequence
    
    Output
    ncon         : number of contacts 
    '''
    from read_files import read_contacts
    cons=read_contacts(contactfile)
    ncon=0
    for con in cons:
        if con.prob >= threshold and abs(con.i - con.j) >= seqsep:
            ncon+=1
    return(ncon)
    

def get_alnlen_rmsd_pid(pdb1,pdb2):
    import re
    import subprocess
    out=subprocess.Popen(['TMalign',pdb1,pdb2],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout,stderr=out.communicate()
    for line in stdout.split("\n"):
        line=line.strip("\n")
        m=re.match('^Aligned\s+length=\s+(\d+),\s*RMSD=\s*(\d+\.\d+),\s*Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)',line)
        if m:
            inalnlen=int(m.group(1))
            inrmsd=float(m.group(2))
            inseqid=float(m.group(3))
            return(inalnlen,inrmsd,inseqid)
            
 
 
