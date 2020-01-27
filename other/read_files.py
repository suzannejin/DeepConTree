
class Seqs:
    def __init__(self,name,seq,seqlen):
        self.name=name
        self.seq=seq
        self.seqlen=seqlen
        

class Contacts:
    def __init__(self,i,j,d1,d2,prob):
        self.i=i
        self.j=j
        self.d1=d1
        self.d2=d2
        self.prob=prob
        

class Contactbench:
    def __init__(self,i,j,dist,prob,true):
        self.i=i
        self.j=j
        self.dist=dist
        self.prob=prob
        self.true=true  # 0 for FALSE, 1 for TRUE



def read_multifasta(filename):
    ''' 
    Input  - filename : multifasta file
    Output - seqs     : dictionary {sequence name : sequence}
    '''
    import collections
    seqs=collections.OrderedDict()
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            if line[0]==">":
                name=line[1:]
                seqs[name]=""
            else:
                seqs[name]+=line
    return(seqs)


def read_multifasta_name_seq_len(filename):
    import collections
    fasta=read_multifasta(filename)
    seqs=[]
    for name,seq in fasta.items():
        seqs.append(Seqs(name,seq,len(seq)))
    return(seqs)


def read_contacts(filename):
    '''
    Input  - filename : contact file
    Output - cons     : a list of class Contacts (i,j,d1,d2,prob)
    '''
    cons=[]
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            fields=line.split(" ")
            i,j=int(fields[0]),int(fields[1])
            d1,d2=int(fields[2]),int(fields[3])
            prob=float(fields[4])
            con=Contacts(i,j,d1,d2,prob)
            cons.append(con)
    return(cons)            


def read_contactbench(filename):
    '''
    Input  - filename : contactbench file
    Output - cons     : a list of class Contactbench (i,j,dist,prob,true)
    '''
    cons=[]
    with open(filename) as f:
        for line in f:
            line=line.strip("\n")
            fields=line.split(" ")
            i,j=int(fields[0]),int(fields[1])
            dist=int(fields[2])
            prob=float(fields[3])
            if fields[4]=="TRUE" or fields[4]=="1":
                true=1
            else:
                true=0
            con=Contactbench(i,j,dist,prob,true)
            cons.append(con)
    return(cons)
            
