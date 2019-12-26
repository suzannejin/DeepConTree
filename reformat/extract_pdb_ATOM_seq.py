

def pdb2seq(pdbname):
    import Bio.PDB.Polypeptide as bio
    seq = ""
    with open(pdbname,"r") as pdb:
       prev_n,n = 0,0
       for line in pdb:
           line = line.strip("\n")
           if line[:4] == "ATOM":
               n = int(line[23:26])
               if n != prev_n:
                   aa = line[17:20]
                   seq += bio.three_to_one(aa)
                   prev_n = n
    return(seq)


def writeATOMseq2fasta(pdbname,outseq):
    import os
    seq = pdb2seq(pdbname)
    out = open(outseq,"w")
    seqname = os.path.basename(pdbname).replace(".pdb","")
    out.write(">" + seqname + "\n" + seq + "\n")
    out.close()


if __name__ == '__main__':
    import argparse
    app = argparse.ArgumentParser(description="Extract the sequence from the ATOM section of a PDB file, and write it to a fasta file.")
    app.add_argument("pdb",type=str,help="PDB file.")
    app.add_argument("out",type=str,help="Output sequence file in fasta format (.fa)")
    args = app.parse_args()
    writeATOMseq2fasta(args.pdb,args.out)
    
    