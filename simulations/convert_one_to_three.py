#
# Convert seq from "one-letter" to "three-letter" code
# Usage: python path_to/convert_one_to_three.py <fasta.fa> <outputfile.seq>
#


def __convert_one_to_three(inputfilename,outputfilename):
  
    from Bio.PDB.Polypeptide import *
    from Bio import SeqIO  

    myfas = SeqIO.read(inputfilename, "fasta")
    f = open(outputfilename,"wt")

    for i in myfas.seq:
        f.write("{} ".format(one_to_three(i)))

    f.close()



if __name__ == '__main__':
    
    import sys
    inputfilename=sys.argv[1]
    outputfilename=sys.argv[2]
      
    __convert_one_to_three(inputfilename,outputfilename)
