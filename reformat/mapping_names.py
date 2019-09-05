
__description__='''

    Change sequence name based on PDB id to UniProt id in a fasta file, and viceversa.
    
    
    The PDB id is used as filename for some programs like DeepConPred2 or Spider3,
    whereas the UniProt id is required for aligners like TMalign.

'''


def __map_names(fasta_file,reference_file,corrected_fasta_file,to):
    ''' Change pdb id to uniprot id, and viceversa.
    
    to = "pdb_to_uniprot" / "uniprot_to_pdb" / ""
    '''
    
    names = {}
    fasta = {}
    out = open(corrected_fasta_file,"wt")
  
  
    # Get fasta lines  {name : seq}
    with open(fasta_file) as fa:
        for line in fa:
            line = line.strip()
            if line[0] == ">":
                ID = line[1:]
                fasta[ID] = ""
            else:
                fasta[ID] += line
  
  
    # Read reference file
    uniprot = []
    pdb = []
    with open(reference_file) as ref:
        for line in ref:
            line = line.strip()
            line = line.split(" ")
            uniprot_id = line[0][1:]
            pdb_id = line[2]
            uniprot.append(uniprot_id)
            pdb.append(pdb_id)
            
    # Check change
    if to != "pdb_to_uniprot" and to != "uniprot_to_pdb":
        if fasta.keys()[0] in uniprot:
            to = "uniprot_to_pdb"
        elif fasta.keys()[0] in pdb:
            to = "pdb_to_uniprot"
    # Get names {original_id : new_id} 
    for i in range(len(uniprot)):
        if to == "pdb_to_uniprot":
            names[pdb[i]] = uniprot[i]
        elif to == "uniprot_to_pdb":
            names[uniprot[i]] = pdb[i]
            

    # Write output
    for original_id in fasta.iterkeys():
        new_id = names[original_id]
        seq = fasta[original_id]
        out.write(">{}\n{}\n".format(new_id,seq))
      
    out.close()
        
  


if __name__ == "__main__":
  
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    # I/O
    parser.add_argument("fasta_file",type=str,help="Fasta file. PDB id is used as sequence name.")
    parser.add_argument("reference_file",type=str,help="Reference file. Each line should contain the UniProt identifier, a separator and the PDB identifier. For example: >AMCY_PARDE/43-131 _P_ 1MDAA-1")
    parser.add_argument("out",type=str,help="Output fasta file with corrected sequence name.")
    parser.add_argument("change",type=str,default=None,nargs='?',help="It can be 'pdb_to_uniprot' or 'uniprot_to_pdb'. Default=None, so it will read the fasta file names and change it automatically.")
    args = parser.parse_args()
  
    __map_names(args.fasta_file,args.reference_file,args.out,args.change)
