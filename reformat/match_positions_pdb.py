


def __match_pdb(pdb,output):
    ''' It corrects the atom/residue positions of a pdb file.
    
    In this way, the first atom/residue that appears in the given pdb file would have position 1.
    '''
  
    import re
    
    out = open(output,"wt")
    
    n = 0
    with open(pdb) as f:
        for line in f:
            line=line.strip("\n")
            
            # ATOM lines
            if line[0:4] == "ATOM":
                matches = re.finditer(r' ([0-9]+) ',line)  
                m=0
                for i in matches:
                    if m==0:   # Get atom number (atom_value_ori), the starting position in the line, the ending position in the line
                        atom_start=i.span()[0]
                        atom_end=i.span()[1]
                        atom_value_ori=int(i.group()[1:-1])
                    if m==1:  # Get residue number (res_value_ori), the starting position in the line, the ending position in the line
                        res_start=i.span()[0]
                        res_end=i.span()[1]
                        res_value_ori=int(i.group()[1:-1])
                    m+=1
                    
                # Get min atom/residue number from the first ATOM line
                if n==0:
                    minatom=atom_value_ori
                    minres=res_value_ori
                
                # New atom/res number
                atom_value=atom_value_ori-minatom+1
                res_value=res_value_ori-minres+1
                
                # Remaining values
                dif_atom=len(str(atom_value_ori))-len(str(atom_value))
                atom_str=" "*dif_atom+str(atom_value)
                dif_res=len(str(res_value_ori))-len(str(res_value))
                res_str=" "*dif_res+str(res_value)
                
                line=line[:atom_start+1]+atom_str+line[atom_end-1:res_start+1]+res_str+line[res_end-1:]
                n+=1
              
            # Write
            out.write(line+"\n")



