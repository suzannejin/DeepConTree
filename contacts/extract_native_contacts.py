
__description__ = "Extract native contacts from PDB file."

def __extract_native_contacts(pdb,out):
    
    import math
    
    tmp = {}
    o = open(out,"wt")

    with open(pdb) as p:
        for line in p:
            line = line.strip("\n")
            fields = line.split()
            fields = filter(None,fields)
            if len(fields) == 0:
                continue
            if ((fields[0] == "ATOM" and fields[2] == "CB" and fields[3] != "GLY") or 
                (fields[0] == "ATOM" and fields[2] == "CA" and fields[3] == "GLY")):
                res = int(fields[5])
                x,y,z = float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())
                tmp[res] = [x,y,z]
         
    for i in tmp.keys():
        for j in tmp.keys():
            if i >= j:
                continue
            coord1 = tmp[i]
            coord2 = tmp[j]
            dist = math.sqrt( (coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2 )
            
            if dist <= 8:
                o.write("{} {} 0 8 1.00\n".format(i,j))
    
    o.close()



if __name__ == '__main__':

    import argparse
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("pdb",type=str,help="PDB file.")
    app.add_argument("out",type=str,help="Output file with native contacts written in RR format.")
    args = app.parse_args()
    
    __extract_native_contacts(args.pdb,args.out)

