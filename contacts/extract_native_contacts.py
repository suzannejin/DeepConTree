
__description__ = "Extract native contacts from PDB file."

def __extract_native_contacts(pdb,out):
    
    import math
    
    tmp = {}
    o = open(out,"wt")

    with open(pdb) as p:
        for line in p:
            line = line.split()
            line = filter(None,line)
            if len(line) == 0:
                continue
            if ((line[0] == "ATOM" and line[2] == "CB" and line[3] != "GLY") or 
                (line[0] == "ATOM" and line[2] == "CA" and line[3] == "GLY")):
                res = int(line[5])
                x,y,z = float(line[6]),float(line[7]),float(line[8])
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

