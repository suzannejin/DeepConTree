

__description__='Transforms a contact map to a contact list.'



def __map_to_rr(map_filename,rr_filename):
  
    import numpy as np
    
    matx=np.loadtxt(map_filename)
    rr_file=open(rr_filename,'wt')
    
    length=len(matx)
    for i in range(length):
        for j in range(length):
            if j <= i:
                continue
            score=matx[i][j]
            rr_file.write('{} {} {} {} {}\n'.format(i+1,j+1,0,8,score))  
            
    rr_file.close()
      

if __name__ == "__main__":
  
    import argparse
    parser=argparse.ArgumentParser(description=__description__)
    # I/O
    parser.add_argument('contact_map',type=str,help='Contact map file')
    parser.add_argument('contact_rr',type=str,help='Contact list in RR format')
    args=parser.parse_args()
    
    __map_to_rr(args.contact_map,args.contact_rr)

  
