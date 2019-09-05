

__description__="Unpickle contact map"



def __pkl_to_map(pkl_input_filename,matrix_output_filename):
  
    import pickle
    
    matrix_out_file=open(matrix_output_filename,'wt')
    
    with open(pkl_input_filename,'rb') as f:
        matx = pickle.load(f)
    
    matx_size=len(matx)
    for i in range(matx_size):
        line=""
        for j in range(matx_size):
            space=" "
            if j == matx_size-1:
                space="\n"
            line=line+str(repr(matx[i][j]))+space
        matrix_out_file.write(line)
  
    matrix_out_file.close()


if __name__ == "__main__":
  
    import argparse
    parser=argparse.ArgumentParser(description=__description__)
    # I/O
    parser.add_argument('pickled_map',type=str,help='Pickled contact map file')
    parser.add_argument('out_map',type=str,help='Output file (unpickled contact map)')
    args=parser.parse_args()
    
    __pkl_to_map(args.pickled_map,args.out_map)

  
