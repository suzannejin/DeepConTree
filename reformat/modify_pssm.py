
__description__ = "It adjusts the format of the original .pssm file, so that it can be used by DeepConPred2."



def __modify_pssm(input_name,output_name):

    import re
  
    in_pssm = open(input_name,'rt')
    out_pssm = open(output_name,'wt')
  
    for line in in_pssm:
        line = line.strip("\n")
        fields = line.split(" ")
        fields = filter(None,fields)
    
        header = re.match(r'(\s+\w\s{3}\w\s{3}\w\s{3}\w)',line)
        if header:
            out_pssm.write(" "*9)
            n=0
            for element in fields:
                if n < 20:
                    out_pssm.write("  {}".format(element))
                else:
                    out_pssm.write("   {}".format(element))
                n+=1
            out_pssm.write("\n")
  
        row = re.match(r'(\s+\w+\s\w\s+)',line)
        if row:
            out_pssm.write("{:>5} {}  ".format(fields[0],fields[1]))
            for i in range(2,22):
                out_pssm.write("{:>3}".format(fields[i]))
            out_pssm.write(" ")
            for i in range(22,len(fields)-2):
                out_pssm.write("{:>4}".format(fields[i]))
            out_pssm.write("{:>6}{:>5}\n".format(fields[-2],fields[-1]))
    
        if not header and not row:
            out_pssm.write(line+"\n")
  
    in_pssm.close()
    out_pssm.close()
  
  
if __name__ == "__main__":
  
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    # I/O
    parser.add_argument('input_pssm',type=str,help='Original pssm file')
    parser.add_argument('output_pssm',type=str,help='Modified pssm file')
    args = parser.parse_args()
    
    __modify_pssm(args.input_pssm,args.output_pssm)
