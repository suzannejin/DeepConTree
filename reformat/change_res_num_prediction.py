#
# This script is employed to change the residue numbers of the predicted structure according the to residue numbers of the experimental structure
# Usage: python3 change_res_num.py <predicted_structure> <experimental_structure> 
#

import sys
import os

rd = open(sys.argv[2],"r")
for exp in rd:
    if exp[:4] == 'ATOM':
        exp_line = [exp[:6], exp[6:11], exp[12:16], exp[17:20], exp[21], exp[22:26], exp[30:38], exp[38:46], exp[46:54]]
        init_num = int(exp_line[5])-1
        break
    
rd.close()

f = open(sys.argv[1],"r")
for line in f:
    if line[:4] == 'ATOM' or line[:6] == "HETATM":
      splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
      splitted_line[5] = int(splitted_line[5])+init_num
# To format again the pdb file with the fields extracted
      output= "%-6s%5s %4s %3s %s%4s    %8s%8s%8s" % tuple(splitted_line)
      print(output)
  else:
      sys.stdout.write(line)
f.close()
