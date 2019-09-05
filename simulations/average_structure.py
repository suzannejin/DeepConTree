
__description__="Compute average pdb file"


def __avg_pdb(pdb_list,out_name):
  
  pdbs = []
  npdbs = len(pdb_list)
  
  for n in range(npdbs):
    pdb = pdb_list[n]
    pdbs.append([])
    with open(pdb) as f:
      for line in f:
        if line[0:4] == "ATOM":
          line = filter(None,line.strip().split(" "))
          pdbs[n].append(line)
          
  avg = []
  nlines = len(pdbs[0])
  for nline in range(nlines):
    sum_x,sum_y,sum_z = 0,0,0
    for npdb in range(npdbs):
      line = pdbs[npdb][nline]
      sum_x += line[5]
      sum_y += line[6]
      sum_z += line[7]
    x = sum_x/npdbs
    y = sum_y/npdbs
    z = sum_z/npdbs
    


if __name__ == '__main__':
  
  import os
  import argparse
  
  app = argparse.ArgumentParser(description=__description__)
  app.add_argument("-pdb","--pdb",nargs="+",required=True,help="One or more pdb files.")
  app.add_argument("-o","--out",required=True,help="Output average pdb file.")
  args = app.parse_args()

  __avg_pdb(args.pdb,args.out)
