
__description__ = "Create a bash script to compute structure based tree using t-coffee : Phylo-3D"


def __print_phylo3d_uw_d4(wd,pfam,msa,mode,mode_exp):
  
  import os
  
  aln = wd + "/in/" + pfam + ".fa"
  ref = wd + "/in/" + pfam + ".template_list"
  pdb_dir = wd + "/pdb"
  out_dir = wd + "/out/tree3d"
  outbash = "{}/{}_{}_phylo3D_unw_d4-{}-{}.sh".format(out_dir,pfam,msa,mode,mode_exp)
  extension = "{}_unw_d4_{}_{}-{}".format(pfam,msa,mode,mode_exp)
  
  with open(outbash,"wt") as out:
    out.write("export THREED_TREE_MODE="+mode+"\n")
    out.write("export THREED_TREE_MODE_EXP="+mode_exp+"\n")
    out.write("export THREED_TREE_NO_WEIGHTS=1\n")
    out.write("cd "+pdb_dir+"\n")
    #out.write("t_coffee -other_pg seq_reformat -in {} -in2 {} \
#-action +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates \
#-output newick > {}/{}.trees\n".format(aln,ref,out_dir,extension))
    out.write("t_coffee -other_pg seq_reformat -in {} -in2 {} \
-action +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates \
-output dm > {}/{}.matrices\n".format(aln,ref,out_dir,extension))

  


if __name__ == '__main__':
  
  import argparse
  
  app = argparse.ArgumentParser(description=__description__)
  app.add_argument("-dir","--main_dir",type=str,help="Main directory with the corresponding folders: in, pdb, out. Fasta alignment and template list file are contained in 'in' folder, pdb files are in 'pdb' folder, all the files generated from the trmsd pipeline are store in 'out' folder.")
  app.add_argument("-pfam","--pfam",type=str,help="Pfam family id.")
  app.add_argument("-msa","--msa_method",type=str,help="MSA methods to apply.")
  app.add_argument("-mode","--mode",type=str,help="Export T-Coffee MODE.")
  app.add_argument("-mexp","--mode_exp",type=str,help="Export T-Coffee MODE EXP (2=square; 3=cubic)")
  args = app.parse_args()
  
  __print_phylo3d_uw_d4(args.main_dir,args.pfam,args.msa_method,args.mode,args.mode_exp)
  
