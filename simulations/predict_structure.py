

__description__ = '''

     Compute protein 3D structure using:
     
     + Xplor-NIH software
    
     + contact restraints
    
     + SSE restraints

'''


    
def __prepare_files(fasta_file,seq_file,con_file,con_new,sse_file,tbl,tbl_sse,args,config):
    ''' Prepare files needed to run simulations '''
  
    # Copy files to keep them inside the output folder
    os.system("cp {} .".format(fasta_file))
    os.system("cp {} .".format(sse_file))
  
    # Convert seq from "one-letter" to "three-letter" code
    __convert_one_to_three(fasta_file,seq_file)

    # Prepare contact files according to given score threshold
    os.system("awk '$2-$1>=" + str(args.seqsep) + "&&$5>=" + str(args.score) + "{print $0}' " + con_file + " | sort -k5gr > " + con_new)

    # Prepare restraint files
    if args.xplor_betaoo == "yes":
        restraint = "o-o"
    else:
        restraint = ""
    __prepare_restraints_contacts(con_new,fasta_file,tbl,restraint)
    __prepare_restraints_dihedral(sse_file,tbl_sse)
  
  
def __run_xplor(ID,seq_file,tbl,tbl_sse,out_dir,bin_dir,args,config):
    ''' Run simulated annealing using Xplor-NIH '''
    
    description_prog = "Running Xplor-NIH"
    annealing_out = ID + ".md1.out"
    

    script = bin_dir + "/simulations/annealing.py"
    command = "{} -smp {} -py -o {} {} {} {} {} {} {} {} {}".format(
                    config.get("xplor","command"),args.cpu,
                    annealing_out,script,
                    seq_file,tbl,tbl_sse,
                    args.xplor_mode,out_dir,
                    args.xplor_nmodels,
                    args.xplor_topavg)
                            
    __run_command(description_prog,command," ")
    
    
  
  
# ================
# ===== MAIN =====


if __name__ == '__main__':
  
    import os
    import argparse
    import configparser
    import shutil
    from simulations.convert_one_to_three import __convert_one_to_three
    from simulations.prepare_restraints import __prepare_restraints_contacts
    from simulations.prepare_restraints import __prepare_restraints_dihedral
    from other.run_command import __run_command

    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("-fa","--fasta_file",type=str,required=True,help="Sequence file in fasta format.")
    app.add_argument("-con","--con_file",type=str,required=True,help="Contact file in RR format.")
    app.add_argument("-sse","--sse_file",type=str,required=True,help="Secondary structure elements in Spider3-local format.")
    app.add_argument("-o","--out_dir",type=str,required=True,help="Output directory.")
    app.add_argument("-cpu",type=int,default=1,help="Number of cpus. Default=1")
    app.add_argument("-score","--score",type=float,default=0,help="Score threshold (fifth field in the RR contact file). Default = 0")
    app.add_argument("-seqsep","--seqsep",type=int,default=8,
        help="Cutoff for the separation of residues i,j through the sequence. \
        Default = 8. Hence, contacts separated by less than 8 residues are removed.")  
    # Simulations - Xplor-NIH parameters
    app.add_argument("-xplor_mode",type=str,default="soft",choices=["soft","hard"])
    app.add_argument("-xplor_nmodels",type=int,default=250,help="Number of decoy models to be generated.")
    app.add_argument("-xplor_topavg",type=float,default=0.10,help="Rate of the top models that should be used to generate the average structure.")
    app.add_argument("-xplor_betaoo",type=str,default="no",choices=["no","yes"],help="Use beta-strand O-O restraints.")
    app.add_argument("-xplor_nexp",type=int,default=1,help="Number of experiments.")
    args = app.parse_args()
    
    # Files
    fasta_file = os.path.abspath(args.fasta_file)
    con_file = os.path.abspath(args.con_file)
    sse_file = os.path.abspath(args.sse_file)
    ID = os.path.basename(fasta_file).split(".")[0]
    seq_file = ID + ".seq"
    con_new = ID + ".con"
    tbl = con_new + ".tbl"
    tbl_sse = ID + "_sse.tbl"

    # Manage output directories
    out_dir = os.path.abspath(args.out_dir) 
    os.chdir(out_dir)  # This script will be working in the output directory

    # Configuration file
    bin_dir = os.path.dirname(os.path.abspath(__file__)) + "/.."
    yaml = "{}/configuration.yaml".format(bin_dir)
    config = configparser.ConfigParser()
    config.read(yaml)


    # Prepare files
    __prepare_files(fasta_file,seq_file,con_file,con_new,sse_file,tbl,tbl_sse,args,config)

    # Run Xplor-NIH
    __run_xplor(ID,seq_file,tbl,tbl_sse,out_dir,bin_dir,args,config)
    



