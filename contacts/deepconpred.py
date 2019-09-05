

__description__='''

    DeepConPred2 pipeline. Prediction of contacts.
    
    It can also be used to predict SSEs by Spider3, if -only_sse is specified.
    
    Remember to modify the configurations in deepconpred.yaml and set_envs.sh before running this script.

'''




# ================================
# ===== GENERATE INPUT FILES =====


def __hhblits(fasta_file,feature_out_prefix,config,cpu):
    ''' Use HHBlits to generate HMM profile.
    
    This information is used in Spider3 for the prediction of secondary structure.
    '''
  
    description_prog = "HHBlits: generation of HMM profile"
    command = "{} -i {} -ohhm {}.hhm -d {} -v{} -maxres {} -cpu {} -Z {}".format(
                config.get("hhblits-deepconpred","command"),
                fasta_file,
                feature_out_prefix,
                config.get("hhblits-deepconpred","database"),
                config.get("hhblits-deepconpred","v"),
                config.get("hhblits-deepconpred","maxres"),
                cpu,
                config.get("hhblits-deepconpred","z"))
  
    __run_command(description_prog,command,separator="="*60)
  

def __psiblast(fasta_file,feature_out_prefix,config,cpu):
    ''' Use PSI-Blast to create a position-specific scoring matrix (PSSM).
    
    This matrix is used in Spider3 for the prediction of secondary structure.
    It is also used in DeepConPred2 for the prediction of contacts. 
    It is modified before being used in DeepConPred2.
    '''
    
    description_prog = "PSI-Blast: generation of PSSM"
    command = "{} -db {} -num_iterations {} -num_alignments {} -num_threads {} -query {} -out {}.bla -out_ascii_pssm {}.pssm".format(
                config.get("psiblast-deepconpred","command"),
                config.get("psiblast-deepconpred","database"),
                config.get("psiblast-deepconpred","niterations"),
                config.get("psiblast-deepconpred","nalignments"),
                cpu,
                fasta_file,
                feature_out_prefix,
                feature_out_prefix)
    __run_command(description_prog,command,separator="="*60)
    
    # Transform PSSM format to be used in DeepConPred2
    original_pssm = "{}.pssm".format(feature_out_prefix)
    modified_pssm = "{}.PSSM".format(feature_out_prefix)
    __modify_pssm(original_pssm,modified_pssm)
  

def __spider3(fasta_file,prot_id,feature_out_prefix,config,env_file,outname,only_sse):
    ''' Use Spider3 to predict secondary structure.'''
    
    # Change to Spider3 directory
    spider_dir = config.get("spider3-deepconpred","path")
    os.chdir(spider_dir)
    
    # Copy input files to Spider3 directory
    a = [(fasta_file,prot_id),(feature_out_prefix+".pssm",prot_id+".pssm"),(feature_out_prefix+".hhm",prot_id+".hhm")]
    for files in a:
        shutil.copy(files[0],files[1])

    # Run Spider3
    os.system("echo {} {}.pssm {}.hhm > tlist".format(prot_id,prot_id,prot_id))
    description_prog = "Spider3: prediction of secondary structure. It is executed from Spider3 directory"
    command = "source {} spider && bash ./scripts/impute_script.sh tlist".format(env_file)
    __run_command(description_prog,command,separator="="*60)
    
    # Organize files
    shutil.move(prot_id+".spd33",feature_out_prefix+".spd33")
    for a in glob.glob(prot_id+"*"):
        os.remove(a)
    os.remove("tlist")
    if only_sse is True:
        shutil.copy(feature_out_prefix+".spd33",outname)

  
  
def __ccmpred(feature_out_prefix,aln_file,config,cpu):
    ''' Use CCMPred (plmDCA) to predict residue contacts.'''
    
    # Conversion from a3m to psicov alignment format (and remove identical rows)
    # This is the alignment format accepted by CCMPred
    modified_aln_file = "{}.psicov".format(feature_out_prefix)
    os.system("egrep -v '^>' {} | sed 's/[a-z]//g' | sort -u > {}".format(aln_file,modified_aln_file))
    
    # Predict contacts
    description_prog = "CCMPred: prediction of residue contacts based on plmDCA"
    command = "{} {} {}.ccmpred -t {}".format(
                config.get("ccmpred-deepconpred","command"),
                modified_aln_file,
                feature_out_prefix,
                cpu)
    __run_command(description_prog,command,separator="="*60)
   



# =======================
# ===== DEEPCONPRED =====


def __deepconpred(deepconpred_dir,fasta_file,prot_id,feature_out_prefix,prediction_map):
    ''' Refinement of contact map.'''
    
    os.chdir(deepconpred_dir)  # DeepConPred2 works in the directory where the program is installed
    data_dir = "{}/data".format(deepconpred_dir)
    result_dir = "{}/result".format(deepconpred_dir)
    
    # Copy input files to DeepConPred2 directory
    in_prefix = "{}/{}".format(data_dir,prot_id)
    shutil.copy(fasta_file,in_prefix)
    for suffix in [".PSSM",".spd33",".ccmpred"]:
        shutil.copy(feature_out_prefix+suffix,in_prefix+suffix)
    
    # Run DeepConPred2
    description_prog = "DeepConPred2: refinement of contact map. It is executed from DeepConPred2 directory"
    command = "python DeepConPred2.py {}".format(prot_id)
    __run_command(description_prog,command,separator="="*60)
    
    # Organize files
    shutil.move(result_dir+"/"+prot_id+".contactP",prediction_map)
    for a in glob.glob(data_dir+"/"+prot_id+"*"):
        os.remove(a)




# ================
# ===== MAIN =====


if __name__ == '__main__':
  
    import os
    import datetime
    import argparse
    import configparser
    import glob
    import shutil
    from other.run_command import __run_command
    from reformat.contact_map_to_list import __map_to_rr
    from reformat.modify_pssm import __modify_pssm
    
    # Parse arguments
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("fasta_file",type=str,help="Fasta file.")
    app.add_argument("alignment",type=str,help="Multiple sequence alignment file (in a3m format).")
    app.add_argument("out",type=str,help="Output file.")
    app.add_argument("-cpu",type=int,default=1,help="Number of cpus. Default=1")
    app.add_argument("-only_sse",action="store_true",help="Only predict SSE using Spider3.")
    args = app.parse_args()
    
    # Configuration file
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    yaml = "{}/../configuration.yaml".format(pipeline_dir)
    config = configparser.ConfigParser()
    config.read(yaml)
    deepconpred_dir = config.get("deepconpred","path")
    
    # Working environment
    env_file = "{}/../set_envs.sh".format(pipeline_dir)
    
    # Protein name 
    prot_id = os.path.basename(args.out).split(".")[0]
    
    # Feature & files
    fasta_file = os.path.abspath(args.fasta_file)
    aln_file = os.path.abspath(args.alignment)
    out_dir = os.path.dirname(os.path.abspath(args.out))
    out = os.path.abspath(args.out)   # contacts/deepconpred/PROT.con | sse/PROT.sse
    feature_dir = "{}/{}_deepconpred_features".format(out_dir,prot_id)  # contacts/deepconpred/PROT_deepconpred_features
    if args.only_sse is True:
        feature_dir = "{}/{}_spider_features".format(out_dir,prot_id)  # sse/PROT_spider_features
    feature_out_prefix = "{}/{}".format(feature_dir,prot_id)
    prediction_map = "{}.mat".format(feature_out_prefix)
    prediction_con = out
    
    
    # Create folders
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(feature_dir):
        os.makedirs(feature_dir)


    
    # HHBlits: generation of HMM profile
    __hhblits(fasta_file,feature_out_prefix,config,args.cpu)
    
    # PSI-Blast: generation of PSSM 
    __psiblast(fasta_file,feature_out_prefix,config,args.cpu)
  
    # Spider3: prediction of secondary structure
    __spider3(fasta_file,prot_id,feature_out_prefix,config,env_file,out,args.only_sse)
  
  
    if args.only_sse is True: # Finish if user only wants to predict SSEs
        exit(0)
  
    # CCMPred: Prediction of residue contacts based on plmDCA
    __ccmpred(feature_out_prefix,aln_file,config,args.cpu)
  
    # DeepConPred: Refinement of contact map
    __deepconpred(deepconpred_dir,fasta_file,prot_id,feature_out_prefix,prediction_map)
    
    # Convert contact map to list (RR format)
    __map_to_rr(prediction_map,prediction_con)
