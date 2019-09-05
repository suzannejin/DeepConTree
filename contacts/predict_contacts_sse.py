

__description__ = '''

    Protein contacts prediction.
    
    It predicts contacts using: DeepContact, DeepConPred and/or DeepCov.
    
    It can also be used to predict SSEs by Spider3, if -only_sse is specified.
    
    Remember to modify the configurations in deepconpred.yaml and set_envs.sh before running this script.

'''



# ========================
# ===== GENERATE MSA =====


def __hhblits_aln(fasta_file,aln_file,modified_aln_file,config,cpu):
    ''' Generate alignments for DeepConPred2 and DeepCov. 
    
    In the case of DeepContact, it already contains a step that 
    generates MSA, although using parameters that are different from 
    the ones we used for DeepConPred2 and DeepCov.
    
    Hence, the DeepContact pipeline, 
    that I build by the moment, 
    takes as input a fasta file instead of an alignment.
    
    Considering the time factor, for the future it would be great 
    to find the optimal parameters to generate MSA for all predictors.
    '''
    
    description_prog = "HHBlits: generate MSA"
    command="{} -i {} -d {} -oa3m {} -n {} -diff {} -cov {} -cpu {}".format(
                    config.get("alignment","command"),
                    fasta_file,
                    config.get("alignment","database"),
                    aln_file,
                    config.get("alignment","niterations"),
                    config.get("alignment","diff"),
                    config.get("alignment","cov"),
                    cpu)
    __run_command(description_prog,command,separator="#"*60)
    
    # Conversion from a3m to psicov alignment format (and remove identical rows)
    # This is the alignment format accepted by DeepCov
    os.system("egrep -v '^>' {} | sed 's/[a-z]//g' | sort -u > {}".format(aln_file,modified_aln_file))
    
  
  
# ============================
# ===== PREDICT CONTACTS =====


def __predict_contacts(method,bin_dir,prot_id,fasta_file,aln_file,modified_aln_file,config,cpu):
    ''' Prediction of protein residue contacts.
    
    It sets the working environment and runs DeepCov/DeepConPred2/DeepContact pipeline. 
    '''
  
    pipeline = "{}/contacts/{}.py".format(bin_dir,method)
    env_file = "{}/set_envs.sh".format(bin_dir)
    out = "{}/{}.con".format(method,prot_id)
    
    description_prog = "Run {} pipeline".format(method).upper()
    
    if method == "deepcontact":
        command = "source {} {} && python {} {} {} -cpu {}".format(env_file,method,pipeline,fasta_file,out,cpu)
    elif method == "deepconpred":
        command = "source {} {} && python {} {} {} {} -cpu {}".format(env_file,method,pipeline,fasta_file,aln_file,out,cpu)
    else:
        pipeline = "{}/deepcov.sh".format(config.get("deepcov","path"))  
        command = "source {} {} && bash {} -m {} -r {} -i {} -o {}".format(                 
                    env_file,
                    method,
                    pipeline,
                    config.get("deepcov","model"),
                    config.get("deepcov","receptive_field"),
                    modified_aln_file,
                    out)
    __run_command(description_prog,command,separator="#"*60)




# ================
# ===== MAIN =====


if __name__ == '__main__':
  
    import os
    import datetime
    import argparse
    import configparser
    from other.run_command import __run_command
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("fasta_file",type=str,help="Sequence fasta file.")
    app.add_argument("out_dir",type=str,help="Output directory.")
    app.add_argument("-run","--run",nargs='+',choices=['deepcov','deepconpred','deepcontact'],help="Contact predictors.")
    app.add_argument("-cpu",type=int,default=1,help="Number of cpus. Default=1")
    app.add_argument("-only_sse",action="store_true",help="To only predict SSEs.")
    args = app.parse_args()
  
    
    # Configuration file
    bin_dir = os.path.dirname(os.path.abspath(__file__)) + "/.."
    yaml = "{}/configuration.yaml".format(bin_dir)
    config = configparser.ConfigParser()
    config.read(yaml)
  
  
    # Manage input & output
    fasta_file = os.path.abspath(args.fasta_file)
    prot_id = os.path.basename(args.fasta_file).split(".")[0]
    out_dir = os.path.abspath(args.out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)  # = Working directory    'contacts' or 'sse'
          
  
  
    # Generate MSA 
    if args.only_sse is True or "deepcov" in args.run or "deepconpred" in args.run:
        if not os.path.exists("aln/"):
            os.makedirs("aln/")
        aln_file = "aln/{}.a3m".format(prot_id)
        modified_aln_file = "aln/{}.psicov".format(prot_id)
        __hhblits_aln(fasta_file,aln_file,modified_aln_file,config,args.cpu)
    else:
        aln_file = ""
        modified_aln_file = ""
    
    
    # Only generate SSEs
    if args.only_sse is True:
        os.system("python {}/contacts/deepconpred.py {} {} {}/{}.spd33 -cpu {} -only_sse".format(bin_dir,fasta_file,aln_file,out_dir,prot_id,args.cpu))
        exit(0)
  
  
    # Run DeepCov, DeepConPred and/or DeepContact
    for method in args.run:
        if not os.path.exists(method):
            os.makedirs(method)
        __predict_contacts(method,bin_dir,prot_id,fasta_file,aln_file,modified_aln_file,config,args.cpu)
 
    
