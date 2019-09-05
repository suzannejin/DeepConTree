

__description__='''

    DeepContact pipeline.
    
    Remember to modify the configurations in deepconpred.yaml and set_envs.sh before running this script.

'''




# =======================
# ===== DEEPCONTACT =====


def __deepcontact():
    ''' DeepContact requires three steps: generation of input files, extract features from files, prediction of contacts.'''
    
    # Generation of input files
    description_prog = "Generation of input files"
    command="python {} {} {} {}".format(run_pipeline,default_yaml,fasta_file,feature_dir)
    __run_command(description_prog,command,separator=" ")
    
    # Extract features from input files
    description_prog = "Extract features from input files"
    command="python {} {} {} {} {}".format(feature_generation,feature_yaml,feature_dir,prot_id,feature_pkl)  # The original script does not need ProtID.
    __run_command(description_prog,command,separator=" ")                                                          # The script is changed to manipulate input filenames.
    
    # Prediction of residue contacts
    description_prog = "Prediction of protein residue contacts"
    command="python {} {} {}".format(main_script,feature_pkl,prediction_pkl)
    __run_command(description_prog,command,separator=" ")



# ================
# ===== MAIN =====


if __name__ == '__main__':
    
    import os
    import datetime
    import argparse
    import configparser  
    import shutil
    from write_deepcontact_parameters import __write_deepcontact_yaml
    from other.run_command import __run_command
    from reformat.pkl_to_matrix import __pkl_to_map
    from reformat.contact_map_to_list import __map_to_rr
    
    # Parse arguments
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("fasta_file",type=str,help="Fasta file.")
    app.add_argument("out",type=str,help="Output contact file (in RR format).")
    app.add_argument("-cpu",type=int,default=1,help="Number of cpus. Default=1")
    args = app.parse_args()
    
    # Configuration file
    bin_dir = os.path.dirname(os.path.abspath(__file__)) + "/.."
    yaml = "{}/configuration.yaml".format(bin_dir)
    config = configparser.ConfigParser()
    config.read(yaml)
    deepcontact_dir = config.get("deepcontact","path")
    
    # Protein name 
    prot_id = os.path.basename(args.fasta_file).split(".")[0]
    
    # Features & files
    prediction_con = os.path.abspath(args.out)
    out_dir = os.path.dirname(prediction_con)
    feature_dir = "{}/{}_deepcontact_features".format(out_dir,prot_id)
    feature_pkl = "{}/{}_features.pkl".format(feature_dir,prot_id)
    prediction_pkl = "{}/{}_predictions.pkl".format(feature_dir,prot_id)
    prediction_map = "{}/{}_predictions.mat".format(feature_dir,prot_id)
    fasta_file = "{}/{}.fasta".format(feature_dir,prot_id)  # Extension '.fasta'. Otherwise, DeepContact will give errors.

    # DeepContact scripts
    run_pipeline = "{}/data-processing/run_pipeline.py".format(deepcontact_dir)
    feature_generation = "{}/deepcontact/feature_gen.py".format(deepcontact_dir)
    main_script = "{}/deepcontact/main.py".format(deepcontact_dir)

    # DeepContact configuration yaml files
    default_yaml = "{}/default.yaml".format(deepcontact_dir)
    feature_yaml = "{}/deepcontact/feature.yaml".format(deepcontact_dir)

  
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(feature_dir):
        os.makedirs(feature_dir)
    shutil.copy(args.fasta_file,fasta_file)
    
    
    # Write default.yaml file
    __write_deepcontact_yaml(default_yaml,config,args.cpu)
    
    
    # Run DeepContact
    __deepcontact()
    
    # Reformat contact files
    __pkl_to_map(prediction_pkl,prediction_map)
    __map_to_rr(prediction_map,prediction_con)

