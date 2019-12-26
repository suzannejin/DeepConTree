

__description__ = '''

    Pipeline: given a set of sequences with known 3D structures, 
    
    protein contacts & structures,
    
    and structure based trees 
    
    are computed.
    
    After all, the RF metric between the generated 3D-trees are calculated.
    
    
    You can also just compute contacts, or structures, or trees separately.
    
    
    Note: remember to modify the configurations in deepconpred.yaml and set_envs.sh before running this script.

'''



# =================================
# ==== Parse Arguments & Paths ==== 

def __parse_arguments():
    
    import sys
    import argparse
    import re
    
    
    app = argparse.ArgumentParser(description=__description__)
    
    # Input & output
    app.add_argument("-seq","--seq",type=str,required=True,help="Multi-fasta file.")
    app.add_argument("-pdbf","--pdb_dir",type=str,required=True,help="Directory where pdb files (suffix .pdb) are stored.")
    app.add_argument("-ref","--ref_file","--template_file",type=str,required=True,help="Template/reference file.")
    app.add_argument("-o","--out_dir",type=str,required=True,help="Output directory.") 
    
    # CPU
    app.add_argument("-cpu",type=int,default=1,help="Number of cpus. Default=1")
    
    # Actions
    app.add_argument("-actions",nargs='+',choices=['contacts','simulations','trees'],required=True)
    
    # Contact predictors
    app.add_argument("-method",nargs='+',choices=['deepcov','deepconpred','deepcontact'],required='contacts' in sys.argv,help="Contact predictors.")
    
    # Filtering of contacts
    app.add_argument("-pdeepcov",type=float,default=0.50,help="Probability score threshold for DeepCov. Default=0.50")
    app.add_argument("-pdeepconpred",type=float,default=0.50,help="Probability score threshold for DeepConPred2. Default=0.50")
    app.add_argument("-pdeepcontact",type=float,default=0.97,help="Probability score threshold for DeepContact. Default=0.97")
    app.add_argument("-seqsep",type=int,default=8,help="Minimum separation of a pair of residues in the sequence. Default=8")
    # Cscore
    app.add_argument("-cscore",type=float,default=None,help="Conservation score threshold. Default=None")
    app.add_argument("-c_aln",type=str,required='-cscore' in sys.argv,help="Alignment for the estimation of conservation score in fasta format (suffix .fa)")
    app.add_argument("-cweight",type=str,default="no",choices=["yes","no"],help="Apply Voronoid weight for the estimation of conservation score. Default=no")
    
    # Provide contacts and SSE files for simulations
    app.add_argument("-con_list",type=str,required='simulations' in sys.argv and 'contacts' not in sys.argv,help="Contact list with three columns separated by space: /absolute/path/to/contact_folder folder_name probability_threshold.")
    app.add_argument("-ssef",type=str,help="Path to SSE predictions (same format as Spider3 output). Each file should be named according to the PDB files, and with extension '.sse'. Eg.1MDAA-1.sse")
             
    # Provide structures for tree reconstruction
    app.add_argument("-str_list",type=str,required='trees' in sys.argv and 'simulations' not in sys.argv,help="Structure list with two columns: /absolute/path/to/structure_folder folder_name.")

    # Simulations - Xplor-NIH parameters
    app.add_argument("-xplor_mode",type=str,default="soft",choices=["soft","hard"],help="Xplor simulations mode. Default=soft")
    app.add_argument("-xplor_nmodels",type=int,default=250,help="Number of decoy models to be generated. Default=250")
    app.add_argument("-xplor_topavg",type=float,default=0.10,help="Rate of the top models that should be used to generate the average structure. Default=0.10")
    app.add_argument("-xplor_betaoo",type=str,default="no",choices=["no","yes"],help="Use beta-strand O-O restraints. Default=no")
    app.add_argument("-xplor_nexp",type=int,default=1,help="Number of experiments. Default=1")

    args = app.parse_args()
    
    
    return args



def __organize_paths(args):
    
    paths = {}
    scripts = {}
    
    
    paths["bin_dir"] = os.path.dirname(os.path.abspath(__file__)) 
    
    # Configuration file
    paths["yaml"] = "{}/configuration.yaml".format(paths["bin_dir"])  

  
    # Output directory
    paths["out_dir"] = os.path.abspath(args.out_dir)
    if not os.path.exists(paths["out_dir"]):
        os.makedirs(paths["out_dir"])
  
    # Input directories & files
    paths["ref_file"] = os.path.abspath(args.ref_file)
    paths["multifasta"] = os.path.abspath(args.seq)
    paths["fasta_name"] = os.path.basename(paths["multifasta"]).split(".")[0]
    paths["fastapdb"] = paths["out_dir"] + "/" + paths["fasta_name"] + "_pdb.fa"
    paths["fastauniprot"] = paths["out_dir"] + "/" + paths["fasta_name"] + "_uniprot.fa"
    paths["pdb_dir"] = os.path.abspath(args.pdb_dir)
    a = [args.con_list,args.str_list,args.ssef,args.c_aln]
    b = ["con_list","str_list","sse_dir","c_aln"]
    for i in range(len(a)):
        if a[i] is not None:
            paths[b[i]] = os.path.abspath(a[i])
        else:
            paths[b[i]] = None


    # Scripts 
    scripts["environment"] = paths["bin_dir"] + "/set_envs.sh"
    scripts["contacts_sse"] = paths["bin_dir"] + "/contacts/predict_contacts_sse.py"
    scripts["xplor"] = paths["bin_dir"] + "/simulations/predict_structure.py" 
    scripts["trees"] = paths["bin_dir"] + "/trees/trees.py"
    scripts["rf"] = paths["bin_dir"] + "/trees/rf.sh"
    scripts["reformat_pdb"] = paths["bin_dir"] + "/reformat/reformat_pdb.sh"
    scripts["cscore"] = paths["bin_dir"] + "/contacts/cscore.sh"
  
  
    return(paths,scripts)
            
        


# ==================
# ==== CONTACTS ====        


def __compute_contacts(contact_script,methods,cpu,ids):
    ''' Extract native contacts and predict contacts using different methodologies specified by user '''
  
  
    # Extract native contacts 
    __extract_native(ids)
  

    # Predict contacts using DeepCov, DeepConPred2 and DeepContact for each protein
    run = ""
    for method in methods:
        run += method + " "
    for prot in ids:
        description_prog = "Predicting contacts for " + prot + " using " + str(methods) + "\n"
        command = "python {} seq/{}.fa contacts/ -run {} -cpu {}".format(contact_script,prot,run,cpu)
        __run_command(description_prog,command,separator="="*60)
    
    
    # Contactbench
    for method in methods:
        for prot in ids:
            prediction = "contacts/" + method + "/" + prot + ".con"
            native = "contacts/native/" + prot + ".con"
            out = "contacts/" + method + "/" + prot + ".contactbench"
            __rr_to_contactbench(prediction,native,out)



def __extract_native(ids):
    ''' Extract native contacts '''
    
    if not os.path.exists("contacts/native"):
        os.makedirs("contacts/native")
    
    for prot in ids:
        pdb = "pdb_matched/" + prot + ".pdb"
        out = "contacts/native/" + prot + ".con"
        __extract_native_contacts(pdb,out)




# ====================
# ==== STRUCTURES ====

def __compute_structures(config,paths,scripts,args,ids,cids):
    ''' Compute 3D structure through simulated annealing 
    
    and constrained by contact restraints and dihedral angle restraints.
    '''
  
    # List of contacts
    if "contacts" in args.actions: 
        paths["con_list"] = paths["out_dir"] + "/contacts/con_list"
        f = open(paths["con_list"],"wt")
        #f.write(paths["out_dir"] + "/contacts/native native 0\n")   # Reconstruct structures from native contacts
        for method in args.method:
            if method == "deepcov":
                score = args.pdeepcov
            elif method == "deepconpred":
                score = args.pdeepconpred
            elif method == "deepcontact":
                score = args.pdeepcontact
            if args.cscore is None:
                line = paths["out_dir"] + "/contacts/" + method + " " + method + "_score" + str(score).replace(".","-") + "_seqsep" + str(args.seqsep) + " " + str(score) + "\n"
            else:
                line = paths["out_dir"] + "/contacts/" + method + " " + method + "_score" + str(score).replace(".","-") + "_seqsep" + str(args.seqsep) + "_cscore" + str(args.cscore).replace(".","-") + " " + str(score) + "\n"
            f.write(line)
        f.close()
  
    
    # Filter by conservation score
    if args.cscore is not None:
        __apply_cscore(args,paths,scripts,ids,cids)
  
  
    # Manage SSE files
    __manage_sse(paths["sse_dir"],scripts["contacts_sse"],ids,args.actions,args.method,args.cpu)


    # List of structures
    if "trees" in args.actions:
        paths["str_list"] = paths["out_dir"] + "/simulations/str_list"
        str_list = open(paths["str_list"],"wt")

    # Run Xplor-NIH n times for each protein & method
    nexperiment = int(args.xplor_nexp)
    for n in range(1,nexperiment+1):
        # Read contacts folder, name, threshold {name : (path/to/folder, threshold)}
        cons = {}
        with open(paths["con_list"]) as f:
            for line in f:
                line = filter(None,line.strip().split(" "))
                if len(line) == 0:
                    continue
                folder = line[0]
                name = line[1]
                threshold = line[2]
                cons[name] = (folder,threshold)
        # Run simulations
        for name in cons.keys():
            confolder = cons[name][0]
            threshold = cons[name][1]
            for prot in ids:
                fasta_file = "seq/" + prot + ".fa"
                con_file = confolder + "/" + prot + ".con"
                sse_file = "sse/" + prot + ".sse"
                out_dir = "simulations/" + name + "/n" + str(n) + "/" + prot + "_models" 
                if not os.path.isfile(con_file):
                    print("Warning: There is no " + con_file)
                    continue
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                description_prog = "Run Xplor-NIH for " + prot + " (" + name + ") score>=" + threshold + " seqsep>=" + str(args.seqsep)
                command = "python {} -fa {} -con {} -sse {} -o {} -cpu {} -score {} -seqsep {} -xplor_mode {} -xplor_nmodels {} -xplor_topavg {} -xplor_betaoo {} -xplor_nexp {}".format(
                                    scripts["xplor"],fasta_file,
                                    con_file,sse_file,out_dir,args.cpu,
                                    threshold,args.seqsep,
                                    args.xplor_mode,args.xplor_nmodels,args.xplor_topavg,args.xplor_betaoo,args.xplor_nexp)
                __run_command(description_prog,command,"="*60)
                a = out_dir + "/annealing_ave.pdb"
                b = "simulations/" + name + "/n" + str(n) + "/" + prot + ".pdb"
                shutil.copy(a,b)
                print(command)
    
            # Write list of structures
            if "trees" in args.actions:
                str_list.write(paths["out_dir"] + "/simulations/" + name + "/n" + str(n) + " " + name + "_n" + str(n) + "\n")
    if "trees" in args.actions:
        str_list.close()



def __manage_sse(sse_dir,script,ids,actions,methods,cpu):    
    ''' Manage SSE files.'''
    
    # Create sse folder
    if not os.path.exists("sse"):
        os.makedirs("sse")
        
    # If SSE files are provided by user    
    if sse_dir is not None:
        for prot in ids:
            a = sse_dir+"/"+prot+".sse"
            b = "sse/"+prot+".sse"
            if os.path.abspath(a) != os.path.abspath(b):
                shutil.copy(a,b)
    else:
        # Copy SSE files predicted by Spider3 
        if "contacts" in actions and "deepconpred" in methods:
            for prot in ids:
                a = "contacts/deepconpred/"+prot+"_deepconpred_features/"+prot+".spd33"
                b = "sse/"+prot+".sse"
                shutil.copy(a,b)
        # Predict SSE
        else:
            for prot in ids:
                description_prog = "Predicting SSEs with Spider3 for " + prot
                command = "python {} seq/{}.fa sse/ -cpu {} -only_sse".format(script,prot,cpu)
                __run_command(description_prog,command,separator="="*60)
                os.rename("sse/"+prot+".spd33","sse/"+prot+".sse")



def __get_cids(args,paths,scripts):
    
    # Manage alignment file
    lines = open(paths["c_aln"]).readlines()
    line = lines[0].strip().split(" ")
    calign = paths["out_dir"] + "/" + os.path.basename(paths["c_aln"]).split(".")[0] + "_calign.fa"
    if not os.path.isfile(paths["pdb_dir"] + "/" + line[0][1:] + ".pdb"):
        __map_names(paths["c_aln"],paths["ref_file"],calign,"uniprot_to_pdb")
    else:
        shutil.copy(paths["c_aln"],calign)
    paths["c_aln"] = calign
        
    
    # Get proteins in the alignment file
    prots = []
    with open(paths["c_aln"]) as f:
        for line in f:
            line = line.strip().split(" ")[0]
            if line[0] == ">":
                prot = line[1:]
                if os.path.isfile(paths["pdb_dir"] + "/" + prot + ".pdb") is True:
                    prots.append(line[1:])
    return(prots)
    
    

def __apply_cscore(args,paths,scripts,ids,cids):
    ''' Calculate the conservation of the predicted contacts and filter them according to the Cscore'''
                

    # Copy contacts and create contactbench files
    if "contacts" not in args.actions:
        __extract_native(cids)
        with open(paths["con_list"]) as f:
            for line in f:
                line = line.strip().split(" ")
                folder = line[0]
                name = line[1]
                if name == "native":
                    continue
                if not os.path.exists("contacts/"+name):
                    os.makedirs("contacts/"+name)
                for prot in cids:
                    pred = folder + "/" + prot + ".con"
                    if os.path.isfile(pred):
                        prediction = "contacts/"+name+"/"+prot+".con"
                        native = "contacts/native/"+prot+".con"
                        contactbench = "contacts/"+name+"/"+prot+".contactbench"
                        shutil.copy(pred,prediction)
                        __rr_to_contactbench(prediction,native,contactbench)

    # Read con_list {name : (path/to/folder, threshold)}
    con_list = {}
    with open(paths["con_list"]) as f:
        for line in f:
            line = line.strip().split(" ")
            name = line[1]
            p = float(line[2])
            if "contacts" in args.actions:
                folder = line[0] 
            else:
                folder = paths["out_dir"] + "/contacts/" + name
            # Create folder
            if not os.path.exists(folder + "/cscore") and name != "native":
                os.makedirs(folder + "/cscore")
            con_list[name] = (folder,p)

    # Calculation of conservation score
    paths["con_list"] = paths["out_dir"] + "/contacts/con_list"
    w = open(paths["con_list"],"wt")
    for name in con_list.keys():
        if name == "native":
            w.write(paths["out_dir"] + "/contacts/native native 0\n")
            continue
        folder = con_list[name][0]
        p = con_list[name][1]
        # Filter contactbench files
        for prot in cids:
            outprefix = folder + "/cscore/" + prot
            contactbench1 = folder + "/" + prot + ".contactbench"
            contactbench2 = outprefix + ".contactbench"
            if not os.path.isfile(contactbench1):
                print("Warning: no contact file for " + prot + " " + name +" is provided")
                print("Warning: The conservation score will be computed without considering the contact predictions of the protein " + prot)
            else:
                cont = open(contactbench2,"wt")
                with open(contactbench1) as f:
                    for line in f:
                        fields = line.strip().split(" ")
                        if int(fields[2]) >= args.seqsep and float(fields[3]) >= p:
                            cont.write(line)
                cont.close()
        # Estimate Cscore
        for prot in ids:
            outprefix = folder + "/cscore/" + prot
            contactbench = outprefix + ".contactbench"
            description_prog = "Estimating the conservation of the predicted contacts"
            command = "bash {} {} {} {} {} {} {}".format(scripts["cscore"],paths["c_aln"],contactbench,outprefix,p,args.cscore,args.cweight)
            __run_command(description_prog,command,separator=" ")
        # Rewrite con_list
        folder = folder + "/cscore"
        w.write(folder + " " + name + " " + str(p) + "\n")
    w.close()
    

        
    
# ==================
# ==== 3D TREES ====

def __compute_trees(config,paths,scripts,cpu,ids):
    ''' Compute structure-based trees '''
    
    
    # Read list of structures into a dic {name : path/to/folder}
    strdic = {}
    with open(paths["str_list"]) as f:
        for line in f:
            line = line.strip().split(" ")
            if line[1] != "experimental":
                strdic[line[1]] = line[0]
    
    # Experimental structures
    strdic["experimental"] = paths["out_dir"] + "/pdb_matched"
    
    # List of trees with two columns: 1) path to tree.nwk  and  2) tree name
    if not os.path.exists("trees"):
        os.makedirs("trees")
    paths["tree_list"] = paths["out_dir"] + "/trees/tree_list"
    f = open(paths["tree_list"],"wt")
    msa,mode,mexp = "tmalign",4,2
    treename = paths["fasta_name"] + "_unw_d4_" + msa + str(mode) + "-" + str(mexp) + ".mat_1.txt.nwk"
    f.write(paths["out_dir"] + "/trees/experimental/out/tree3d/" + treename + " experimental\n")
    
    # Compute 3D-tree for every set of structures
    for name in strdic.keys():
        orifolder = strdic[name]
        out_dir = paths["out_dir"] + "/trees/" + name
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for prot in ids:
            if name == "experimental":
                shutil.copy(orifolder+"/"+prot+".pdb",out_dir+"/"+prot+".pdb")
            else:
                # Reformat pdb files
                prediction = orifolder+"/"+prot+".pdb"
                experimental = "pdb_matched/" + prot + ".pdb"
                out = out_dir+"/"+prot+".pdb"
                os.system("bash {} {} {} {}".format(scripts["reformat_pdb"],prediction,experimental,out))
        description_prog = "Compute trees for " + name
        command = "source {} tree && python {} {} {} {} {}".format(scripts["environment"],scripts["trees"],paths["fastauniprot"],paths["ref_file"],out_dir,cpu)
        __run_command(description_prog,command,"="*60)
        
        # Write list of tree with two columns: 1) path to tree.nwk  and  2) tree name
        if name != "experimental":
            f.write(paths["out_dir"] + "/trees/" + name + "/out/tree3d/" + treename + " " + name + "\n")
    f.close()
    
    # Calculate RF score
    description_prog = "Calculate RF score"
    command = "bash {} {} {} {} {}".format(
                    scripts["rf"],
                    paths["fasta_name"],
                    config.get("tree3d","method"),
                    config.get("tree3d","mode"),
                    config.get("tree3d","mexp"))
    __run_command(description_prog,command,"="*60)
    
                            



# ================
# ===== MAIN =====


if __name__ == '__main__':
  
    import os,sys
    import glob
    import configparser
    import shutil
    from contacts.extract_native_contacts import __extract_native_contacts
    from reformat.match_positions_pdb import __match_pdb
    from reformat.mapping_names import __map_names
    from reformat.rr_to_contactbench import __rr_to_contactbench
    from other.run_command import __run_command


    # Parse arguments
    args = __parse_arguments()
    
    # Organize paths
    paths,scripts = __organize_paths(args)

    # Configuration file
    config = configparser.ConfigParser()
    config.read(paths["yaml"])

    # Change working directory to output directory
    os.chdir(paths["out_dir"]) 
  
    # Create subdirectories for easy use
    for i in ["seq","pdb","pdb_matched"]:
        if not os.path.exists(i):
            os.makedirs(i)
    for i in args.actions:
        if not os.path.exists(i):
            os.makedirs(i)


    # Manage files for easy use 
    shutil.copy(paths["ref_file"],".")
    lines = open(paths["multifasta"]).readlines()
    name1 = lines[0].strip().split(" ")[0][1:]  # Name of the first sequence
    refnames = {}  # Ref file in a dictionary
    with open(paths["ref_file"]) as ref:
        for line in ref:
            fields = line.strip("\n").split(" ")
            refnames[fields[0][1:]] = fields[1]
    #If multifasta named with uniprot id
    if name1 in refnames.keys():
        # Create multi-fasta file named with PDB ID
        __map_names(paths["multifasta"],paths["ref_file"],paths["fastapdb"],"uniprot_to_pdb")
        shutil.copy(paths["multifasta"],paths["fastauniprot"])
    # If multifasta named with pdb id
    else:
        __map_names(paths["multifasta"],paths["ref_file"],paths["fastauniprot"],"pdb_to_uniprot")
        shutil.copy(paths["multifasta"],paths["fastapdb"])
        
    # Extract sequence files
    ids = []  # List of prot pdb ids
    with open(paths["fastapdb"]) as f:
        for line in f:
            if line == "\n":
                continue
            line = line.strip()
            if line[0] == ">":
                name = line[1:]
                ids.append(name)
                o = open("seq/{}.fa".format(name),"wt")
                o.write(line + "\n")
            else:
                o.write(line + "\n")
                o.close()
                
    # Manage PDB files
    if args.cscore is not None:
        cids = __get_cids(args,paths,scripts)
        a = cids
    else:
        cids = None
        a = ids
    for prot in a:
        # Copy original files
        shutil.copy(paths["pdb_dir"]+"/"+prot+".pdb","pdb/.")    
        # Match atom/residue positions starting with 1
        __match_pdb("pdb/"+prot+".pdb","pdb_matched/"+prot+".pdb")  
    

    
    # Predict contacts
    if "contacts" in args.actions:   
        __compute_contacts(scripts["contacts_sse"],args.method,args.cpu,ids)
    
    
    # Predict 3D structures
    if "simulations" in args.actions:
        __compute_structures(config,paths,scripts,args,ids,cids)
        
        
    # Reconstruct structural trees
    if "trees" in args.actions:
        __compute_trees(config,paths,scripts,args.cpu,ids)
        
