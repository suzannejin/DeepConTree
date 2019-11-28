

__description__='''

    Compute 1d & 3d trees.

'''


def __alignment(aln_script,fasta_file,ref_file,pfam,thread):
    ''' Compute alignment using TMalign '''
    os.system("bash {} {} {} {} {}".format(aln_script,fasta_file,ref_file,pfam,thread))
  
  
def __1d_tree(pfam,config):
    ''' Create alignment based tree '''  ## Question: if the alignment used is built considering structural information, then is it still 1d-tree?
    
    command = "{} -i out/aln/{}_tmalign.ph -o out/tree1d/{}.nwk -m {} {} out/tree1d/{}.boostrap".format(
                    config.get("tree1d","command"),
                    pfam,pfam,
                    config.get("tree1d","method"),
                    config.get("tree1d","parameters"),
                    pfam)
    os.system(command)
  
  
def __3d_tree(tree3d_script,pfam,ref_file):
    ''' Create Phylo-3d tree 
    
    
    msa:   MSA method
    mode:  export T-Coffee MODE
    mexp:  export T-Coffee MODE $EXP (2=square; 3=cubic)
    '''
    
    msa,mode,mexp = config.get("tree3d","method"),config.get("tree3d","mode"),config.get("tree3d","mexp")
    os.system("bash {} {} {} {} {} {}".format(tree3d_script,pfam,ref_file,msa,mode,mexp))
    
  


# ================
# ===== MAIN =====


if __name__ == '__main__':
  
    import os
    import datetime
    import argparse
    import configparser
    from other.run_command import __run_command
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("fa",type=str,help="Multi-fasta file.")
    app.add_argument("ref_file",type=str,help="Reference/template file. Each line from the reference file should contain the UniProt identifier, a separator and the PDB identifier. For example: >AMCY_PARDE/43-131 _P_ 1MDAA-1 .")
    app.add_argument("pdb_dir",type=str,help="Directory where all the pdb files are stored. Output files will be created here as well.")
    app.add_argument("thread",type=int,help="Number of threads running T-coffee for the alignment of sequences.")
    args = app.parse_args()
    
    # Scripts
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    aln_script = pipeline_dir + "/compute_alignment.sh"
    tree3d_script = pipeline_dir + "/compute_3dtree.sh"
    
    # Input files & folders
    fasta_file = os.path.abspath(args.fa)
    ref_file = os.path.abspath(args.ref_file)
    pdb_dir = os.path.abspath(args.pdb_dir)
    pfam = os.path.basename(fasta_file).split(".")[0].replace("_uniprot","").replace("_pdb","")  
    
    
    # Configuration file
    yaml = pipeline_dir + "/../configuration.yaml"
    config = configparser.ConfigParser()
    config.read(yaml)
    
    
    os.chdir(pdb_dir)  # This script works in the input directory where pdb files are stored
    
    # Create output directories
    for i in ["out/aln","out/tree1d","out/tree3d/"+pfam]:
        if not os.path.exists(i):
            os.makedirs(i)
    
    # Compute MSA
    __alignment(aln_script,fasta_file,ref_file,pfam,args.thread)
    
    # Compute trees
    __1d_tree(pfam,config)
    __3d_tree(tree3d_script,pfam,ref_file)
  
  
