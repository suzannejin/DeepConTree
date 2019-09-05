

__description__ = "Write the 'default.yaml' configuration file into DeepContact directory."


def __write_deepcontact_yaml(default_yaml,config,cpu):

    out = open(default_yaml,"wt")

    # First lines
    out.write( "gen_feature_with_jackhmmer: " + config.get("aln-deepcontact","gen_feature_with_jackhmmer") + "\n\n")
    out.write( "path:\n        output_base: output\n\n" )
    
    for section in config.sections():
        if section.find("-deepcontact") == -1 or section in ["deepcontact","aln-deepcontact"]:
            continue        
        
        # Write section name without "-deepcontact" extension. Eg: cmmpred-deepcontact --> ccmpred
        out.write(section.replace("-deepcontact","") + ":\n")
        
        for key,value in config.items(section):
            
            
            # Convert DCA & DCB to upper case
            if key == "dca":
                key = "DCA"
            elif key == "dcb":
                key = "DCB"
                
                
            # Write keys & parameters
            out.write("        " + key + ": " + value + "\n")
        
        
        # Write about the number of threads and processes
        # n_threads
        if section in ["hhblits-deepcontact","jackhmmer-deepcontact","ccmpred-deepcontact","freecontact-deepcontact","blast-deepcontact"]:
            out.write("        n_threads: " + str(cpu) + "\n")
        # n_processes
        elif section in ["alnstats-deepcontact","makemat-deepcontact","psipred-deepcontact","psipred_pass2-deepcontact","solvpred-deepcontact","hhmake-deepcontact"]:
            out.write("        n_processes: " + str(cpu) + "\n")


        out.write("\n")
            
    out.close()




if __name__ == '__main__':

    
    import os
    import argparse
    import configparser
    import math
    
    app = argparse.ArgumentParser(description=__description__)
    app.add_argument("cpu",type=int,default=1,help="Number of cpus. Default=1")
    args = app.parse_args()
    

    # Configuration file
    bin_dir = os.path.dirname(os.path.abspath(__file__)) 
    yaml = "{}/../configuration.yaml".format(bin_dir)
    config = configparser.ConfigParser()
    config.read(yaml)
    
    
    # Output file: default.yaml
    deepcontact_dir = config.get("deepcontact","path")
    default_yaml = deepcontact_dir + "/default.yaml"
  
    print("Writing "  + default_yaml)
    
    __write_deepcontact_yaml(default_yaml,config,args.cpu)
    
            
    
