#!/bin/bash

#
# This script is employed to configure the working environment,
# so that the corresponding programs {DeepCov, DeepConPred2, DeepContact, Spider3, etc}
# would be computed with the proper paths and dependencies, etc.
#
# Please configurate each function according to your working environment so that the corresponding programs can be used properly.
# The main script 'pipeline.py' and other scripts will set the environment by calling this script before running the programs.
#
# The usage is: bash set_envs.sh method 
# For example: bash set_envs.sh deepcontact
#



function set_env_deepcov () {
  
    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/users/cn/sjin/programs/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else 
        if [ -f "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh" ]; then
            . "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh"  
        else
            export PATH="/users/cn/sjin/programs/miniconda3/bin:$PATH"  
        fi
    fi
    unset __conda_setup
    # <<< conda initialize <<<
  
  
    conda activate deepcov-env
  
}


function set_env_deepconpred () {
  
  module load Python/2.7.14-foss-2018a   # Tensorflow 1.6.0 is installed using this version of python in our cluster

}


function set_env_deepcontact () {

    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/users/cn/sjin/programs/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else 
        if [ -f "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh" ]; then
            . "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh"  
        else
            export PATH="/users/cn/sjin/programs/miniconda3/bin:$PATH"  
        fi
    fi
    unset __conda_setup
    # <<< conda initialize <<<
    
    
    # Load modules required by FreeContact
    module load BLAST/2.2.26-Linux_x86_64
    module load Boost/1.63.0-foss-2017a
    module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
    
    # Other variables
    conda activate deepcontact-env
    export BLASTDB=/db/ncbi/201903/blastdb/db
    export BLASTMAT=/software/as/el7.2/EasyBuild/CRG/software/BLAST/2.2.26-Linux_x86_64/data
    export PYTHONPATH=/users/cn/sjin/programs/DeepContact/deepcontact  # This enables DeepContact to import python modules from other of their scripts
    
}


function set_env_spider3 () {

    # Note that Spider3 uses Tensorflow 0.12, whereas DeepConPred2 uses Tensorflow 1.6.0
    export PYTHONPATH=/users/cn/sjin/.local2

}


function set_env_tree () {

    # >>> conda initialize >>>
    # !! Contents within this block are managed by 'conda init' !!
    __conda_setup="$('/users/cn/sjin/programs/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else 
        if [ -f "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh" ]; then
            . "/users/cn/sjin/programs/miniconda3/etc/profile.d/conda.sh"  
        else
            export PATH="/users/cn/sjin/programs/miniconda3/bin:$PATH"  
        fi
    fi
    unset __conda_setup
    # <<< conda initialize <<<
    
    
    conda activate trmsd-env
    
    module load R/3.6.0-foss-2019a
  
}


method=$1

[[ $method == "deepcov" ]] && set_env_deepcov
[[ $method == "deepconpred" ]] && set_env_deepconpred
[[ $method == "deepcontact" ]] && set_env_deepcontact
[[ $method == "spider" ]] && set_env_spider3
[[ $method == "tree" ]] && set_env_tree

export PYTHONPATH=~/DeepConTree:$PYTHONPATH  
