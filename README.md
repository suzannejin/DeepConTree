# DeepConTree pipeline

## Installation
1.Install the following dependencies:
  - DeepCov
  - DeepConPred2
  - DeepContact
  - Xplor-NIH 
  - T-coffee 
2. Install R packages:
  - phangorn
  - geiger
  - xlsx
3. Python modules:
  - argparse
  
  Other dependencies, such as Tensorflow, Theano, Lasagne, etc. are specified by each of the software above. They should be installed as well.

4. Modify configuration.yaml and set_envs.sh files with the correct paths
5. Set environment variable `export PYTHONPATH=/path/to/DeepConTree:$PYTHONPATH`

## Usage
This pipeline can be used for the prediction of **protein contacts**, the prediction of protein **3D structures** through a contact-driven protein folding approach, and/or the reconstruction of **structure-based trees**.
You need to specify the following in order to run the pipeline:
1. `-seq`    Multi-fasta file (with suffix .fa) with all the sequences on which you want to predict the contacts/structures or create 3D-trees. The sequences can be named with either PDB ID or UniProt ID, as long as you provide the reference template file.
2. `-pdbf`   The folder where the pdb files (with suffix .pdb) are stored. The pipeline only considers the PDB files of the sequences specified in the multi-fasta file.
3. `-ref`    The reference template file. Each line from the reference file should contain three columns: 1) the UniProt identifier (sequence name), 2) a separator and 3) the PDB identifier. For example: >AMCY_PARDE/43-131 _P_ 1MDAA-1 .
4. `-o`      Output directory. All the files generated will be stored here. The pipeline will also copy the original fasta & pdb files to here for easy usage. If you don't want to keep these files in the output directory, you can just add some lines at the end of the pipeline to remove these files.

You have to specify what do you want to do:

`-actions`   Three options are provided: **contacts, simulations, trees**. You have to specify at least one action to run the pipeline. 

Since these processes are computationally demanding, you can use several CPUs:

`-cpu`       Number of CPUs. Default is 1.   


### Prediction of contacts 
If you want to predict contacts, please specify `-action contacts`
and the contact predictors that you want to use:

`-method`    Three options are provided: **deepcov, deepconpred, deepcontact**. 

For example, if you want to predict contacts using deepconpred and deepcontact, you can do:
```bash
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions contacts -method deepconpred deepcontact
```

### Prediction of 3D structures
If you want to predict structures, please specify `-action simulations`.
You can run the pipeline in a straightforward way by predicting both contacts and structures. Optionally, you can specify the filtering thresholds you want. If the thresholds are not specified, default values will be used (check `python path/to/pipeline.py -h`).
```bash
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions contacts simulations -method deepconpred deepcontact -pdeepconpred 0.35 -pdeepcontact 0.95
```
Another option (recommended) is:
1. First, predict the **contacts for all the proteins** you want.
2. Check the predictions and **decide the probability threshold** based on precision-recall metrics, etc. If you want you can also set a threshold based on conservation score.
3. **Run simulations** using the already predicted contacts and the selected filtering thresholds. In this case, you have to:
    a) Provide a list `-con_list path/to/con_list` with three columns: /absolute/path/to/contact_folder folder_name probability_threshold. In this way, you can provide multiple folders, where contacts predicted by different methods, etc. are stored. The resulting structures will be stored in output_dir/simulations/folder_name. As you might want to use a different threshold of probability score for each subset of contacts, this should be specified in the third column. 
    b) If you want to filter the contacts according to both pscore and cscore, you have to specify -cscore as well. In this case, you need to provide the alignment that will be considered for the calculation of conservation score. Make sure that the pdb files of all the proteins in that alignment are stored in the pdb folder you provide, and the contact files of all the proteins are provided. Note that the -seq multifasta file that you provide can contain only few proteins that you are interested on predicting the structure, but the alignment might contain other homologous sequences that are not in the -seq multifasta file.
    c) Optionally, you can provide the folder where the SSE (secondary structure element) predictions are stored. If not, the pipeline will use Spider3 to predict the SSEs before running the simulations. The SSE files (with suffix .sse) should have the same format as the output of Spider3.
For example:
```bash
# 1 Predict contacts
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions contacts -method deepconpred deepcontact 
# 2 Run simulations once you have decided the thresholds
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions simulations -con_list path/to/con_list -ssef path/to/ssef
# 3 Run simulations once you have decided the thresholds, and you also want to apply a threshold based on the conservation score
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions simulations -con_list path/to/con_list -ssef path/to/ssef -cscore cscore -c_aln path/to/alignment_for_cscore.fa 
```
Since simulations need **large computational requirements**. It is recommended to run simulations for few proteins each time (using a multifasta file that only contains the sequences of these few proteins), and specifying the same output directory so that you can have all the results in the same folder. 

### Reconstruction of structure-based trees
If you want to reconstruct the trees based on structures, please specify `-action trees`. 
Moreover, the normalized Robinson-Foulds metric is estimated for each pair of trees. An excel file with a matrix of the pairwise normalized RF score will be generated. In addition, a text file with the normalized RF score obtained when comparing the predicted tree to the one created based on experimental structures will be created. 
As in the previous case, you can run the pipeline in a straightforward way:
```bash
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions contacts simulations trees -method deepconpred deepcontact -pdeepconpred 0.35 -pdeepcontact 0.95
```
Or run the pipeline once you have all the structures predicted. In this case, you have to provide a list `-str_list path/to/str_list` with two columns: /absolute/path/to/structure_folder folder_name. In this way, a tree will created for each folder provided.
```bash
python path/to/pipeline.py -seq path/to/multifasta.fa -pdbf path/to/pdbf -ref path/to/ref -o path/to/output_directory -cpu 16 -actions trees -str_list path/to/str_list
```



