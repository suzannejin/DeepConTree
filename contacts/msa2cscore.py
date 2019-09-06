#!/usr/bin/env python3


#
## The original script comes from https://github.com/CarlosLElorduy/msa2cscore
## Usage: python msa2cscore.py <fasta with pdb id> <contactbench> <output> <pscore> <cscore> <yes/no considerWeights>
#


import sys, os
from os import listdir

# Checks if num is a number
def is_number(num):
    try:
        float(num)
        return True
    except:
        return False

# ==== PROTEIN NAME FROM PATH ====
# Returns name of the protein from its path
#   "/.../protname.contactbench" --> "protname"
def getProteinNamePath(protein_filef):
    protein_namef=''
    prot_posf,f0=-1,False
    while protein_filef[prot_posf]!='/':
        if f0:
            protein_namef = protein_filef[prot_posf]+protein_namef
        elif protein_filef[prot_posf]=='.':
            f0=True
        prot_posf-=1
    return (protein_namef,prot_posf)


# ==== PROTEIN NAME FROM ".contactbench" FILE ====
# Returns name of the protein from its contact file (NOT PATH)
#   "protname.contactbench" --> "protname"
def getProteinNameFile(protein_contbenchf):
    prot_posff=-1
    while protein_contbenchf[prot_posff]!='.':
        prot_posff-=1
    protein_nameff = protein_contbenchf[:prot_posff]
    return protein_nameff


# ==== FAMILY DIRECTORIES + NAME ====
# Returns directory where the family directory is located and the name of the family
#   "/path/to/FAMILY/" --> ("/path/to/","FAMILY")
def getFamilyDirName(sequences_pathf):
    seq_pos=-1
    while sequences_pathf[seq_pos]!='/':
        seq_pos-=1
    family_dirf      = sequences_pathf[:seq_pos]+'/'         # Path of the directory where all the families are located.
    family_namef     = sequences_pathf[seq_pos+1:]
    return (family_dirf,family_namef)


# ==== READ ALIGNMENT ====
# Returns:| aligf:         dictionary relating the name of the sequence with its respective alignment.
#         | normal2aligff: for each protein, it relates each contact with the new contact according to the alignment.
#         | pfam2pdb:      dictionary relating the pfam identifier with the pdb one. If pdb identifiers are already in the alignment, it relates the pdb identifier with itself
def readAlignment(file,temp_listf,considerdifferentf):
    print('Reading alignment...')
    pfam2pdb={}
    if considerdifferent:
        try:
            with open(temp_listf,'r') as reference:
                for line in reference.readlines():
                    line=line.split()
                    pfam2pdb[line[0][1:]]=line[2]
        except:
            print('Error! Problem with',temp_listf,'. Possible errors: Non-existent file | Not a separator between identifiers (see README.md)')
            sys.exit(1)

    with open(file,"r") as file_alignment:
        aligf={}
        main_string=''
        name=''
        prev_name=''
        first=True
        for n in file_alignment.readlines():
            if n[0]=='>':
                name=n[1:-1]
                if not considerdifferentf:
                    pfam2pdb[name]=name

                if first:
                    prev_name=name
                    first=False
                else:
                    aligf[pfam2pdb[prev_name]]=main_string
                    main_string=''
                    prev_name=name
            else:
                main_string+=n[0:-1]
        aligf[pfam2pdb[prev_name]]=main_string

        normal2aligff={}
        for prot in aligf.keys():
            c_normal=0
            c_alig=0
            normal2aligff[prot]={}
            for n in aligf[prot]:
                c_alig+=1
                if n!='-':
                    c_normal+=1
                    normal2aligff[prot][c_normal]=c_alig
    print('\tAlignment obtained succesfully')
    return (aligf,normal2aligff,pfam2pdb)


# ==== SEQUENCE_CONTACT2WEIGHT ====
# Returns: | sequence_contact2pscoref: dictionary containing the contacts with their pscore. 
#          | dic_newcontactsf:         a dictionary with only the changed contacts (with the alignment).
#          | ToF_indf:                 flag of True or False to know if the output should still include the True and False column ("correct" column)
def assignSequenceContactPscore(familydirf,aligf,normal2aligf,threshold_pf,ToF_indf):
    print('Obtaining pscores higher than',threshold_pf,'...')
    sequence_contact2pscoref={}
    dic_newcontactsf={}
    for filef in listdir(familydirf):
        if filef[-13:]=='.contactbench':
            prot_namef=getProteinNameFile(filef)
            if prot_namef!='2B8MA-1':
                dic_newcontactsf[prot_namef]=set()
                sequence_contact2pscoref[prot_namef]={}
                with open(familydirf+filef,'r') as contact_filef:
                    for line in contact_filef.readlines():
                        columns=line.split()
                        if is_number(columns[0]):
                            if float(columns[3])>=threshold_pf:
                                contact=( int(columns[0]), int(columns[1]) )

                                newc=changeContactValues(contact,normal2aligf,prot_namef)

                                if ToF_indf:
                                    try:
                                        columns[4]
                                        if not (columns[4]=='TRUE' or columns[4]=='FALSE'):
                                            ToF_indf=False
                                    except:
                                        ToF_indf=False

                                if ToF_ind:
                                    newc=list(newc)
                                    if columns[4]=='TRUE':
                                        newc.append(1)
                                    else:
                                        newc.append(0)
                                    newc=tuple(newc)

                                sequence_contact2pscoref[prot_namef][newc]=columns[3]
                                dic_newcontactsf[prot_namef].add((newc[3],newc[4]))

                        else: break
    print('\tPscores obtained succesfully')
    return (sequence_contact2pscoref,dic_newcontactsf,ToF_indf)


# ==== OBTAIN WEIGHTS ====
# Returns dictionary relating the name of the sequence with its weight
#   prot2weight[protein_name] = weight
def getWeights(alig_pathf,pfam2pdbf):
    print('Calculating Voronoid weights...')
    protein2weightf={}
    max_val=-999
    temp_protein2weightf    ={}

    current_file=cwd+'/'+'weights.tempfile'
    if os.path.exists(current_file):
        os.remove(current_file)
    
    os.system('touch '+current_file)
    try:
        os.system("weight -v "+alig_pathf+' > '+current_file)
    except:
        print('Error! Problem calculating the Voronoid weights.')
        print('Possible solutions: SQUID not installed (see README.md)')
        sys.exit(1)
        
    with open(current_file,'r') as weigh_file:
        for line in weigh_file.readlines():
            line=line.split()
            if line != []:
                if line[0]=='#=GS':
                    temp_protein2weightf[ pfam2pdbf[line[1]] ]=float(line[3])
                    if float(line[3])>max_val: max_val=float(line[3])
    os.system('rm '+current_file)
    for prot in temp_protein2weightf.keys():
        protein2weightf[prot]=temp_protein2weightf[prot]/max_val # Normalized value
    print('\tWeights calculated succesfully')
    return protein2weightf


# ==== CHANGE CONTACT VALUES ====
# Changes the contact from (old1,old2) --> (seqname,old1,old2,new1,new2) according to their alignment
def changeContactValues(contactf,norm2al,protf):
    changed_contact=( protf ,contactf[0], contactf[1], norm2al[protf][contactf[0]]   ,   norm2al[protf][contactf[1]]  )
    #               ( prot ,   old1    ,    old2    ,              new1             ,               new2             )
    return changed_contact



# +--------------+
# |     MAIN     |
# +--------------+
if __name__=='__main__':
    cwd=os.getcwd()
    
    
    
    
    # alignment_fasta=cwd+'/Example/ExAlignment/files/FAMILY1_pdb.fa'
    # protein_path=cwd+'/Example/method1/FAMILY1/seq1.contactbench'
    # name_output_file=cwd+'/EXAMPLE1_pdb.file'
    # threshold_p=0.3
    # considerWeights=True
    # considerdifferent=False

  # ==== INPUTS ====
    alignment_fasta  = sys.argv[1]                          # Path of the file containing the alignment                             /.../FAMILY.correct_aln.fa
    protein_path     = sys.argv[2]                          # Path of the file containing the pscores and contacts of the proteins  /.../PROTEIN.contactbench
    name_output_file = sys.argv[3]                          # Name of the file where the cscores will be written
    #try:
        #ref_temp_list    = cwd+'/'+sys.argv[4]                      # File containing the mapping of the two types of identifiers used
        #considerdifferent= True
    #except:
        #ref_temp_list    = None
        #considerdifferent= False
    ref_temp_list = None
    considerdifferent = False
    
    threshold_p      = float(sys.argv[4])    # Threshold for the pvalue
    # Checking threshold being between 0-1
    while not 0<=threshold_p<=1:
        threshold_p  = float(input('Select threshold for pscore (0-1): '))

    threshold_c      = float(sys.argv[5])    # Threshold for the pvalue
    # Checking threshold being between 0-1
    while not 0<=threshold_c<=1:
        threshold_c  = float(input('Select threshold for cscore (0-1): '))


    #considerWeights=None                                            # Option to consider the weights of the sequences.
    #weight_input=input('Calculate and consider Voronoid weights (y/n)? ')
    #while considerWeights==None:
        #if weight_input.upper()=='Y':
            #considerWeights=True
        #elif weight_input.upper()=='N':
            #considerWeights=False
        #else:
            #weight_input=input('Consider Voronoid weights (y/n)? ')
    v_considerWeights=sys.argv[6]
    if v_considerWeights=="yes":
      considerWeights=True
    else:
      considerWeights=None

  # ==== OBTAINING DIRECTORY NAMES AND NAMES ====
    protein_name, prot_position=getProteinNamePath(protein_path)
    sequences_path     = protein_path[:prot_position]               # Path of the family, which has all the sequences inside.
    
    family_dir,family_name=getFamilyDirName(sequences_path)         # Family_dir: directory where the family is located. | Family_name: name of the family
    sequences_path    +='/'
    
    
  # =============================================


    ToF_ind=True                                                    # Flag to check if TRUE and FALSE contacts will be written in the output file

    alig,normal2alig,pfam2pdb=readAlignment(alignment_fasta,ref_temp_list,considerdifferent)    # Alig:        dictionary relating the name of the sequence with its alignment.
                                                                                                # Normal2alig: dictionary relating the original contact with the changed one according to the alignment.
    if considerWeights:
        protein2weight=getWeights(alignment_fasta,pfam2pdb)
        

    sequence_contact2pscore,dic_newcontacts,ToF_ind=assignSequenceContactPscore(sequences_path,alig,normal2alig,threshold_p,ToF_ind)
                                                                                        # Sequence_contact2pscore: dictionary relating the complete contact with its pscore.
                                                                                        # Dic_newcontacts:         dictionary relating each sequence of the family with only their changed contacts (according to their alignment)
                                                                                        # ToF_ind:                 boolean that checks if we still have to take into account the TRUE and FALSE contacts when writting the final output (if, when reading the contact files, it doesn't detect any TRUE or FALSE, the boolean will become "False")

   # ==== CALCULATING CSCORES ====
    print('Obtaining cscores per contact of',protein_name,'...')
    contact2cscore={}                                                                   # Contact2cscore: dictionary that relates the complete contact with its cscore
    for contact_list in sequence_contact2pscore[protein_name].keys():
        oldcontact=(contact_list[1],contact_list[2])
        newcontact=(contact_list[3],contact_list[4])

        currentWeights=float(0)
        totalWeights  =float(0)

        for second_file in sequence_contact2pscore.keys():
            if second_file!=protein_name:
                if second_file in alig.keys():
                    if newcontact in dic_newcontacts[second_file]:
                        if considerWeights:
                            currentWeights+=protein2weight[second_file]
                        else:
                            currentWeights+=1
                            
                    if considerWeights:
                        totalWeights+=protein2weight[second_file]
                    else:
                        totalWeights+=1
        current_cscore=currentWeights/totalWeights
        if current_cscore>=threshold_c:
            contact2cscore[contact_list]=currentWeights/totalWeights                        # Cscore calculated by the proportion of the addition of the weights of the sequences that contain the contact by the addition of the weights of all the sequences. If weights are not taken into account, all weights become 1.0 by default
    print('\tCscores obtained succesfully')
   # =============================

    if os.path.exists(name_output_file):                # If the output file already exists, it deletes it in order to write it from zero
        os.remove(name_output_file)

    print('Writting',sys.argv[3],'...')
    with open(name_output_file,'w+') as output:
        if ToF_ind:                                                                         # If the boolean/flag is True, it will add an extra column specifying if the contact is true, being 1 if "True" and 0 if "False"
            output.write('name\tleft\tright\tpscore\tcscore\tcorrect\n')
            for contact in contact2cscore.keys():
                line_write=  ( contact[0]+'\t'+
                                str(contact[1])+'\t'+str(contact[2])+'\t'+
                                str(sequence_contact2pscore[protein_name][contact])+'\t'+
                                str(contact2cscore[contact])+'\t'+
                                str(contact[5])+'\n')
                output.write(line_write)
        else:
            output.write('name\tleft\tright\tpscore\tcscore\n')
            for contact in contact2cscore.keys():
                line_write=  ( contact[0]+'\t'+
                                str(contact[1])+'\t'+str(contact[2])+'\t'+
                                str(sequence_contact2pscore[protein_name][contact])+'\t'+
                                str(contact2cscore[contact])+'\n')
                output.write(line_write)
    print('\tWritten succesfully')

