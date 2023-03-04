import math
import os
import shutil
import sys
sys.path.append('src/')


def generate_file_for_identity_score_calculation():

    print("\n ------- Generating necessary files for calculating RNA sequence alignment ------- \n")

    f_organism_name = open("Organism_list/Current_Organism_list","r")
    organism_name = []
    organism_count = 0


    while(True):
        line = f_organism_name.readline()
        
        if line == "":
            break
        
        line = line.strip("\n")
        organism = line.replace(" ", "_")
        organism_name.append(organism)
        organism_count += 1

    f_organism_name.close()

  
    path_identity = "Organism_identity_score/identity_score_list_temp"

    if(os.path.isdir(path_identity)):
        shutil.rmtree(path_identity)    
    os.mkdir(path_identity)

    path_nodelist = "Organism_group/node_list_temp"

    if(os.path.isdir(path_nodelist)):
        shutil.rmtree(path_nodelist)    
    os.mkdir(path_nodelist)

    path_nodelist = "Organism_group/group_based_on_identity_score"

    if(os.path.isdir(path_nodelist)):
        shutil.rmtree(path_nodelist)    
    os.mkdir(path_nodelist)
    ### Done Creating Temporary Folders ###



    ### Creating Temporary Folders for RNA chains with 100 identiy score###
    path_100_identity = "Organism_group/" + "identity_100_list"

    if(os.path.isdir(path_100_identity)):
        shutil.rmtree(path_100_identity)    
    os.mkdir(path_100_identity)
    ### Done Creating Temporary Folders ###



    for i in range(0, len(organism_name)):

        f_open = open("Organism_chains/" + organism_name[i] + "_RNA_chains", "r")
        
        chain_list = {}
        line = f_open.readline()
        key_val = 0

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = line.split('\t')
           
            pdb_id, chain, RNA_name, sequence = line[0], line[1].strip(), line[9], line[13]
            chain = chain.replace("\"", "")
            chain_list[key_val] = [pdb_id, chain, RNA_name, sequence]
            key_val += 1
        
        f_open.close()
        

        ### Writing List of Nodes to temp folder 'node_list_temp' ###
        f_nodelist_out = open("Organism_group/node_list_temp/" + organism_name[i] + "_nodelist", "w")

        
        for key in chain_list:
            chain_id = chain_list[key][1]
            chain_id = chain_id.replace(" ", "")
            chain_id_list = chain_id.split(",")

            for ch in range(0, len(chain_id_list)):
                node = chain_list[key][0] + "_" + chain_id_list[ch]
                f_nodelist_out.write(node+"\n")
      
        f_nodelist_out.close()
        ### Done Writing ###
    

        


