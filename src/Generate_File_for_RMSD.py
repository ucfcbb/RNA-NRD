import math
import os
import shutil
import sys
sys.path.append('src/')


def generate_file_for_RMSD_calculation():

    print("\n ------- Generate necessary files for aligning RNA 3D strctures ------- \n")

    f_organism_name = open("Organism_list/Current_Organism_list", "r")
    organism_name = []
  
    while(True):
        line = f_organism_name.readline()
        
        if line == "":
            break
        
        line = line.strip("\n")
        organism = line.replace(" ", "_")
        organism_name.append(organism)
      
    f_organism_name.close()
   

    ### Creating Temporary Folders for Calculating RMSD from Complete Graph Parallely###
    path = "Organism_group/input_file_list_temp"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)


    path_rmsd = "Organism_RMSD/rmsd_cluster_tmp"

    if(os.path.isdir(path_rmsd)):
        shutil.rmtree(path_rmsd)    
    os.mkdir(path_rmsd)


    path_rmsd_cluster = "Organism_group/cluster_based_on_RMSD_temp"

    if(os.path.isdir(path_rmsd_cluster)):
        shutil.rmtree(path_rmsd_cluster)    
    os.mkdir(path_rmsd_cluster)


    path_log = "Organism_group/Log_file"

    if(os.path.isdir(path_log)):
        shutil.rmtree(path_log)    
    os.mkdir(path_log)
    ### Done Creating Temporary Folders ###



    for i in range(0, len(organism_name)):

        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        
        chain_list = {}
        key_val = 0
        list_count = 0

        # For each organism, write the nodes in each group to different files
        while(True):
            
            line = f_open.readline()
            
            if(line == ''):
                break
            
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line_list = line.split(',')
            list_count += 1

            f_open_input_file = open("Organism_group/input_file_list_temp/" + organism_name[i] + "_cluster_" + str(list_count), "w")
            
            for ch in range(0, len(line_list)):
                line_list[ch] = line_list[ch].strip().lstrip("'").rstrip("'")
                f_open_input_file.write(line_list[ch] + "\n")   

            f_open_input_file.close()

        f_open.close()


