import subprocess
import os.path
import glob
import time
import sys
sys.path.append('src/')


def generate_RMSD_preprocess():

    print("\n ------- Preprocessing PDB files for generating 3D strcutre alignment ------- \n")

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
    


    for i in range(0, len(organism_name)):

        ######################### For each organism, read number of groups based on identity score similarity
        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        
        chain_list = {}
        key_val = 0
        list_count = 0

        
        while(True):
            
            line = f_open.readline()
            
            if(line == ''):
                break
            
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line_list = line.split(',')
            list_count += 1


            organism_cluster_file_path = 'Organism_group/input_file_list_temp/' + organism_name[i] + "_cluster_" + str(list_count)
            f_all_chains_in_cluster = open(organism_cluster_file_path, "r")
            log_file = open('Organism_group/Log_file/Preprocess_' + organism_name[i], "w")
            start_time = time.time()


            ######################### List of all chains in the cluster of an organism
            Total_chain_list = {}
            tot_key_val = 0
            while(True):
                
                line = f_all_chains_in_cluster.readline()
                if(line == ''):
                    break
                line = line.split('_')
                pdb_id, chain,  = line[0], line[1].strip()
                Total_chain_list[tot_key_val] = [pdb_id, chain]
                
                tot_key_val += 1

            len_tot_chain_list = len(Total_chain_list)

            
            ######################### Preprocessing pdb/cif files using STAR3D_dssr
            for main_key in Total_chain_list:
                
                comparing_chain = Total_chain_list[main_key]
                pdb_id = comparing_chain[0]         
                chain = comparing_chain[1]
                print("Preprocessing " + pdb_id + ' ' + chain)
                
                STAR3D_path = 'STAR3D_source_dssr/'
                p1 = subprocess.call('java -cp STAR3D.jar Preprocess ' + pdb_id + ' ' + chain, shell = True, cwd = STAR3D_path)


            f_all_chains_in_cluster.close()
            run_time = round((time.time() - start_time)/3600,2)
            log_file.write("--- %s hours ---" % (str(run_time)))
            log_file.close()
        

        f_open.close()
