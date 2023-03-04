import subprocess
import os.path
import glob
import time
import sys
sys.path.append('src/')


def generate_RMSD_align():

    print("\n ------- Generate RNA 3D strcture alignment ------- \n")
    
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
       
        ######################### Reading list of RNA chain pairs having identity score 100
        organism_100_identical_file_path = 'Organism_group/identity_100_list/' + organism_name[i] + '_identity_100'
        f_100_identity = open(organism_100_identical_file_path, "r")

        identity_100_dict = {}
        while(True):
            line = f_100_identity.readline()
            if(line == ''):
                break

            node1, node2 = line.strip().split('\t')
            if(node1 in identity_100_dict.keys()):
                identity_100_dict[node1].append(node2)
            else:
                identity_100_dict[node1] = [node2]


        ######################### Reading length of RNA chains
        chain_len_file_path = 'Organism_chains/' + organism_name[i] + '_RNA_chains'
        chain_len = open(chain_len_file_path, "r")
        
        line = chain_len.readline()
        seq_len_dic = {}
        
        while(True):
            line = chain_len.readline()
            if(line == ''):
                break

            line = line.strip().split('\t')
            PDB, chain, length = line[0], line[1], line[4]
            chain = chain.replace(" ", "")
            chain = chain.strip('"')
            chain_list = chain.split(",")
            for j in range(0, len(chain_list)):
                pdb_id = PDB + '_' + chain_list[j]
                seq_len_dic[pdb_id] = [length]


        ######################### Reading num of residues with coordinates of RNA chains
        chain_len_file_path = 'Organism_chains/' + organism_name[i] + '_RNA_chains'
        chain_len = open(chain_len_file_path, "r")

        res_len_dic = {}
        
        while(True):
            line = chain_len.readline()
            if(line == ''):
                break

            line = line.strip().split('\t')
            PDB, chain, length = line[0], line[1], line[11]
            chain = chain.replace(" ", "")
            chain = chain.strip('"')
            chain_list = chain.split(",")
            for j in range(0, len(chain_list)):
                pdb_id = PDB + '_' + chain_list[j]
                res_len_dic[pdb_id] = length

    
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
            f_out = open('Organism_RMSD/rmsd_cluster_tmp/' + organism_name[i] + "_cluster_" + str(list_count), "w")
            log_file = open('Organism_group/Log_file/Align_' + organism_name[i], "w")
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


            ######################### If cluster contains only 1 chain, then just print the chain name in f_out file, also consider this as one cluster and store in the cluster_based_on_RMSD_temp foler
            if(len_tot_chain_list == 1):
                f_out.write("Only chain in cluster: " + str(Total_chain_list[0][0]) + "_" + str(Total_chain_list[0][1]) + "\n")
            else:
                f_out.write("PDB_id1\tchain1\tSeq1\tPDB_id2\tchain2\tSeq2\tRMSD_Value\tAlignment_len\tRes_cor1\tRes_cor2\n")


            ######################### Running 3D structural alignment using STAR3D_dssr
            for key in Total_chain_list:

                current_chain = Total_chain_list[key]
                pdb_id1 = current_chain[0]
                chain1 = current_chain[1]
                comparing_key = key

                for key2 in range(comparing_key + 1, len_tot_chain_list):
                    
                    comparing_chain = Total_chain_list[key2]
                    pdb_id2 = comparing_chain[0]         
                    chain2 = comparing_chain[1]

                    id1 = pdb_id1 + '_' + chain1
                    id2 = pdb_id2 + '_' + chain2
                    seq1 = 0
                    seq2 = 0

                    for key in seq_len_dic:
                        if(key == id1):
                            seq1 = seq_len_dic[key]
                            seq1 = seq1[0]
                        if(key == id2):
                            seq2 = seq_len_dic[key]
                            seq2 = seq2[0]
                    
                    print("Aligning Pair: " + pdb_id1 + '_' + chain1 + "||" + pdb_id2 + '_' + chain2)
                   
                    node1 = pdb_id1 + "_" + chain1
                    node2 = pdb_id2 + "_" + chain2
                    res_cor1 = res_len_dic[node1]
                    res_cor2 = res_len_dic[node2]

                    #RMSD_calc_path = 'RNAMotifContrast/src/scripts/'
                    RMSD_calc_path = 'Identical_Align/src/scripts/'

                    if node1 in identity_100_dict.keys():
                        if node2 in identity_100_dict[node1]:
                            #print(node1 + " has 100 identity with " + node2 + "\n")
                            p3 = subprocess.check_output('python3 -W ignore RMSD_calculation_pipeline.py '+ node1 + ' ' + node2, shell = True, cwd = RMSD_calc_path)
                            p3 = p3.decode("utf-8")
                            p3 = p3.strip().split(" ")
                            RMSD = p3[2]+'A'
                            align_len = p3[6]
                            f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id1, chain1, seq1, pdb_id2, chain2, seq2, RMSD, align_len, res_cor1, res_cor2))
                            print("RMSD: ", RMSD) 
                            continue
                            
                    if node2 in identity_100_dict.keys():
                        if node1 in identity_100_dict[node2]:
                            #print(node1 + " has 100 identity with " + node2 + "\n")
                            p3 = subprocess.check_output('python3 -W ignore RMSD_calculation_pipeline.py '+ node1 + ' ' + node2, shell = True, cwd = RMSD_calc_path)
                            p3 = p3.decode("utf-8")
                            p3 = p3.strip().split(" ")
                            RMSD = p3[2]+'A'
                            align_len = p3[6]
                            f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id1, chain1, seq1, pdb_id2, chain2, seq2, RMSD, align_len, res_cor1, res_cor2))
                            print("RMSD: ", RMSD) 
                            continue
                        
                    STAR3D_path = 'STAR3D_source_dssr/'
                    p3 = subprocess.call('java -jar STAR3D.jar -o structure_output/'+pdb_id1+'_'+chain1+'_'+pdb_id2+'_'+chain2+'.aln '+pdb_id1+' '+chain1+' '+pdb_id2+' '+chain2, shell = True, cwd = STAR3D_path)
                    
                    STAR3D_path = 'STAR3D_source_dssr/'
                    RMSD = '0A'
                    align_len = 0
                    if(os.path.isfile(STAR3D_path + 'structure_output/' + pdb_id1 + '_' + chain1+ '_' + pdb_id2 + '_' + chain2 + '.aln')):
                        alignment_open = open(STAR3D_path + 'structure_output/' + pdb_id1 + '_' + chain1 + '_' + pdb_id2 + '_' + chain2 + '.aln', "r")

                        for line_skip in range(0,13):
                            cur_line = alignment_open.readline()

                        cur_line = alignment_open.readline()
                        cur_line = cur_line.split()
                        align_len = cur_line[2]
                        cur_line = alignment_open.readline()
                        cur_line = cur_line.split()
                        RMSD = cur_line[2]
                    else:
                        RMSD = 'No similar stacks'

                    print("RMSD: ", RMSD)   
                    f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id1, chain1, seq1, pdb_id2, chain2, seq2, RMSD, align_len, res_cor1, res_cor2))            
                        


            f_all_chains_in_cluster.close()
            f_out.close()
            run_time = round((time.time() - start_time)/3600,2)
            log_file.write("--- %s hours ---" % (str(run_time)))
