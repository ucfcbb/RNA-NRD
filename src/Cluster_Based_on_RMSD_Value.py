import math
import os
import shutil
import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
sys.path.append('src/')


def generate_cluster_based_on_RMSD(rmsd_threshold, align_ratio_threshold):

    print("\n ------- Generate RNA chain clusters based on RNA 3D strcture similarity ------- \n")

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



    for i in range(0, len(organism_name)):

        ######################### For each organism, read number of groups based on identity score similarity
        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        no_of_input_files = 0
        cluster_length_list = []

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line_list = line.split(',')
            cluster_length_list.append(len(line_list))
            no_of_input_files += 1
        
       
        for input_file_ind in range(0, no_of_input_files):
            
            G = nx.Graph()
            
            ### Read list of nodes from sequence identity based group divisions
            f_open_cluster = open("Organism_group/input_file_list_temp/" + organism_name[i] + "_cluster_" + str(input_file_ind+1), "r")

            while(True):
                line = f_open_cluster.readline()
                line = line.strip()
                if(line == ''):
                    break
                
                G.add_node(line)

            f_open_cluster.close() 
            no_of_chains = cluster_length_list[input_file_ind]

            ### Read RMSD values from alignment output files
            f_open_RMSD = open("Organism_RMSD/rmsd_cluster_tmp/" + organism_name[i] + "_cluster_" + str(input_file_ind + 1), "r")
            line2 = f_open_RMSD.readline()
            
            while(True):
                line2 = f_open_RMSD.readline()
                
                if(line2 == ''):
                    break

                line2 = line2.strip()
                line2 = line2.split('\t')
                pdb_id1 = line2[0]
                chain1 = line2[1]
                seq_len1 = int(line2[2])
                pdb_id2 = line2[3]
                chain2 = line2[4]
                seq_len2 = int(line2[5])
                RMSD = line2[6]
                align_len = int(line2[7])
                res_cor1 = int(line2[8])
                res_cor2 = int(line2[9])

                    
                ### Calculate Align Percent using Number of Residues with Coordinates                    
                align_percen = 0.0
                if(res_cor1 < res_cor2):
                    align_percen = (align_len/res_cor1)*100
                else:
                    align_percen = (align_len/res_cor2)*100

                  
                if(RMSD != "No similar stacks"):
                    RMSD = RMSD.rstrip('A')
                    if(float(RMSD) < rmsd_threshold and align_percen >= align_ratio_threshold):

                        pdb1 = pdb_id1 + "_" + chain1
                        pdb2 = pdb_id2 + "_" + chain2
                        G.add_edge(pdb1, pdb2)
                        
            f_open_RMSD.close()
                

            f_open_cluster_list = open("Organism_group/cluster_based_on_RMSD_temp/" + organism_name[i] + "_cluster_" + str(input_file_ind + 1), "w")
            connected_graph_list = list(nx.connected_components(G))
           
            degree_list = []

            for component in range(0, len(connected_graph_list)):
                f_open_cluster_list.write("Group no: "+str(component) + "\n")
                f_open_cluster_list.write(str(connected_graph_list[component])+ "\n")
                node_list = list(connected_graph_list[component])
                degree_list = []
                
                for node in node_list:
                    degree_list.append(G.degree(node))
                    
                f_open_cluster_list.write(str(degree_list)+ "\n")

            f_open_cluster_list.close()
                
        f_open.close()
        


