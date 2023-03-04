import math
import os
import shutil
import sys
import networkx as nx
import random
import numpy as np
import subprocess
import os.path
import glob
import time



### Global Variables
Chain_dic = {}
Representative_dic = {}
Cluster_dic = {}
RMSD_dic = {}
RP_list = []
merged_clus_id_list = []

G = nx.Graph()



### Read Cluster ID and their representatives
def read_NR_Table(output_with_org):
    f_RPtable = open("../Output/" + output_with_org, "r")

    line = f_RPtable.readline()
    while(True):
        line = f_RPtable.readline()
        if line == "":
            break
        line = line.strip('\n').split('\t')
        clus_id, rp, members, organism, macro_name, fam = int(line[0]), line[1], line[2], line[3], line[4], line[5]
        
        members = members.lstrip('[').rstrip(']').replace("'", "").replace(" ", "").split(",")
        macro_name = macro_name.lstrip('[').rstrip(']').replace('"', '').replace("'", "").replace(" ", "").split(",")
        fam = fam.lstrip('[').rstrip(']').replace('"', '').replace("'", "").replace(" ", "").replace("[", "").replace("]", "").split(",")
        
        Representative_dic[rp] = [clus_id, members, organism, macro_name, fam]
        RP_list.append(rp)
        
   

def merge_clusters():
    
    f_grp = open("Organism_merge/Merge_based_on_rmsd", "r")
    clus_id = 1

    while(True):
        line = f_grp.readline()
        if line == "":
            break

        line = f_grp.readline()
        line = line.lstrip("{").replace("}", "").replace("'", "").replace(" ", "").strip("\n")
        RNA_list = line.split(",")
        line = f_grp.readline()
       
        no_RNA = len(RNA_list)
        if no_RNA > 1: merged_clus_id_list.append(clus_id)
        member_list = []

        for RNA in RNA_list:
            cur_member = Representative_dic[RNA][1]

            for mem in cur_member:
                member_list.append(mem)
                G.add_node(mem)

        Cluster_dic[clus_id] = member_list
        clus_id += 1
                  
    f_grp.close()

    ##### Add edge between RNA chains that already belonged together in a cluster that was not merged with any other cluster:
    for key in Cluster_dic:
        if key not in merged_clus_id_list:

            clus_mem = Cluster_dic[key]
            clus_mem_len = len(clus_mem)
            
            for i in range(clus_mem_len):
                
                chain1 = clus_mem[i]
                
                for j in range(i+1, clus_mem_len):

                    chain2 = clus_mem[j]
                    G.add_edge(chain1, chain2)



def read_existing_RMSD():

    path = "Organism_merge/Mem_RMSD"

    if os.path.exists(path):
    
        f_RMSD = open("Organism_merge/Mem_RMSD", "r")
        f_RMSD.readline()
        
        while(True):
            line = f_RMSD.readline()
            if line == "":
                break

            RNA1, RNA2, RMSD, align_percen = line.strip("\n").split("\t")
            pair = RNA1 + "_" + RNA2
            RMSD_dic[pair] = [RMSD, align_percen]
            
        f_RMSD.close()




def read_res():

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


    ### Reading RNA chain sequence for each organism
    for i in range(0, len(organism_name)):

        f_open = open("Organism_chains/" + organism_name[i] + "_RNA_chains", "r")
        
        line = f_open.readline()

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = line.split('\t')
            
            pdb_id, chain, res = line[0], line[1].strip(), line[11].strip()
            chain = chain.replace("\"", "")
            pdb = pdb_id + "_" + chain
            Chain_dic[pdb] = res

          
                
    
def generate_RMSD(rmsd_threshold, align_ratio_threshold):

    f_RMSD_temp = open("Organism_merge/RMSD_temp", "w")
    f_RMSD_temp.write("RNA1\tRNA2\tRMSD\tAlign_percen\n")
    f_RMSD_temp.close()
    

    for cl in merged_clus_id_list:

        chain_members = Cluster_dic[cl]
        chain_members_len = len(chain_members)

        f_RMSD_temp = open("Organism_merge/RMSD_temp", "a")

        for i in range(chain_members_len):
            
            RNA1 = chain_members[i]
            res1 = int(Chain_dic[RNA1][0])
            pdb_id1, chain1 = RNA1.split("_")
            
            for j in range(i+1, chain_members_len):

                RNA2 = chain_members[j]
                res2 = int(Chain_dic[RNA2][0])
                pair1 = RNA1 + "_" + RNA2
                pair2 = RNA2 + "_" + RNA1

                ######### Check if RMSD already has been calculated #########
                if pair1 in RMSD_dic:

                    RMSD = RMSD_dic[pair1][0]
                    align_percen = float(RMSD_dic[pair1][1])
                    
                    
                    f_RMSD_temp.write("%s\t%s\t%s\t%s\n" % (RNA1, RNA2, RMSD, align_percen))
                   
                    if RMSD != 'No similar stacks':
                        RMSD = float(RMSD.strip("A"))
                        if RMSD < rmsd_threshold and align_percen >= align_ratio_threshold:
                            G.add_edge(RNA1, RNA2)
                    continue


                if pair2 in RMSD_dic:

                    RMSD = RMSD_dic[pair2][0]
                    align_percen = float(RMSD_dic[pair2][1])
                    
                        
                    f_RMSD_temp.write("%s\t%s\t%s\t%s\n" % (RNA1, RNA2, RMSD, align_percen))
                  
                    if RMSD != 'No similar stacks':
                        RMSD = float(RMSD.strip("A"))
                        if RMSD < rmsd_threshold and align_percen >= align_ratio_threshold:
                            G.add_edge(RNA1, RNA2)
                    continue
                ######### Done ######### 


                pdb_id2, chain2 = RNA2.split("_")
                RMSD = '0A'
                align_len = 0
                
                STAR3D_path = 'STAR3D_source_dssr/'
                p3 = subprocess.call('java -jar STAR3D.jar -o structure_output/'+pdb_id1+'_'+chain1+'_'+pdb_id2+'_'+chain2+'.aln '+pdb_id1+' '+chain1+' '+pdb_id2+' '+chain2, shell = True, cwd = STAR3D_path)
            
                if(os.path.isfile(STAR3D_path + 'structure_output/' + pdb_id1 + '_' + chain1+ '_' + pdb_id2 + '_' + chain2 + '.aln')):
                    alignment_open = open(STAR3D_path + 'structure_output/' + pdb_id1 + '_' + chain1 + '_' + pdb_id2 + '_' + chain2 + '.aln', "r")

                    for line_skip in range(0,13):
                        cur_line = alignment_open.readline()

                    cur_line = alignment_open.readline()
                    cur_line = cur_line.split()
                    align_len = int(cur_line[2])
                    cur_line = alignment_open.readline()
                    cur_line = cur_line.split()
                    RMSD = cur_line[2]
                else:
                    RMSD = 'No similar stacks'


                ### Calculate Align Percent using Number of Residues with Coordinates                    
                align_percen = 0.0
                if(res1 < res2):
                    align_percen = (align_len/res1)*100
                else:
                    align_percen = (align_len/res2)*100


                f_RMSD_temp.write("%s\t%s\t%s\t%s\n" % (RNA1, RNA2, RMSD, align_percen))
                
                #print(RNA1, "||", RNA2, " RMSD: ", RMSD)
                if RMSD != 'No similar stacks':
                    RMSD = float(RMSD.strip("A"))
                    if RMSD < rmsd_threshold and align_percen >= align_ratio_threshold:
                        G.add_edge(RNA1, RNA2)
                
        f_RMSD_temp.close()

    ######### Copy RP_RMSD_temp to RP_RMSD and delete RP_RMSD_temp #########
    shutil.copyfile('Organism_merge/RMSD_temp','Organism_merge/Mem_RMSD')
    os.remove("Organism_merge/RMSD_temp")
    ######### Done #########

    
        

def connected_graph():
    
    f_grp_open = open("Organism_merge/RNA_clusters", "w")
    connected_graph_list = list(nx.connected_components(G))

    for component in range(0, len(connected_graph_list)):
        f_grp_open.write("Cluster no: "+ str(component+1) + "\n")
        f_grp_open.write(str(connected_graph_list[component])+ "\n")

        node_list = list(connected_graph_list[component])
        degree_list = []
        
        for node in node_list:
            
            degree_list.append(G.degree(node))
                
        f_grp_open.write(str(degree_list)+ "\n")

    f_grp_open.close()

            

def align_member_struct(rmsd_threshold, align_ratio_threshold, output_with_org):

    print("\n ------- Generating structure alignement among cluster members ------- \n")
    
    start_time = time.time()
    
    ### Create log file
    path = "Organism_merge/Log_file"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)
    ###
    
    read_NR_Table(output_with_org)
    read_res()
    merge_clusters()
    read_existing_RMSD()
    generate_RMSD(rmsd_threshold, align_ratio_threshold)
    connected_graph()

    run_time = round((time.time() - start_time)/3600, 2)
    
    log_file = open('Organism_merge/Log_file/RNA_clusters_runtime', "w")
    log_file.write("--- %s hours ---" % (str(run_time)))
    log_file.close()





        

