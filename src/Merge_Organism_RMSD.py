import math
import os
import shutil
import sys
##from Bio.Seq import Seq 
##from Bio.pairwise2 import format_alignment
##from Bio import pairwise2
import networkx as nx
import random
import numpy as np
import subprocess
import os.path
import glob
import time
import shutil


### Global Variables
Representative_dic = {}
RP_list = []
##organism_division = ['ABC', 'DEFGH', 'IJKLM', 'NOPQR', 'STUV', 'WXYZ']
G = nx.Graph()
RMSD_dic = {}


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
        G.add_node(rp)




def read_existing_RMSD():
    
    path = "Organism_merge/RP_RMSD"

    if os.path.exists(path):

        f_RMSD = open("Organism_merge/RP_RMSD", "r")
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

            if pdb in Representative_dic:
                Representative_dic[pdb].append(res)




                

def generate_RMSD(rmsd_threshold, align_ratio_threshold):
    
    f_grp = open("Organism_merge/Merge_based_on_identity_score", "r")

    f_RMSD_temp = open("Organism_merge/RP_RMSD_temp", "w")
    f_RMSD_temp.write("RNA1\tRNA2\tRMSD\tAlign_percent\n")
    f_RMSD_temp.close()
        
    
    while(True):
        line = f_grp.readline()
        if line == "":
            break
        line = f_grp.readline()
        line = line.lstrip("{").replace("}", "").replace("'", "").replace(" ", "").strip("\n")
        RNA_list = line.split(",")
       
        no_RNA = len(RNA_list)

        if no_RNA > 1:

            for RNA in RNA_list:
                pdb_id, chain = RNA.split("_")
                #print(RNA)
                STAR3D_path = 'STAR3D_source_dssr/'
                p1 = subprocess.call('java -cp STAR3D.jar Preprocess ' + pdb_id + ' ' + chain, shell = True, cwd = STAR3D_path)

            f_RMSD_temp = open("Organism_merge/RP_RMSD_temp", "a")

            for i in range(no_RNA):

                RNA1 = RNA_list[i]
                pdb_id1, chain1 = RNA1.split("_")
                res1 = int(Representative_dic[RNA1][5])
                
                for j in range(i+1, no_RNA):
                    
                    RNA2 = RNA_list[j]
                    res2 = int(Representative_dic[RNA2][5])
                    pair1 = RNA1 + "_" + RNA2
                    pair2 = RNA2 + "_" + RNA1

                    ######### Check if RMSD already has been calculated #########
                    if pair1 in RMSD_dic:

                        RMSD = RMSD_dic[pair1][0]
                        align_percen = float(RMSD_dic[pair1][1])
                        #print("RMSD Exists: ", RMSD)
                        f_RMSD_temp.write("%s\t%s\t%s\t%s\n" % (RNA1, RNA2, RMSD, align_percen))
                        if RMSD != 'No similar stacks':
                            RMSD = float(RMSD.strip("A"))
                            if RMSD < rmsd_threshold and align_percen >= align_ratio_threshold:
                                G.add_edge(RNA1, RNA2)
                        continue


                    if pair2 in RMSD_dic:

                        RMSD = RMSD_dic[pair2][0]
                        align_percen = float(RMSD_dic[pair2][1])
                        #print("RMSD Exists: ", RMSD)
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

                    #print("RMSD: ", RMSD)
                    


                    ### Calculate Align Percent using Number of Residues with Coordinates                    
                    align_percen = 0.0
                    if(res1 < res2):
                        align_percen = (align_len/res1)*100
                    else:
                        align_percen = (align_len/res2)*100

                    f_RMSD_temp.write("%s\t%s\t%s\t%s\n" % (RNA1, RNA2, RMSD, align_percen))

                        
                    if RMSD != 'No similar stacks':
                        RMSD = float(RMSD.strip("A"))
                        if RMSD < rmsd_threshold and align_percen >= align_ratio_threshold:
                            G.add_edge(RNA1, RNA2)
                            
            f_RMSD_temp.close()
            
    f_grp.close()
    

    ######### Copy RP_RMSD_temp to RP_RMSD and delete RP_RMSD_temp #########
    shutil.copyfile('Organism_merge/RP_RMSD_temp','Organism_merge/RP_RMSD')
    os.remove("Organism_merge/RP_RMSD_temp")
    ######### Done ######### 




def connected_graph():
    
    f_grp_open = open("Organism_merge/Merge_based_on_rmsd", "w")
    connected_graph_list = list(nx.connected_components(G))

    for component in range(0, len(connected_graph_list)):
        f_grp_open.write("Cluster no: "+ str(component+1) + "\n")
        f_grp_open.write(str(connected_graph_list[component])+ "\n")

        cluster_no_list = []
        for RP in connected_graph_list[component]:
            clus_id = Representative_dic[RP][0]
            cluster_no_list.append(clus_id)

        f_grp_open.write(str(cluster_no_list)+ "\n")

    f_grp_open.close()

            




def align_RP_struct(rmsd_threshold, align_ratio_threshold, output_with_org):

    print("\n ------- Generating structure alignement among cluster representatives ------- \n")

    start_time = time.time()

    ### Create log file
    path = "Organism_merge/Log_file"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)
    ###

    
    read_NR_Table(output_with_org)
    read_res()
    read_existing_RMSD()
    generate_RMSD(rmsd_threshold, align_ratio_threshold)
    connected_graph()

    run_time = round((time.time() - start_time)/3600, 2)
    
    log_file = open('Organism_merge/Log_file/RMSD_runtime', "w")
    log_file.write("--- %s hours ---" % (str(run_time)))
    log_file.close()







        

