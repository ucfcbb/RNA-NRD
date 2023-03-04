import math
import os
import shutil
import sys
from Bio.Seq import Seq 
##from Bio.pairwise2 import format_alignment
##from Bio import pairwise2
from Bio.Align import PairwiseAligner
import networkx as nx
import random
import numpy as np
import time


### Global Variables
Representative_dic = {}
RP_list = []
##organism_division = ['ABC', 'DEFGH', 'IJKLM', 'NOPQR', 'STUV', 'WXYZ']
G = nx.Graph()
SI_dic = {}


### Read Cluster ID and their representatives
def read_NR_Table(output_with_org):

    #print(G.nodes())
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

    #print(Representative_dic)
    #print(RP_list)
   


def read_existing_SI():

    path = "Organism_merge/RP_SI"

    if os.path.exists(path):
    
        f_SI = open("Organism_merge/RP_SI", "r")
        f_SI.readline()
        
        while(True):
            line = f_SI.readline()
            if line == "":
                break

            RNA1, RNA2, SI = line.strip("\n").split("\t")
            pair = RNA1 + "_" + RNA2
            SI_dic[pair] = SI
            
        f_SI.close()


    

def read_RP_sequence():

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
            
            pdb_id, chain, sequence = line[0], line[1].strip(), line[13].strip("\n")
            chain = chain.replace("\"", "")
            pdb = pdb_id + "_" + chain

            if pdb in Representative_dic:
                Representative_dic[pdb].append(sequence)

                    





def align_sequence(identity_threshold):


    f_SI_temp = open("Organism_merge/RP_SI_temp", "w")
    f_SI_temp.write("RNA1\tRNA2\tSI\n")
    no_RP = len(RP_list)

    for i in range(no_RP):
        RP1 = RP_list[i]
        seq1 = Representative_dic[RP1][5]
        seq1_len = len(seq1)
        
        for j in range(i+1, no_RP):
            RP2 = RP_list[j]
            seq2 = Representative_dic[RP2][5]
            seq2_len = len(seq2)

            pair1 = RP1 + "_" + RP2
            pair2 = RP2 + "_" + RP1
           
            if((seq1_len < 2*seq2_len) and (seq2_len < 2*seq1_len)):

                ######### Check if Sequence identity already has been calculated #########
                if pair1 in SI_dic:

                    SI = SI_dic[pair1]
                    #print("Sequence identity exists: ", SI)
                    f_SI_temp.write("%s\t%s\t%s\n" % (RP1, RP2, SI))
                    SI = float(SI)
                    if SI >= identity_threshold:
                        G.add_edge(RP1, RP2)   
                    continue


                if pair2 in SI_dic:

                    SI = SI_dic[pair2]
                    #print("Sequence identity exists: ", SI)
                    f_SI_temp.write("%s\t%s\t%s\n" % (RP1, RP2, SI))
                    SI = float(SI)
                    if SI >= identity_threshold:
                        G.add_edge(RP1, RP2)   
                    continue
                ######### Done ######### 

                
             
                identity = 0.0
                #alignments = pairwise2.align.globalxx(seq1, seq2)
                aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=0)
                alignments = aligner.align(seq1, seq2)
                
                for alignment in alignments:

                    ### Alignement format: '('UAAUUUCUACUCUUGUAGAUGAGAAGUCAUUUAA--UAAG-GC----CA-CUC', '-AAUUUCUACUCUUGUAGAUG-GAA---A-UU-AGGU--GCGCUUGGCAAC-C', 35.0, 0, 53)'
                    ### Alignement format:  Alignment(seqA='UAAUUUCUACUCUUGUAGAUGAGAAGUCAUUUAA--UAAG-GC----CA-CUC', seqB='-AAUUUCUACUCUUGUAGAUG-GAA---A-UU-AGGU--GCGCUUGGCAAC-C', score=35.0, start=0, end=53)
##                    score = str(alignment).split()
##                    score = score[2]
                    score = float(alignment.score)
##                    if score[:6] == 'score=': score = score.replace("score=", "")
##                    score = float(score.strip(','))

                                    
                    if(seq1_len >= seq2_len):
                        identity = (score/seq1_len) * 100
                    else:
                        identity = (score/seq2_len) * 100

                    f_SI_temp.write("%s\t%s\t%s\n" % (RP1, RP2, identity))
                    
                    if identity >= identity_threshold:
                        G.add_edge(RP1, RP2)
                        #print(identity)
                
                    break
                
    f_SI_temp.close()

    ######### Copy RP_SI_temp to RP_SI and delete RP_SI_temp #########
    shutil.copyfile('Organism_merge/RP_SI_temp','Organism_merge/RP_SI')
    os.remove("Organism_merge/RP_SI_temp")
    ######### Done #########
    



def connected_graph():
    
    f_grp_open = open("Organism_merge/Merge_based_on_identity_score", "w")
    connected_graph_list = list(nx.connected_components(G))

    for component in range(0, len(connected_graph_list)):
        f_grp_open.write("Group no: "+ str(component+1) + "\n")
        f_grp_open.write(str(connected_graph_list[component])+ "\n")

    f_grp_open.close()

            




def align_RP_seq(identity_threshold, output_with_org):

    print("\n ------- Removing organism based division ------- \n")

    print("\n ------- Generating sequence alignement among cluster representatives ------- \n")
    
    start_time = time.time()
    
    ### Create log file
    path = "Organism_merge/Log_file"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)
    ###
    
    read_NR_Table(output_with_org)
    read_RP_sequence()
    read_existing_SI()
    align_sequence(identity_threshold)
    connected_graph()
    
    run_time = round((time.time() - start_time)/3600,2)
    
    log_file = open('Organism_merge/Log_file/SI_runtime_emni', "w")
    log_file.write("--- %s hours ---" % (str(run_time)))
    log_file.close()







        

