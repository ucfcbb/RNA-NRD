import sys
import os
import logging
import time
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import multiprocessing as mp
from prepare_loops import *
from pymol_helper import *
from superimposition_generator import *
# python 3 compatibility


pdbx_dir = 'input_files/'
pdb_fasta_mapping_dir = 'pdb_fasta/'
#pdb_chains = {'6ypu': ['2','4'], '6ytf': ['3']}
#pdb_chains = {'6ypu': ['2','4']}
#pdb_chains = {'6ypu': ['2'], '6ytf': ['3']}


def calculate_rmsd_value(pdb_chains):

    #pdb_chains = {'5xh6': ['B'], '5xh7': ['B']}
    pdb_list = pdb_chains.keys()
    pdb_count = len(pdb_list)


    ########################## Download PDBX and FASTA files ##########################
    try:
        get_pdbx_and_fasta_files(pdb_list)
    except:
        #print("File not found")
        return -11


    ########################## Aligned Positions of Fasta Sequences ##########################
    final_ref_seq_dic = {}

    #print(pdb_chains)

    for key in pdb_chains:
        pdb_id = key
        chains = pdb_chains[key]
        ref_seq_dic = get_fasta_sequence(pdb_id, chains)

        for ref in chains:
            ref_key = pdb_id + '_' + ref
            final_ref_seq_dic[ref_key] = ref_seq_dic[ref]

    keylist = []
    for key in final_ref_seq_dic:
        keylist.append(key) 


    pbd1, chain1 = str(keylist[0]).split('_')
    pbd2, chain2 = str(keylist[1]).split('_')

##    if((pdb1 == "6HRM" and chain1 == "3" ) or (pdb2 == "6HRM" and chain2 == "3")):
##        return -66


    fasta_seq = final_ref_seq_dic[keylist[0]]
    ref_seq = final_ref_seq_dic[keylist[1]]
    aln = pairwise2.align.globalms(fasta_seq, ref_seq, 5, -3, -10, -1)
    (aln_fasta, aln_ref, _, _, _) = aln[0]
    ref_seq_replaced = replace_unknown_letter_in_ref(aln_ref, aln_fasta)
    fasta_to_ref_mapping, ref_to_fasta_mapping = get_aln_mapping(aln_fasta, ref_seq_replaced)
    #print(fasta_to_ref_mapping)



    ########################## Map FASTA sequences to PDB chain sequences ##########################
    generate_pdbx_fasta_mapping_files(pdb_chains)



    ########################## Store mapped positions of fasta sequences to PDBX chain in dictionaries ##########################
    ## pdb chain 1
    map_file1 = open(pdb_fasta_mapping_dir  + keylist[0] + ".rmsx.nch", "r")
    fasta_pdb_map_dic1 = {}

    while(True):
        line = map_file1.readline()
        if(line == ''):
            break

        residue, seq_pos = line.strip().split("\t")
        icode = ' '

        if "'" in residue:
            if "." in residue:
                residue, icode = residue.split(".")
            residue = residue.split("'")
            chain_id = residue[1]
            chain_pos = residue[2]
        else:
            chain_id = ''
            chain_pos = ''
            if "." in residue:
                residue, icode = residue.split(".")

            for i in range(0,len(residue)):
                if(residue[i].isalpha()):
                    chain_id = chain_id + residue[i]
                    continue
                chain_pos = chain_pos + residue[i] 
            chain_pos = int(chain_pos)

        fasta_pdb_map_dic1[int(seq_pos)] = (chain_id, int(chain_pos), icode)
    
    #print(fasta_pdb_map_dic1)

    # pdb chain2
    map_file2 = open(pdb_fasta_mapping_dir  + keylist[1] + ".rmsx.nch", "r")
    fasta_pdb_map_dic2 = {}

    while(True):
        line = map_file2.readline()
        if(line == ''):
            break

        residue, seq_pos = line.strip().split("\t")
        icode = ' '

        if "'" in residue:
            if "." in residue:
                residue, icode = residue.split(".")
            residue = residue.split("'")
            chain_id = residue[1]
            chain_pos = residue[2]
        else:
            chain_id = ''
            chain_pos = ''
            if "." in residue:
                residue, icode = residue.split(".")

            for i in range(0,len(residue)):
                if(residue[i].isalpha()):
                    chain_id = chain_id + residue[i]
                    continue
                chain_pos = chain_pos + residue[i] 
            chain_pos = int(chain_pos)

        fasta_pdb_map_dic2[int(seq_pos)] = (chain_id, int(chain_pos), icode)

    #print(fasta_pdb_map_dic2)


    ########################## Create a dictionary of aligned postions of the pdbx RNA chains ##########################
    # aligned_pdb_chain_pos_map = {}
    # for key in fasta_pdb_map_dic1:
    #     res_pos1 = fasta_pdb_map_dic1[key]

    #     if key in fasta_to_ref_mapping:
    #         key2 = fasta_to_ref_mapping[key]
    #         print(key2)
    #         res_pos2 = fasta_pdb_map_dic2[key]
    #         aligned_pdb_chain_pos_map[res_pos1] = res_pos2


    aligned_pdb_chain_pos_map = {}
    for key in fasta_to_ref_mapping:
        key2 = fasta_to_ref_mapping[key]

        if key in fasta_pdb_map_dic1 and key2 in fasta_pdb_map_dic2:
            res_pos1 = fasta_pdb_map_dic1[key]
            res_pos2 = fasta_pdb_map_dic2[key2]
            aligned_pdb_chain_pos_map[res_pos1] = res_pos2

    #print(aligned_pdb_chain_pos_map)

  
    ########################## Fetching residue coordinates from cif files od chain1 and chain2 ##########################
    chains = [chain1]
    chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict = get_pdbx_and_mapping_data(pbd1, chains)


    pdb_index_list = []
    for r in residue_dict[chain1]:
        (a, b, c) = (r.index.chain_id, r.index.seqnum, r.index.icode)
        pdb_index_list.append((a, b, c))
           # print(a, b, c)

    try:
        coord_backbone_chain1, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(pdbx_dir, pbd1 + '.cif'), pdb_index_list)
    except:
        #print("Encountered That error for PDBD id: " + pdb1 + " " + chain1)
        return -1
    #print(coord_backbone_chain1)

    chains = [chain2]
    chains, residue_dict, ref_seq_dict, res_to_ref, ref_to_res, missing_residue_dict = get_pdbx_and_mapping_data(pbd2, chains)

    pdb_index_list = []
    for r in residue_dict[chain2]:
        (a, b, c) = (r.index.chain_id, r.index.seqnum, r.index.icode)
        pdb_index_list.append((a, b, c))
           # print(a, b, c)
    try:
        coord_backbone_chain2, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(pdbx_dir, pbd2 + '.cif'), pdb_index_list)
    except:
        #print("Encountered That error for PDBD id: " + pdb2 + " " + chain2)
        return -2
    #print(coord_backbone_chain2)


    ########################## Keeeping only coordinates of aligned residues for both chain1 and chain2 ##########################
    aln_chain1_coord = {}
##    print(coord_backbone_chain2)
    for (a, b, c, d) in coord_backbone_chain1:
        if coord_backbone_chain1[(a, b, c, d)] != 0.0 and (b, c, d) in aligned_pdb_chain_pos_map.keys():
            (b2, c2, d2) = aligned_pdb_chain_pos_map[(b, c, d)]
            a2 = pbd2
            #b2 = chain2

            ### Special case for 6ZSC [6ZSC XA 2361 A]
            if a2 == '6ZSC':
                d2 = ' '
                if c2 == 5048 or c2 == 5055:
                    continue
            if a2 == '6ZSC' and b2 == 'XA' and c2 == 2360:
                continue
            
            ### Special case for 6ZSC [6ZSD XA 2359 B]
            if a2 == '6ZSD' and b2 == 'XA':
                d2 = ' '
                
            ### Special case for 6ZSC ['6ZSA', 'XA', 2359, 'B']
            if a2 == '6ZSA' and b2 == 'XA' and c2 == 2359:
                d2 = ' '

            ### Special case for 6ZSC ['6ZSG', 'XA', 2359, 'B']
            if a2 == '6ZSG' and b2 == 'XA' and c2 == 2359:
                d2 = ' '

            ### Special case for 6ZSC ['6ZSB', 'XA', 2359, 'B']
            if a2 == '6ZSB' and b2 == 'XA' and c2 == 2359:
                d2 = ' '
            
                
            if(coord_backbone_chain2[(a2, b2, c2, d2)] != 0.0):
                aln_chain1_coord[(a, b, c, d)] = coord_backbone_chain1[(a, b, c, d)]
    #print(len(aln_chain1_coord))
    #print(aln_chain1_coord)
           
            
    aln_chain2_coord = {}
    for (a, b, c, d) in coord_backbone_chain2:
        if coord_backbone_chain2[(a, b, c, d)] != 0.0 and (b, c, d) in aligned_pdb_chain_pos_map.values():
            key_list = list(aligned_pdb_chain_pos_map.keys())
            val_list = list(aligned_pdb_chain_pos_map.values())
            position = val_list.index((b, c, d))
            (b2, c2, d2) = key_list[position]
            a2 = pbd1
            #b2 = chain1


            ### Special case for 6ZSC [6ZSC XA 2361 A]
            if a2 == '6ZSC' and b2 == 'XA':
                d2 = ' '
                if c2 == 2360 or c2 == 5055 or c2 == 5048:
                    continue
            ### Special case for 6ZSC ['6ZSD', 'XA', 2359, 'A']
            if a2 == '6ZSD' and b2 == 'XA':
                d2 = ' '

                
            if(coord_backbone_chain1[(a2, b2, c2, d2)] != 0.0):
                aln_chain2_coord[(a, b, c, d)] = coord_backbone_chain2[(a, b, c, d)]
    #print(len(aln_chain2_coord))
    #print(aln_chain2_coord)



    ########################## RMSD Calculation ##########################
    coord1 = []
    coord2 = []

    for key in  aln_chain1_coord:
        coord1.append(aln_chain1_coord[key])

    for key in  aln_chain2_coord:
        coord2.append(aln_chain2_coord[key])

    coord1_list = list(filter(None, coord1))
    coord2_list = list(filter(None, coord2))

    X, Y = convert_array(coord1_list, coord2_list)

    # print(X)
    # print(Y)

    if len(X) != len(Y):
        logger.warning('WARNING: Corresponding co-ordinates for alignments not found! rmsd = 20 assigned.')
        rmsd = 20.
    elif len(X) == 0:
        logger.warning('WARNING: Co-ordinates for alignments not found! rmsd = 20 assigned.')
        rmsd = 20.
    else:
        XC = sum(X)/len(X)
        YC = sum(Y)/len(Y)
        # calculating relative co-ordinate using mean as reference
        X -= XC
        Y -= YC
        # time_s = time.time()
        rmsd = kabsch_rmsd(X, Y)

    #print(rmsd)
    return rmsd, len(coord1)
    #rmsd_file.write("%s\t%s\t%s\t%s\t%f\n" % (pdb1, chain1, pdb2, chain2, rmsd))



if __name__=="__main__":

    node1 = sys.argv[1]
    node2 = sys.argv[2]

    pdb1, chain1 = node1.split("_")
    pdb2, chain2 = node2.split("_")
    
    pdb_chains = {}
    if(pdb1 != pdb2):
        pdb_chains[pdb1] = [chain1]
        pdb_chains[pdb2] = [chain2]
    else:
        pdb_chains[pdb1] = [chain1, chain2]

        
##    print(pdb_chains)
    rmsd, align_len = calculate_rmsd_value(pdb_chains)
    print("RMSD: " , rmsd, " Align_len: ", align_len)

    #return rmsd
    
    #rmsd_file.write("%s\t%s\t%s\t%s\t%f\n" % (pdb1, chain1, pdb2, chain2, rmsd))
    #rmsd_file.write(pdb1 + "\t" + chain1 + "\t" + pdb2 + "\t" + chain2 + "\t" + str(rmsd) + "\n")





   


