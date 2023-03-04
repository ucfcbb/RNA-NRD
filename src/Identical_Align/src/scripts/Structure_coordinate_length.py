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



def calculate_coordinate(pdb_chains):

    #pdb_chains = {'5xh6': ['B'], '5xh7': ['B']}
    pdb_list = pdb_chains.keys()
    pdb_count = len(pdb_list)
    keys = list(pdb_chains.keys())
    values = pdb_chains.values()
    pbd1 = keys[0]
    chain1 = pdb_chains[pbd1][0]


    ########################## Download PDBX and FASTA files ##########################
    try:
        get_pdbx_and_fasta_files(pdb_list)
    except:
        #print("File not found")
        return -11

    ########################## Map FASTA sequences to PDB chain sequences ##########################

    generate_pdbx_fasta_mapping_files(pdb_chains)


  
    ########################## Fetching residue coordinates from cif files of chain1 ##########################
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
        return -1


    ########################## Counting number of residues having coordinate ##########################
    coordinate_count = 0

    for key in coord_backbone_chain1:
        if coord_backbone_chain1[key] != 0.0:
            coordinate_count += 1
    return coordinate_count
    



if __name__=="__main__":

    node1 = sys.argv[1]
    pdb1, chain1 = node1.split("_")
    pdb_chains = {}
    pdb_chains[pdb1] = [chain1]
    

    coordinate_count = calculate_coordinate(pdb_chains)
    print("Coordinate_count: " , coordinate_count)






   


