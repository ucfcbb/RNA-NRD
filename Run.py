import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time

from Select_RNA_chains import *
from Organism_specific_RNA_chains import *
from Generate_File_for_Identity_Score import *
from Generate_Identity_Score import *
from Generate_Connected_Graph_from_identity_score import *
from Generate_File_for_RMSD import *
from Generate_RMSD_Preprocess import *
from Generate_RMSD_Align import *
from Cluster_Based_on_RMSD_Value import *
from RNA_NRD_with_org_div import *
from Merge_Organism_SI import *
from Merge_Organism_RMSD import *
from Merge_Organism_Members import *
from RNA_NRD_without_org_div import *



def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate Non-redundant Database')
    parser.add_argument('-i', nargs='?', default='Data.txt', const='Data.txt', help="Input file name containing RNA chains. Default: 'Data.txt'.")
    parser.add_argument('-o1', nargs='?', default='RNA-NRD_Dataset', const='RNA-NRD_Dataset', help="Output file name with organism division. Default: 'RNA-NRD_Dataset'.")
    parser.add_argument('-o2', nargs='?', default='RNA-NRD_without_Organism_Division_Dataset', const='RNA-NRD_without_Organism_Division_Dataset', help="Output file name without organism division. Default: 'RNA-NRD_without_Organism_Division_Dataset'.")
    parser.add_argument('-t', nargs='?', default=80, const=80, help='Provide sequence identity threshold. Default: 80.')
    parser.add_argument('-r', nargs='?', default=4, const=4, help='Provide RMSD threshold. Default: 4.')
    parser.add_argument('-a', nargs='?', default=80, const=80, help='Provide percent of structural alignment ratio threshold. Default: 80.')
    parser.add_argument('-org', nargs='?', default=True, const=True, help='Generate RNA-NRD dataset without organism based division. Default: True.')
    

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()

    user_input_fname = args.i
    output_with_org = args.o1
    output_without_org = args.o2
    identity_threshold = float(args.t)
    rmsd_threshold = float(args.r)
    align_ratio_threshold = float(args.a)
    organism_division = args.org
    


    select_RNA_chains(user_input_fname)
    divide_chains_based_on_organism()
    generate_file_for_identity_score_calculation()
    generate_identity_score()
    generate_connected_graph_from_identity_score(identity_threshold)
    generate_file_for_RMSD_calculation()
    generate_RMSD_preprocess()
    generate_RMSD_align()
    generate_cluster_based_on_RMSD(rmsd_threshold, align_ratio_threshold)
    generate_output_with_org(output_with_org)

    if organism_division == True:
        align_RP_seq(identity_threshold, output_with_org)
        align_RP_struct(rmsd_threshold, align_ratio_threshold, output_with_org)
        align_member_struct(rmsd_threshold, align_ratio_threshold, output_with_org)
        generate_output_without_org(output_without_org)

    


if __name__ == '__main__':

    os.system('chmod +x src/STAR3D_source_dssr/tools/RemovePseudoknots')
    os.system('chmod +x src/STAR3D_source_dssr/tools/MC-Annotate')
    
    main()
