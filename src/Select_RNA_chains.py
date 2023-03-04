import subprocess
import os.path
import sys
import glob
from os import path
import wget
import requests


def select_RNA_chains(input_file_name):
    
    print(" ------- Selecting RNA chains ------- \n")                
    f_family = open("Rfam_family/PDB_family_name.txt", "r")
    family_list = {}
    fam_len = 0
    new = 0
    old = 0

    while(True):
        line = f_family.readline()
        if line == '':
            break
        if line[0:7] == 'Line_no':
            continue
        PDB, CHAIN, FAMILY = line.strip('\n').split('\t')
        key = PDB + '_' + CHAIN
        if key not in family_list:
            family_list[key] = []
            family_list[key].append(FAMILY)
            new += 1
        else:
            family_list[key].append(FAMILY)

            old += 1
            
        fam_len += 1

    f_family.close()


    f1_open = open("Data/" + input_file_name, "r")
    f_out = open("RNA_chain_list_final","w")
    f_out.write("PDB_ID\tChain_ID\tEntity_Polymer_Type\tSource_Organism\tChain_Length\tResolution\t#_RNA_Entities\tRelease_Date\tExperimental_Method\tMacromolecule_Name\tFamily_name\tNum_res_with_coor\tNum_base_pairs\tSequence\n")


    resolution = '-'
    number_of_rna = '-'
    release = '-'
    experimental_method = '-'
    pdb_id = '-'

    line = f1_open.readline()
    RNA_chain_count = 0

    while(True):
        
        line = f1_open.readline()
        if(line == ''):
            break

        line = line.split('\t')

        entry_id, cur_pdb_id, chain_id, source_organism, chain_length, cur_resolution, cur_number_of_rna, cur_release, cur_experimental_method, entity_polymer_type, macromolecule_name, sequence  = line[0], line[12].strip('\n'), line[11].strip('"'), line[9], line[8], line[5], line[4], line[2], line[1], line[7], line[10], line[6]

        
        if(entry_id == ''):
            entry_id = '-'
        if(pdb_id == ''):
            pdb_id = '-'
        if chain_id == '':
           chain_id = '-'
        if(source_organism == ''):
            source_organism = '-'
        if(chain_length == ''):
            chain_length = '-'
        if cur_pdb_id == '':
           cur_pdb_id = '-'
        if cur_resolution == '':
           cur_resolution = '-'
        if cur_number_of_rna == '':
           cur_number_of_rna = '-'
        if(cur_release == ''):
            cur_release = '-'
        if(cur_experimental_method == ''):
            cur_experimental_method = '-'
        if entity_polymer_type == '':
           entity_polymer_type = '-'     
        if macromolecule_name == '':
            macromolecule_name = '-'
        if sequence == '':
            sequence = '-'

        if(cur_pdb_id != '-'):
            pdb_id = cur_pdb_id

        if(entry_id == pdb_id):
            resolution = cur_resolution
            number_of_rna = cur_number_of_rna
            release = cur_release
            experimental_method = cur_experimental_method      

            
        if(entity_polymer_type == 'RNA' and int(chain_length) >= 20):
            chain_id = chain_id.replace('"', '').strip('\n')
            chain_id = chain_id.replace('_', '')
            chain_id_list = chain_id.split(',')
            
            for c in range(0,len(chain_id_list)):

                RNA_chain_count += 1
                family_name = 'other'
                pdb_key = pdb_id + '_' + str(chain_id_list[c])

                
                coord_num = '-'
                bp_num = 0
                
                if pdb_id != '-' and chain_id_list[c] != '-':

                    ### Adding number of residues having coordinates in a chain using modified code from RNAmotifContrast
                    #coord_calc_path = 'src/RNAMotifContrast/src/scripts/'
                    coord_calc_path = 'src/Identical_Align/src/scripts/'
                    p3 = subprocess.check_output('python3 -W ignore Structure_coordinate_length.py '+ pdb_key, shell = True, cwd = coord_calc_path)
                    p3 = p3.decode("utf-8")
                    p3 = p3.strip().split(" ")
                    coord_num = p3[2]

                    ### Adding number of base pairs in a chain from a chain using DSSR
                    DSSR_path = 'DSSR/'
                    pdb_id_low = pdb_id.lower()
                    annotation_path = DSSR_path + pdb_id_low + ".out"
                    link = 'http://skmatic.x3dna.org/pdb/' + pdb_id_low + '/' + pdb_id_low + '.out'

                   
                    if os.path.exists(annotation_path):
                        print(pdb_id_low + " : annotation exists")
                    else:
                        try:
                            print(pdb_id_low + " : downloading annotation")
                            #wget.download(link, DSSR_path, bar=None)
                            annofile = requests.get(link, allow_redirects=False)
                            #print(annofile.status_code)
                            if annofile.status_code == 200:
                                open(annotation_path, 'wb').write(annofile.content)
                            else:
                                print("Missing annotation for PDB file " + pdb_id_low)
                                print("Solution: Generate DSSR annotation for the PDB file using tool DSSR or from website 'http://skmatic.x3dna.org/'. Then put inside the folder 'RNA_NRD_source_code/DSSR'. Finally run the code again!")
                                sys.exit()                                
                        except Exception as e:
                            print("Missing annotation for PDB file " + pdb_id_low)
                            print("Solution: Generate DSSR annotation for the PDB file using tool DSSR or from website 'http://skmatic.x3dna.org/'. Then put inside the folder 'RNA_NRD_source_code/DSSR'. Finally run the code again!")
                            sys.exit()
                        
                    
                    if(os.path.isfile(DSSR_path + pdb_id_low + '.out')):
                        bp_open = open(DSSR_path + pdb_id_low + '.out', "r")

                        ### Skipping lines
                        while(True):
                            cur_line = bp_open.readline()

                            if cur_line == '':
                                break

                            cur_line = cur_line.split(" ")

                            ### Reading base-pairs
                            if cur_line[0] == "List" and cur_line[1] == "of" and cur_line[3] == "base" and cur_line[4] == "pairs":

                                line = bp_open.readline()
                                line = bp_open.readline()
                        
                                while(True):
                                    line = bp_open.readline()
                                    
                                    if(line[:3] == '   '):
                                        break
                                    
                                    line = line.split()
                                    nt1, nt2 = line[1].split('.')[0], line[2].split('.')[0]
                                    
                                    if nt1 == chain_id_list[c] and nt2 == chain_id_list[c]:
                                        bp_num += 1
                                
                                
                ### Adding R-fam Family Info
                if pdb_key in family_list:
                    family_name = family_list[pdb_key]                
                f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id, str(chain_id_list[c]), entity_polymer_type, source_organism, chain_length, resolution, number_of_rna, release, experimental_method, macromolecule_name, family_name, coord_num, str(bp_num), sequence))

                      
    f1_open.close()
    f_out.close()
