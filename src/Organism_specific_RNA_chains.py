import os
import sys

def divide_chains_based_on_organism():

    print("\n ------- Dividing RNA chains based on organism ------- \n")

    os.chdir('src/')
##    path = os.getcwd()
##    print(path)
    
    All_organism_list = []
    All_organism_len = []
    Current_organism_list = []
    New_organism_list = []
    Organism_dic = {}
    Tot_chain_count = 0


    ### Read all organism list
    f_org = open("Organism_list/All_Organism_list", "r")

    while(True):
        line = f_org.readline()
        if line == "":
            break
        line = line.strip("\n")
        line = line.replace(" ", "_")
        All_organism_list.append(line)
        All_organism_len.append(len(line))
        org_file = line + '_f_out'
        Organism_dic[line] = [0, org_file]
        
    f_org.close()



    ### Read current organism list
    f_open = open("../RNA_chain_list_final", "r")

    line = f_open.readline()

    while(True):
        line = f_open.readline()
        
        if line == "":
            break

        line = line.strip("\n").split("\t")
        pdb_id, chain_id, entity_polymer_type, organism, chain_length, resolution, number_of_rna, release, experimental_method, macromolecule_name, family_name, num_coord, num_bp, sequence  = line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13]
        KNOWN_ORG_FLAG = 0
        
        if organism == "-": organism = "others"

        for i in range(len(All_organism_list)):
            
            org = All_organism_list[i]
            org_len = All_organism_len[i]

            if org == organism[:org_len]:

                file_name = Organism_dic[org][1]
                
                if org not in Current_organism_list:
                    
                    Current_organism_list.append(org)
                    file_name = open("Organism_chains/" + org + "_RNA_chains","w")
                    file_name.write("PDB_ID\tChain_ID\tEntity_Polymer_Type\tSource_Organism\tChain_Length\tResolution\t#_RNA_Entities\tRelease_Date\tExperimental_Method\tMacromolecule_Name\tFamily_name\tNum_res_with_coor\tNum_base_pairs\tSequence\n")
                    file_name.close()
                    
                file_name = open("Organism_chains/" + org + "_RNA_chains","a")
                file_name.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id, chain_id, entity_polymer_type, organism, chain_length, resolution, number_of_rna, release, experimental_method, macromolecule_name, family_name, num_coord, num_bp, sequence))
                file_name.close()
                
                KNOWN_ORG_FLAG = 1
                break

        if KNOWN_ORG_FLAG == 0:
            
            if organism not in Current_organism_list:
                
                print("Found new organism: " + organism)
                New_organism_list.append(organism)
                Current_organism_list.append(organism)
                
                f_org = open("Organism_list/All_Organism_list", "a")
                f_org.write(organism + "\n")
                f_org.close()
                
                file_name = open("Organism_chains/" + organism + "_RNA_chains","w")
                file_name.write("PDB_ID\tChain_ID\tEntity_Polymer_Type\tSource_Organism\tChain_Length\tResolution\t#_RNA_Entities\tRelease_Date\tExperimental_Method\tMacromolecule_Name\tFamily_name\tNum_res_with_coor\tNum_base_pairs\tSequence\n")
                file_name.close()

            file_name = open("Organism_chains/" + organism + "_RNA_chains","a")
            file_name.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pdb_id, chain_id, entity_polymer_type, organism, chain_length, resolution, number_of_rna, release, experimental_method, macromolecule_name, family_name, num_coord, num_bp, sequence))
            file_name.close()
             
    f_open.close()


    ### Write current organism names
    Current_organism_list.sort()

    f_cur_org = open("Organism_list/Current_Organism_list", "w")

    for og in Current_organism_list:
        f_cur_org.write(og + "\n")
        
    f_cur_org.close()

    print("No of organisms currently considered: ", len(Current_organism_list))
    print("Organism list: ", ', '.join(Current_organism_list))



    ### Write new organism names in the file 'All_Organism_list'
    for new_org in New_organism_list:
        All_organism_list.append(new_org)

    All_organism_list.sort()

    f_up_org = open("Organism_list/All_Organism_list", "w")

    for og in All_organism_list:
        f_up_org.write(og + "\n")
        
    f_up_org.close()    


