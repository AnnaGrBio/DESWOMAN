import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_transcripts_data_and_threeshold import *



def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def reshuffle_prot_file(link_prot1):
    # Opening the protein file for reading
    my_file1 = openFile(link_prot1)
    # Creating a new file name for the reshuffled protein file
    name_prot1_reshuffled = link_prot1.split(".")[0] + "_new.fa"
    
    # Opening the new file for writing
    my_new_file1 = open(name_prot1_reshuffled, "w")
    
    # Iterating through each line of the protein file
    for line in my_file1:
        my_new_file1.write(line) # new
        #if line[0] == ">":  #old # If the line starts with ">", it's a header line
            #elts_line = line.split()
            #print (elts_line)
            #my_gene = ">" + elts_line[1].split("=")[1] + "\n"  # Extracting gene name
            #my_new_file1.write(my_gene)  # Writing gene name to the new file
        #else:
            #if "." in line:  # If there is a dot in the line, likely a version number
                #new_line = line.split(".")[0] + "\n"  # Removing version number
                #my_new_file1.write(new_line)  # Writing modified line to the new file
            #else:
                #my_new_file1.write(line)  # Writing unchanged line to the new file
    
    # Closing the new file
    my_new_file1.close()
    
    # Returning the name of the reshuffled protein file
    return name_prot1_reshuffled


import os  # Importing the os module for system commands

def extract_all_prots_genome(pop_ind_name, path_to_genome):
    # Getting the current directory
    current_path = os.getcwd()

    # Constructing the path to the protein file to assess
    file_to_assess = current_path + "/old_prot_blast/" + pop_ind_name + "_prot_new.fa"

    # Checking if the file exists
    if os.path.isfile(file_to_assess) == False:
        print("File does not exist")
        
        # Constructing paths to required files
        link_output = path_to_genome + "/" + pop_ind_name + "_prot.fa"
        link_genome = get_path_to_fasta(path_to_genome, pop_ind_name)
        link_gff = get_path_gff(path_to_genome, pop_ind_name)
        
        # Constructing the command to extract proteins from GFF annotation
        command = "gffread -y " + link_output + " -g " + link_genome + " " + link_gff
        # Executing the command
        os.system(command)
        
        # Deleting temporary files
        file_to_delete = link_genome + ".fai"
        os.system("rm " + file_to_delete)
        print ("lili")
        
        # Reshuffling the protein file
        link_final_output = reshuffle_prot_file(link_output)
        print ("lala")
        os.system("rm " + link_output)
        
        # Creating a directory if it doesn't exist
        path_to_create = current_path + "/old_prot_blast"
        if os.path.isdir(path_to_create) == False:
            print("Creating path")
            os.system("mkdir old_prot_blast")
        
        # Moving the final output file to the appropriate directory
        os.system("mv " + link_final_output + " old_prot_blast")
    else:
        print("File exists")



def make_reciprocal_prot_blasts(query_name, target_name):
    # Getting the current directory
    current_path = os.getcwd()

    # Constructing paths to query and target protein files
    link_prot_query = current_path + "/old_prot_blast/" + query_name + "_prot_new.fa"
    link_prot_target = current_path + "/old_prot_blast/" + target_name + "_prot_new.fa"

    # Creating Diamond databases for query and target proteins
    command1 = "./diamond makedb --in " + link_prot_query + " -d " + query_name
    command2 = "./diamond makedb --in " + link_prot_target + " -d " + target_name
    os.system(command1)
    os.system(command2)

    # Performing reciprocal protein BLAST searches
    command_blast1 = "./diamond blastp -d " + target_name + " -q " + link_prot_query + " -o " + "old_prot_blast/prot_query_blast_out.txt"
    command_blast2 = "./diamond blastp -d " + query_name + " -q " + link_prot_target + " -o " + "old_prot_blast/prot_target_blast_out.txt"
    os.system(command_blast1)
    os.system(command_blast2)

    # Removing temporary Diamond database files
    os.system("rm *.dmnd")



def make_dico_blast():
    # Creating an empty dictionary for the first protein BLAST results
    dico_blast1 = {}
    # Opening the first protein BLAST output file for reading
    file_blast1 = openFile("old_prot_blast/prot_query_blast_out.txt")
    # Iterating through each line in the first BLAST output file
    for line in file_blast1:
        if line[0] != "#":  # Skipping comment lines
            elts_line = line.split()
            gene_query = elts_line[0]  # Extracting query gene
            gene_target = elts_line[1]  # Extracting target gene
            # Adding target gene to dictionary under query gene key
            if gene_query not in dico_blast1.keys():
                dico_blast1[gene_query] = [gene_target]
            else:
                dico_blast1[gene_query].append(gene_target)
    # Creating an empty dictionary for the second protein BLAST results
    dico_blast2 = {}
    # Opening the second protein BLAST output file for reading
    file_blast2 = openFile("old_prot_blast/prot_target_blast_out.txt")
    # Iterating through each line in the second BLAST output file
    for line in file_blast2:
        if line[0] != "#":  # Skipping comment lines
            elts_line = line.split()
            gene_query = elts_line[0]  # Extracting query gene
            gene_target = elts_line[1]  # Extracting target gene
            # Adding target gene to dictionary under query gene key
            if gene_query not in dico_blast2.keys():
                dico_blast2[gene_query] = [gene_target]
            else:
                dico_blast2[gene_query].append(gene_target)
    # Removing temporary BLAST output files
    os.system("rm old_prot_blast/prot_query_blast_out.txt")
    os.system("rm old_prot_blast/prot_target_blast_out.txt")
    # Returning both dictionaries
    return dico_blast1, dico_blast2
