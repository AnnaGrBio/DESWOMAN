import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def validate_dict_content_genome(dico_final):
    list_line_to_delete = []  # Initialize an empty list to store keys to delete
    for line_species in dico_final.keys():
        # Check if the number of items (values) for each species is not equal to 2
        if len(dico_final[line_species]) != 2:
            # If it's not equal to 2, add the species key to the list for deletion
            list_line_to_delete.append(line_species)
    
    # Check if there are any species to delete
    if len(list_line_to_delete) > 0:
        # If there are species to delete, iterate through the list
        for i in list_line_to_delete:
            # Delete the species key from the dictionary
            del dico_final[i]


def validate_dict_content_transcriptome(dico_final):
    list_line_to_delete = []  # Initialize an empty list to store keys to delete
    for line_species in dico_final.keys():
        # Check if the number of items (values) for each species is equal to 3
        if len(dico_final[line_species]) == 3:
            # If it's equal to 3, add the species key to the list for further inspection
            list_line_to_delete.append(line_species)
    
    # Check if there are any species to delete
    if len(list_line_to_delete) > 0:
        # If there are species to delete, iterate through the list
        for i in list_line_to_delete:
            # Check if "transcriptome_gtf" key exists in the species dictionary
            if "transcriptome_gtf" in dico_final[i].keys():
                # If "transcriptome_gtf" key exists, delete it
                del dico_final[i]["transcriptome_gtf"]
            
            # Check if "transcriptome_fasta" key exists in the species dictionary
            if "transcriptome_fasta" in dico_final[i].keys():
                # If "transcriptome_fasta" key exists, delete it
                del dico_final[i]["transcriptome_fasta"]
    

def create_outgroup_dict_strategy1(dico_variables):
    dico_final = {}  # Initialize an empty dictionary to store final data
    
    # Extracting variables from the input dictionary
    path_to_genome = dico_variables["path_to_genome_repository"]
    path_to_transcriptome = dico_variables["path_to_transcriptome_repository"]
    name_query = dico_variables["query"]
    
    # Processing genome files
    elts_dir = os.listdir(path_to_genome)
    for file in elts_dir:
        lst = file.split(".")
        target_name = lst[0]
        
        # Checking if the file is not the query and has a valid extension
        if target_name != name_query:
            if lst[1] in ["fna", "fa", "fasta"]:  # Valid genome file extensions
                path_to_query_fasta_genome = path_to_genome + "/" + file
                if target_name not in dico_final.keys():
                    dico_final[target_name] = {"genome_fasta": path_to_query_fasta_genome}
                else:
                    dico_final[target_name]["genome_fasta"] = path_to_query_fasta_genome
            elif lst[1] in ["gff", "gff3"]:  # Valid GFF file extensions
                path_to_genome_gff = path_to_genome + "/" + file
                if target_name not in dico_final.keys():
                    dico_final[target_name] = {"genome_gff": path_to_genome_gff}
                else:
                    dico_final[target_name]["genome_gff"] = path_to_genome_gff
    
    # Validating genome data in the dictionary
    validate_dict_content_genome(dico_final)
    
    # Processing transcriptome files
    elts_dir2 = os.listdir(path_to_transcriptome)
    for file in elts_dir2:
        lst = file.split(".")
        target_name = lst[0]
        
        # Checking if the transcriptome file corresponds to a species in dico_final
        if target_name in dico_final.keys():
            if lst[1] == "gtf":  # Valid GTF file extension
                path_to_transcriptome_gtf = path_to_transcriptome + "/" + file
                dico_final[target_name]["transcriptome_gtf"] = path_to_transcriptome_gtf
            elif lst[1] in ["fna", "fa", "fasta"]:  # Valid transcriptome file extensions
                path_to_transcriptome_fasta = path_to_transcriptome + "/" + file
                dico_final[target_name]["transcriptome_fasta"] = path_to_transcriptome_fasta
    
    # Validating transcriptome data in the dictionary
    validate_dict_content_transcriptome(dico_final)
    
    return dico_final
