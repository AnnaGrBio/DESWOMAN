import os
from module_colors import *


def remove_intermediates(name_intermediate_directory):
    file_to_test = name_intermediate_directory + "/output_last_blast_denovo_to_homolog.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_last_blast_to_genome.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_last_diamond_gene_rec_blast_1.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)
    
    file_to_test = name_intermediate_directory + "/output_last_diamond_gene_rec_blast_2.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_last_diamond_gene_rec_blast_makedb_1.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_last_diamond_gene_rec_blast_makedb_2.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_step2_prot_diamond.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)

    file_to_test = name_intermediate_directory + "/output_step2_prot_makeblastdb.txt"
    if os.path.isfile(file_to_test):
        os.remove(file_to_test)