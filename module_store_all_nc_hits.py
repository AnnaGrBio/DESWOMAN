import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord




def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def make_dico_len_unspliced_denovo():
    my_denovo_file = "DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa"
    dico_name_size = {}
    for seq_record in SeqIO.parse(my_denovo_file, "fasta"):
        dico_name_size[str(seq_record.id)] = len(str(seq_record.seq))
    return dico_name_size
    
    
def extract_all_blast_outputs(sp_pop_ind_name):
    dict_my_blast_selected_outputs = {}
    my_file_name = "blast_denovo_to_target_genome/" + sp_pop_ind_name + "_all_bast_output.txt"
    my_blast_results = openFile(my_file_name)
    for line in my_blast_results:
        if line[0] != "#":
            list_elts = line.split()
            denovo_name = list_elts[0]
            if int(list_elts[5]) < int(list_elts[6]):
                my_hit_name = list_elts[5]+"-"+list_elts[6]+"-"+list_elts[1] + "-f"
            else:
                my_hit_name = list_elts[6]+"-"+list_elts[5]+"-"+list_elts[1] + "-r"
            if denovo_name not in dict_my_blast_selected_outputs.keys():
                dict_my_blast_selected_outputs[denovo_name] = [my_hit_name]
            else:
                dict_my_blast_selected_outputs[denovo_name].append(my_hit_name)
    print (len(dict_my_blast_selected_outputs))
    return dict_my_blast_selected_outputs


def extract_best_blast_outputs(sp_pop_ind_name):
    dict_my_blast_selected_outputs = {}
    my_file_name = "blast_denovo_to_target_genome/" + sp_pop_ind_name + "_all_bast_output.txt"
    my_blast_results = openFile(my_file_name)
    for line in my_blast_results:
        if line[0] != "#":
            list_elts = line.split()
            denovo_name = list_elts[0]
            if int(list_elts[5]) < int(list_elts[6]):
                my_hit_name = list_elts[5]+"-"+list_elts[6]+"-"+list_elts[1] + "-f"
            else:
                my_hit_name = list_elts[6]+"-"+list_elts[5]+"-"+list_elts[1] + "-r"
            if denovo_name not in dict_my_blast_selected_outputs.keys():
                dict_my_blast_selected_outputs[denovo_name] = my_hit_name
    print (len(dict_my_blast_selected_outputs))
    return dict_my_blast_selected_outputs


def store_de_novo_informations():
    dico_denovo_info = {}
    my_file = openFile("DESMAN_denovo_output/information_file.txt")
    for line in my_file[1:]:
        line = line.split("\n")[0]
        name_denovo = line.split(",")[9]
        position = line.split(",")[12]
        chromosome = line.split(",")[1]
        start = line.split(",")[6]
        stop = line.split(",")[7]
        dico_denovo_info[name_denovo] = [position,chromosome,start,stop]
    return dico_denovo_info
