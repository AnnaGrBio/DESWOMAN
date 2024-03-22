import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def openFile(NameFile):
    
    """
    open file in lecture mode 
    """
    
    F=open(NameFile, "r")
    L=F.readlines()
    return L  

# Same here (est ce que je specifie un perc id?)
#parameters_db3 = {"type" : "blastp", "e_value" : "0.05", "coverage" : None, "strand" : None}
parameters_db2 = {"type" : "nucl", "e_value" : "0.05", "coverage" : 80, "strand" : "plus"}
parameters_db1 = {"type" : "nucl", "e_value" : "0.05", "coverage" : 80, "strand" : None}


def get_file_size(my_file):
    dico_name_size = {}
    for seq_record in SeqIO.parse(my_file, "fasta"):
        dico_name_size[str(seq_record.id)] = len(str(seq_record.seq))
    return dico_name_size


def retrieve_name_hit_nucl_blast(parameters_db):
    list_correct_hits = []
    name_output = "DESMAN_denovo_output/denovo_blast_output_nucl.txt"
    my_file = openFile(name_output)
    if parameters_db["coverage"] == None:
        for line in my_file:
            elts_line = line.split()
            if elts_line[0] != "#":
                list_correct_hits.append(elts_line[0])
    else:
        my_coverage = int(parameters_db["coverage"])
        for line in my_file:
            elts_line = line.split()
            if elts_line[0] != "#" and int(elts_line[7]) >= my_coverage:
                list_correct_hits.append(elts_line[0])
    list_correct_hits = list(set(list_correct_hits))
    return list_correct_hits


def retrieve_name_hit_prot_blast():
    list_correct_hits = []
    name_output = "DESMAN_denovo_output/denovo_blast_output_prot.txt"
    my_file = openFile(name_output)
    for line in my_file:
        elts_line = line.split()
        if elts_line[0] not in list_correct_hits:
            list_correct_hits.append(elts_line[0])
    list_correct_hits = list(set(list_correct_hits))
    return list_correct_hits


def merge_hit_lists(list1, list2):
    list_final = []
    if len(list1) > 0:
        for name in list1:
            list_final.append(name)
    if len(list2) > 0:
        for name in list2:
            list_final.append(name)
    if len(list_final) > 0:
        list_final = list(set(list_final))
    return list_final



def reshufe_files_in_denovo(list_denovo_to_delete):
    # denovo nucl
    os.system ("cat DESMAN_denovo_output/denovo_nucl.fa | grep \">\" | wc -l")
    my_file = open("denovo_nucl.fa", "w")
    for seq_record in SeqIO.parse("DESMAN_denovo_output/denovo_nucl.fa", "fasta"):
        if str(seq_record.id) not in list_denovo_to_delete:
            my_file.write(">"+str(seq_record.id)+"\n")
            my_file.write(str(seq_record.seq)+"\n")
    my_file.close()
    os.system ("mv denovo_nucl.fa DESMAN_denovo_output")
    os.system ("cat DESMAN_denovo_output/denovo_nucl.fa | grep \">\" | wc -l")

    # denovo nucl lowered intron
    if os.path.isfile("DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa"):
        os.system ("cat DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa | grep \">\" | wc -l")
        my_file = open("denovo_unspliced_lowered_introns.fa", "w")
        for seq_record in SeqIO.parse("DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa", "fasta"):
            if str(seq_record.id) not in list_denovo_to_delete:
                my_file.write(">"+str(seq_record.id)+"\n")
                my_file.write(str(seq_record.seq)+"\n")
        my_file.close()
        os.system ("mv denovo_unspliced_lowered_introns.fa DESMAN_denovo_output")
        os.system ("cat DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa | grep \">\" | wc -l")

    # denovo prot
    os.system ("cat DESMAN_denovo_output/denovo_protein.fa | grep \">\" | wc -l")
    my_file = open("denovo_protein.fa", "w")
    for seq_record in SeqIO.parse("DESMAN_denovo_output/denovo_protein.fa", "fasta"):
        if str(seq_record.id) not in list_denovo_to_delete:
            my_file.write(">"+str(seq_record.id)+"\n")
            my_file.write(str(seq_record.seq)+"\n")
    my_file.close()
    os.system ("mv denovo_protein.fa DESMAN_denovo_output")
    os.system ("cat DESMAN_denovo_output/denovo_protein.fa | grep \">\" | wc -l")

    # info_file
    my_file = open("information_file.txt", "w")
    study_file = openFile("DESMAN_denovo_output/information_file.txt")
    for line in study_file:
        name_orf = line.split(",")[9]
        if name_orf not in list_denovo_to_delete:
            my_file.write(line)
        #else:
            #print (name_orf)
    my_file.close()
    os.system ("mv information_file.txt DESMAN_denovo_output")
    



