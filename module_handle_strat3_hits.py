import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_target_prot_list(ID_transcriptome):
    dico_my_query = {}
    dico_target = {}
    new_path = "Intermediate_prot_BLAST"
    if os.path.exists(new_path) == False:
        os.system ("mkdir " + new_path)
    for seq_record in SeqIO.parse("DESMAN_denovo_output/denovo_protein.fa", "fasta"):
        name_pop = str(seq_record.id).split("_")[0]
        if name_pop == ID_transcriptome:
            dico_my_query[str(seq_record.id)] = str(seq_record.seq)
        else:
            dico_target[str(seq_record.id)] = str(seq_record.seq)
    query_prot_file = open("query_prot.fa", "w")
    for prot_name in dico_my_query.keys():
        query_prot_file.write(">" + prot_name + "\n")
        query_prot_file.write(dico_my_query[prot_name] + "\n")
    query_prot_file.close()
    target_prot_file = open("target_prot.fa", "w")
    for prot_name in dico_target.keys():
        target_prot_file.write(">" + prot_name + "\n")
        target_prot_file.write(dico_target[prot_name] + "\n")
    target_prot_file.close()

    os.system ("mv query_prot.fa Intermediate_prot_BLAST")
    os.system ("mv target_prot.fa Intermediate_prot_BLAST")
    

def handle_blast_inputs_transc():
    dico_denovo_prot_hits = {}
    file_blast = openFile("Intermediate_prot_BLAST/diamond_transc_prot.out")
    for line in file_blast:
        elts_line = line.split()
        query = elts_line[0]
        target = elts_line[1]
        pop_target = target.split("_")[0]
        if query not in dico_denovo_prot_hits.keys():
            dico_denovo_prot_hits[query] = {pop_target:[target]}
        else:
            if pop_target not in dico_denovo_prot_hits[query].keys():
                dico_denovo_prot_hits[query][pop_target] = [target]
            else:
                dico_denovo_prot_hits[query][pop_target].append(target)
    os.system ("rm Intermediate_prot_BLAST/diamond_transc_prot.out")
    os.system ("rm Intermediate_prot_BLAST/query_prot.fa")
    os.system ("rm Intermediate_prot_BLAST/target_prot.fa")
    return dico_denovo_prot_hits


def extract_coord(link_file):
    dico = {}
    info_file = openFile(link_file)
    for line in info_file[1:]:
        elts_line = line.split(",")
        chrom = elts_line[1]
        name = elts_line[9]
        start = int(elts_line[10])
        stop = int(elts_line[11])
        dico[name] = [chrom,start,stop]
    return dico


def assess_overlap(l1_old,l2_old):
    if l1_old[1] > l1_old[2]:
        l1 = [l1_old[0], l1_old[2], l1_old[1]]
    else:
        l1 = [l1_old[0], l1_old[1], l1_old[2]]
    if l2_old[1] > l2_old[2]:
        l2 = [l2_old[0], l2_old[2], l2_old[1]]
    else:
        l2 = [l2_old[0], l2_old[1], l2_old[2]]
    overlap = False
    if l1[0] == l2[0]:
        start1 = l1[1]
        start2 = l2[1]
        stop1 = l1[2]
        stop2 = l2[2]
        if start1 >= start2 and start1 <= stop2:
            overlap = True
        elif stop1 >= start2 and stop1 <= stop2:
            overlap = True
        elif start1 <= stop1 and stop1 >= stop2:
            overlap = True
        elif start1 >= start2 and stop1 <= stop2:
            overlap = True
        elif start1 == start2 and stop1 == stop2:
            overlap = True
    #if overlap == False:
        #print (l1)
        #print (l2)
        #print ("****")
    return overlap


def reduce_by_location(dico_homologs):
    correct_dico_homologs = {}
    dico_all_coord = extract_coord("DESMAN_denovo_output/information_file.txt")
    for denovo in dico_homologs.keys():
        subdico_homologs = dico_homologs[denovo]
        for pop in subdico_homologs.keys():
            list_homologs = subdico_homologs[pop]
            for target in list_homologs:
                overlap = assess_overlap(dico_all_coord[denovo], dico_all_coord[target])
                if overlap == True:
                    if denovo not in correct_dico_homologs.keys():
                        correct_dico_homologs[denovo] = [target]
                    else:
                        correct_dico_homologs[denovo].append(target)
    return correct_dico_homologs


def fill_dico_orthogroups(dico_my_hits, dico_orthogroups):
    nb = len(dico_orthogroups) + 1
    for denovo in dico_my_hits.keys():
        my_orthogroup = dico_my_hits[denovo]
        my_orthogroup.append(denovo)
        if nb == 1:
            newname = "orthogroup" + str(nb)
            dico_orthogroups[newname] = my_orthogroup
        else:
            already_present = False
            for orthogroup_name in dico_orthogroups.keys():
                orthogroup = dico_orthogroups[orthogroup_name]
                if denovo in orthogroup:
                    all_present = True
                    for name in my_orthogroup:
                        if name not in orthogroup:
                            all_present = False
                            break
                    if all_present == True:
                        already_present = True
                        break
            if already_present == False:
                newname = "orthogroup" + str(nb)
                dico_orthogroups[newname] = my_orthogroup
        nb += 1

def build_final_file_strat2(dico_orthogroups):
    os.system ("rm -r Intermediate_prot_BLAST")
    new_file = open("Orthogroup_output_step3.txt", "w")
    for orthogroup in dico_orthogroups.keys():
        new_file.write(orthogroup)
        for names in dico_orthogroups[orthogroup]:
            if names == dico_orthogroups[orthogroup][0]:
                new_file.write(":" + names)
            else:
                new_file.write("," + names)
        new_file.write("\n")
    new_file.close()
    os.system ("mv Orthogroup_output_step3.txt DESMAN_denovo_output")
