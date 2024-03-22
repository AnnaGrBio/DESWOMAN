import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 

def search_surrounding_gene(start,stop,chrom,dico_gff_focal_sorted,window):
    list_surrounding = []
    #print (chrom)
    #print (start)
    if chrom in dico_gff_focal_sorted.keys():
        matrix = dico_gff_focal_sorted[chrom]
        classe = False
        if len(matrix) >= window*2:
            for i in range(0,len(matrix)-1):
                if int(start) >= int(matrix[i][0]) and int(start) <= int(matrix[i+1][0]):
                    window_start = i+1-window
                    window_stop = i+1+window
                    if window_start < 0:
                        window_start = 0
                    if window_stop > len(matrix):
                        window_stop = len(matrix)
                    for j in range(window_start,window_stop):
                        list_surrounding.append(matrix[j])
                        #print ("************")
                        classe = True
                    break
            if classe == False:
                if int(start) < int(matrix[0][0]):
                    for i in range(window):
                        list_surrounding.append(matrix[i])
                elif int(start) > int(matrix[len(matrix)-1][0]):
                    for i in range(1,window+1):
                        list_surrounding.append(matrix[len(matrix)-i])
        else:
            #print ("lala")
            list_surrounding = matrix
    else:
        list_surrounding = []
    #print (list_surrounding)
    #print ("***************")
    return list_surrounding


def search_overlapping_gene(start,stop,chrom,dico_gff_focal_sorted):
    list_overlap = []
    if chrom in dico_gff_focal_sorted.keys():
        matrix = dico_gff_focal_sorted[chrom]
        if len(matrix) == 1:
            list_overlap.append(matrix[0])
        else:
            for sublist in matrix:
                start_target = int(sublist[0])
                stop_target = int(sublist[1])
                if int(start) <= start_target:
                    if int(stop) >= start_target:
                        list_overlap.append(sublist)
                        #break
                elif int(start) >= start_target and int(stop) <= stop_target:
                    #print (start)
                    #print (stop)
                    #print (start_target)
                    #print (stop_target)
                    #print ("***************")
                    list_overlap.append(sublist)
                    #break
                elif int(start) <= stop_target and int(stop) >= stop_target:
                    list_overlap.append(sublist)
                    #break
    else:
        list_overlap = []
    return list_overlap





def determine_neigboor_function(denovo_genomic_infos,dico_gff_focal_sorted,window):
    start = denovo_genomic_infos[2]
    stop = denovo_genomic_infos[3]
    chrom = denovo_genomic_infos[1]
    if denovo_genomic_infos[0] == "intergenic_transcripts":
        list_neighboor = search_surrounding_gene(start,stop,chrom,dico_gff_focal_sorted,window)
    else:
        list_neighboor = search_overlapping_gene(start,stop,chrom,dico_gff_focal_sorted)
    return list_neighboor


def get_synteny_score(listA,listB, dico_blastA, dico_blastB):
    new_listA = []
    for list_gene in listA:
        new_listA.append(list_gene[2])
    new_listB = []
    for list_gene in listB:
        #print (list_gene)
        new_listB.append(list_gene[2])
    nb_reciprocal_hit = 0
    nb_gene_checked = 0
    for geneA in new_listA:
        nb_gene_checked += 1
        if geneA in dico_blastA.keys():
            list_hitsA = dico_blastA[geneA]
            for hitA in list_hitsA:
                if hitA in new_listB and hitA in dico_blastB.keys():
                        if geneA in dico_blastB[hitA]:
                            nb_reciprocal_hit += 1
                            break

    score = 100*nb_reciprocal_hit/nb_gene_checked
    return score



def search_all_denovo_hits(dico_denovo_informations,dict_my_blast_selected_outputs,dico_gff_focal_sorted,dico_gff_outgroup_sorted, window, dico_blast1, dico_blast2):
    dico_denovo_best_hit = {}
    list_best_score = []
    if window == 0:
        for denovo in dico_denovo_informations.keys():
            if denovo in dict_my_blast_selected_outputs.keys():
                best_hit = dict_my_blast_selected_outputs[denovo][0]
                dico_denovo_best_hit[denovo] = best_hit
    else:
        for denovo in dico_denovo_informations.keys():
            denovo_genomic_infos = dico_denovo_informations[denovo]
            if denovo in dict_my_blast_selected_outputs.keys():
                list_denovo_hits = dict_my_blast_selected_outputs[denovo]
                list_neighboor_denovo = determine_neigboor_function(denovo_genomic_infos,dico_gff_focal_sorted,window)
                best_hit = []
                best_score = -1
                for hit in list_denovo_hits:
                    list_elts_hit = hit.split("-")
                    hit_genomic_infos = [denovo_genomic_infos[0],list_elts_hit[2],list_elts_hit[0],list_elts_hit[1]]
                    list_neighboor_hit = determine_neigboor_function(hit_genomic_infos,dico_gff_outgroup_sorted,window)
                    score_synteny = get_synteny_score(list_neighboor_denovo,list_neighboor_hit, dico_blast1, dico_blast2)
                    if score_synteny > best_score:
                        best_hit = hit
                        best_score = score_synteny
                if best_score > 0:
                    dico_denovo_best_hit[denovo] = best_hit
                    list_best_score.append(1)
                else:
                    list_best_score.append(0)
            else:
                list_denovo_hits = []
        avg = sum(list_best_score)/len(list_best_score)
        print (avg)
    print (len(dico_denovo_informations))
    print (len(dico_denovo_best_hit))
    return dico_denovo_best_hit
    