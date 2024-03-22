import os
from Bio import SeqIO
from module_transcripts_data_and_threeshold import *


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def store_genome_informations(species_name, path_genome):
    link_to_gff = get_path_gff(path_genome,species_name)
    my_file = openFile(link_to_gff)
    dico_gff = {}
    for line in my_file:
        list_elts = line.split()
        if len(list_elts) > 2 and list_elts[2] == "mRNA":
            chromosome = list_elts[0]
            sublist = list_elts[8].split(";")[0]
            gene_ID = sublist.split("=")[1]
            start = list_elts[3]
            stop = list_elts[4]
            if chromosome not in dico_gff.keys():
                dico_gff[chromosome] = [[start,stop,gene_ID]]
            else:
                dico_gff[chromosome].append([start,stop,gene_ID])
    return (dico_gff)


def sort_gff_dic(dico_gff):
    new_dico_gff = {}
    for chromosome in dico_gff.keys():
        if len(dico_gff[chromosome]) > 1:
            list_ordered_start = []
            list_ordered_stops = []
            list_ordered_genes = []
            list_genes_pos = dico_gff[chromosome]
            for elt in list_genes_pos:
                gene = elt[2]
                gene_start = int(elt[0])
                gene_stop = int(elt[1])
                if len(list_ordered_start) == 0:
                    list_ordered_start.append(gene_start)
                    list_ordered_stops.append(gene_stop)
                    list_ordered_genes.append(gene)
                else:
                    place = False
                    for i in range(0,len(list_ordered_start)):
                        if list_ordered_start[i] > gene_start:
                            list_ordered_start.insert(i,gene_start)
                            list_ordered_stops.insert(i,gene_stop)
                            list_ordered_genes.insert(i,gene)
                            place == True
                            break
                    if place == False:
                        list_ordered_start.append(gene_start)
                        list_ordered_stops.append(gene_stop)
                        list_ordered_genes.append(gene)
            new_dico_gff[chromosome] = []
            for i in range(0,len(list_ordered_start)):
                sublist = [list_ordered_start[i],list_ordered_stops[i],list_ordered_genes[i]]
                new_dico_gff[chromosome].append(sublist)
        else:
            new_dico_gff[chromosome] = dico_gff[chromosome]
    #print (new_dico_gff["X_Chromosome"])
    return new_dico_gff
