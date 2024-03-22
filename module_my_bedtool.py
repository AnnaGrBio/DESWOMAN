import os
from module_sort_gff_gene_order import *


def openFile(NameFile):
    
    """
    open file in lecture mode 
    """
    
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def assess_gff_markers(my_gff):
    list_all_elts_2 = []
    for line in my_gff:
        if line[0] != "#":
            elts_2 = line.split()[2]
            if elts_2 not in list_all_elts_2:
                list_all_elts_2.append(elts_2)
    gene_marker = "gene"
    if "exon" in list_all_elts_2:
        exon_marker = "exon"
    else:
        exon_marker = "CDS"
    return gene_marker, exon_marker


def build_object_ID(elts):
    chrom_scaffold = elts[0]
    start = int(elts[3])
    stop = int(elts[4])
    orientation = elts[6]
    sublist = [chrom_scaffold, start, stop, orientation]
    return sublist


def extract_introns(list_exons):
    list_introns = []
    if list_exons[0][3] == "+":
        for my_index in range(len(list_exons) - 1):
            start_intron = list_exons[my_index][2] + 1
            stop_intron = list_exons[my_index + 1][1] - 1
            sublist_intron = [list_exons[my_index][0], start_intron, stop_intron, list_exons[my_index][3]]
            list_introns.append(sublist_intron)
    else:
        for my_index in range(len(list_exons) - 1):
            stop_intron = list_exons[my_index][1] - 1
            start_intron = list_exons[my_index + 1][2] + 1
            sublist_intron = [list_exons[my_index][0], start_intron, stop_intron, list_exons[my_index][3]]
            list_introns.append(sublist_intron)

    return list_introns


def dico_object_per_chrom(list_elts):
    dico = {}
    for sublist in list_elts:
        new_sublist = [sublist[1],sublist[2],sublist[3]]
        if sublist[0] not in dico.keys():
            dico[sublist[0]] = [new_sublist]
        else:
            dico[sublist[0]].append(new_sublist)
    return dico


def extract_gff_elts(link_to_gff):
    list_genes_data = []
    list_intron_data = []
    list_exon_data = []
    my_gff = openFile(link_to_gff)

    list_exons_my_gene = []
    gene_marker, exon_marker = assess_gff_markers(my_gff)

    for line in my_gff:
        if line[0] != "#":
            elts = line.split()
            if elts[2] == gene_marker:
                if len(list_exons_my_gene) > 1:
                    list_my_gene_introns = extract_introns(list_exons_my_gene)
                    for intron_list in list_my_gene_introns:
                        list_intron_data.append(intron_list)
                list_exons_my_gene = []
                #sublist_gene = build_object_ID(elts)
                #list_genes_data.append(sublist_gene)
            elif elts[2] == exon_marker:
                sublist_exon = build_object_ID(elts)
                list_exon_data.append(sublist_exon)
                list_exons_my_gene.append(sublist_exon)
    if len(list_exons_my_gene) > 1:
        list_my_gene_introns = extract_introns(list_exons_my_gene)
        for intron_list in list_my_gene_introns:
            list_intron_data.append(intron_list)
    dico_chrom_exon_pos = dico_object_per_chrom(list_exon_data)
    dico_chrom_intron_pos = dico_object_per_chrom(list_intron_data)
    return dico_chrom_exon_pos,dico_chrom_intron_pos


def my_bedtool_genes_introns_ordered(link_to_query_genome_gff):

    dico_chrom_exon_pos,dico_chrom_intron_pos = extract_gff_elts(link_to_query_genome_gff)
    dico_sorted_chrom_exon_pos = sort_gff_dic(dico_chrom_exon_pos)
    dico_sorted_chrom_intron_pos = sort_gff_dic(dico_chrom_intron_pos)

    return dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos


def find_intronic_overlap(dico_intron, dict_transcripts_properties, unknown_transcripts, antisense_transcript):
    list_intronic_transcripts = []
    for transcript in unknown_transcripts:
        if transcript in dict_transcripts_properties.keys():
            chrom_transcript = dict_transcripts_properties[transcript][0]
            orientation_transcript = dict_transcripts_properties[transcript][1]
            start_transcript = int(dict_transcripts_properties[transcript][2])
            stop_transcript = int(dict_transcripts_properties[transcript][3])
            if chrom_transcript in dico_intron.keys():
                #print (chrom_transcript)
                list_intron_chrom = dico_intron[chrom_transcript]
                intronic = False
                antisense = False
                for intron_spe_list in list_intron_chrom:
                    #print (intron_spe_list)
                    if intron_spe_list[0] <= start_transcript and intron_spe_list[1] >= stop_transcript:
                        if orientation_transcript == intron_spe_list[2]:
                            intronic = list_intronic_transcripts.append(transcript)
                        else:
                            antisense = antisense_transcript.append(transcript)
                        break
                #if intronic == True:
                    #list_intronic_transcripts.append(transcript)
                #else:
                    #if antisense == True:
                        #antisense_transcript.append(transcript)


    return list_intronic_transcripts, antisense_transcript


def find_genic_overlap(dico_exon, dict_transcripts_properties, list_accepted_transcript):
    list_genic = []
    list_antisense = []
    for transcript in list_accepted_transcript:
        if transcript in dict_transcripts_properties.keys():
            chrom_transcript = dict_transcripts_properties[transcript][0]
            orientation_transcript = dict_transcripts_properties[transcript][1]
            start_transcript = int(dict_transcripts_properties[transcript][2])
            stop_transcript = int(dict_transcripts_properties[transcript][3])
            if chrom_transcript in dico_exon.keys():
                list_exon_chrom = dico_exon[chrom_transcript]
                genic = False
                antisense = False
                for exon_spe_list in list_exon_chrom:
                    #print (intron_spe_list)
                    if start_transcript <= exon_spe_list[0] and stop_transcript > exon_spe_list[0]:
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                    if start_transcript < exon_spe_list[1] and stop_transcript >= exon_spe_list[1]:
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                    if start_transcript >= exon_spe_list[0] and stop_transcript <= exon_spe_list[1]:
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                if genic == True:
                    list_genic.append(transcript)
                else:
                    if antisense == True:
                        list_antisense.append(transcript)
    return list_genic, list_antisense


    


def my_bedtool_extract_overlap(dict_transcripts_properties,dict_gene_transcripts, dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos):
    print (len(dict_gene_transcripts))
    # Initialize a list to store all transcripts that will be used, based on if they are present in the dict_gene_transcripts, which contain only the transcripts validated previously with a specific TPM value.
    list_accepted_transcript = []

    # Iterate through genes and their transcripts to build the list of accepted transcripts
    for gene in dict_gene_transcripts.keys():
        for transcript in dict_gene_transcripts[gene]:
            list_accepted_transcript.append(transcript)

    # Initialize lists for different types of transcripts
    intergenic_transcripts = []
    genic_transcripts = [] # These one ovrlap with gens in the same direction and are not inside an intron
    antisense_transcript = []
    intronic_transcripts = []
    print("Total number of transcripts above the specified TPM thresshold : "+str(len(list_accepted_transcript)))
    genic_transcripts, antisense_transcript = find_genic_overlap(dico_sorted_chrom_exon_pos, dict_transcripts_properties, list_accepted_transcript)
    unknown_transcripts = []
    for transcript in list_accepted_transcript:
        if transcript not in genic_transcripts and transcript not in antisense_transcript:
            unknown_transcripts.append(transcript)
    intronic_transcripts, antisense_transcript = find_intronic_overlap(dico_sorted_chrom_intron_pos, dict_transcripts_properties, unknown_transcripts, antisense_transcript)
    intergenic_transcripts = []
    for transcript in unknown_transcripts:
        if transcript not in genic_transcripts and transcript not in antisense_transcript and transcript not in intronic_transcripts:
            intergenic_transcripts.append(transcript)
    # Print the counts of different types of transcripts
    print("Number of genic transcripts (not denovo) : " + str(len(genic_transcripts)))
    print("Number of intergenic transcripts : " + str(len(intergenic_transcripts)))
    print("Number of antisense transcripts : " + str(len(antisense_transcript)))
    print("Number of intronic transcripts (in gene frame) : " + str(len(intronic_transcripts)))
    return [intergenic_transcripts, intronic_transcripts, antisense_transcript, genic_transcripts]





def index_choice_transcript_overlap(list_choice_denovo):
    """
    This function maps user choices for transcript overlap types to corresponding indices.
    It takes a list of user choices and returns a list of indices based on the presence of specific choices.
    """

    # Initialize an empty list to store indices
    list_index = []

    # Check each type of transcript overlap and append the corresponding index if present in the user's choices
    if "intergenic" in list_choice_denovo:
        list_index.append(0)
    if "intronic" in list_choice_denovo:
        list_index.append(1)
    if "antisense" in list_choice_denovo:
        list_index.append(2)
    if "genic" in list_choice_denovo:
        list_index.append(3)

    # Return the list of indices
    return list_index



def extract_denovo_transcripts(sublist_transcripts_overlap, list_choice_denovo, dict_transcript_fasta, dict_gene_transcripts):
    """
    This function extracts denovo transcripts based on user-defined choices of transcript overlap types.
    It creates dictionaries for denovo gene transcripts, denovo transcript FASTA sequences, transcript status, and gene status.
    The function also generates a FASTA file containing denovo transcripts.
    """

    # Initialize dictionaries to store denovo gene transcripts, denovo transcript FASTA sequences, transcript status, and gene status
    dict_denovo_gene_transcripts = {}
    dict_denovo_transcript_fasta = {}
    dict_transcrip_status = {}
    dict_gene_status = {}

    # Open a FASTA file to write denovo transcripts
    fasta_denovo = open("denovo_transcript.fa", "w")

    # Map user choices to corresponding indices
    dict_index_value = {0: "intergenic", 1: "intronic", 2: "antisense", 3: "genic"}
    list_indexe_overlap = index_choice_transcript_overlap(list_choice_denovo)

    # Iterate through each chosen transcript overlap type
    for index in list_indexe_overlap:
        my_transcript_list = sublist_transcripts_overlap[index]

        # Iterate through transcripts in the chosen overlap type
        for transcript_name in my_transcript_list:
            gene_name = transcript_name.split(".")[0] + "." + transcript_name.split(".")[1]

            # Update dictionaries with denovo gene transcripts and denovo transcript FASTA sequences
            if gene_name in dict_denovo_gene_transcripts.keys():
                dict_denovo_gene_transcripts[gene_name].append(transcript_name)
            else:
                dict_denovo_gene_transcripts[gene_name] = [transcript_name]

            dict_denovo_transcript_fasta[transcript_name] = dict_transcript_fasta[transcript_name]
            dict_transcrip_status[transcript_name] = dict_index_value[index]

            # Write denovo transcripts to the FASTA file
            fasta_denovo.write(">" + transcript_name + "\n")
            fasta_denovo.write(dict_denovo_transcript_fasta[transcript_name] + "\n")

    # Close the FASTA file
    fasta_denovo.close()

    # Update gene status based on the number of denovo transcripts
    for gene in dict_denovo_gene_transcripts.keys():
        if len(dict_denovo_gene_transcripts[gene]) < len(dict_gene_transcripts[gene]):
            dict_gene_status[gene] = "genic"
        else:
            dict_gene_status[gene] = "denovo"

    # Return dictionaries containing denovo gene transcripts, denovo transcript FASTA sequences, transcript status, and gene status
    return dict_denovo_gene_transcripts, dict_denovo_transcript_fasta, dict_transcrip_status, dict_gene_status