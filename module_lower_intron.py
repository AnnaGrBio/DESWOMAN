import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_orfs import * 
from module_unspliced_orfs import *



def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


### Functions part 5


def extract_intron_exon_pos(my_list_of_exon_pos):
    """
    This function generates a dictionary that associates each genomic position with an indication of whether it is
    part of an exon (E) or an intron (I) based on the provided list of exon positions.

    Parameters:
    - my_list_of_exon_pos (list): A list of exon positions, where each exon position is represented as a sublist
      with start and stop positions.

    Returns:
    dict: A dictionary where each genomic position is used as a key, and the corresponding value is "E" if the
      position is part of an exon or "I" if it is part of an intron.
    """
    dict_positions = {}

    # Iterate through each exon in the list
    for i in range(len(my_list_of_exon_pos)):
        list_exon_pos = my_list_of_exon_pos[i]

        # Mark each position within the exon as "E"
        for j in range(list_exon_pos[0] - 2, list_exon_pos[1] + 1):
            dict_positions[j] = "E"

        # If it's not the last exon, mark the positions between exons as "I"
        if i != len(my_list_of_exon_pos) - 1:
            for k in range(list_exon_pos[1] - 1, my_list_of_exon_pos[i + 1][0] ):
                dict_positions[k] = "I"

    return dict_positions 


def build_dico_intron_exon_specified(dict_gene_exon_pos, transcript_name):
    #print (transcript_name)
    """
    This function constructs dictionaries to store information about intron and exon positions for each transcript
    based on the provided dictionary of gene exon positions.

    Parameters:
    - dict_gene_exon_pos (dict): A dictionary where each transcript name is associated with genomic information,
      including the chromosome, direction, and a list of exon start and stop positions.

    Returns:
    tuple: A tuple containing two dictionaries:
      - dict_intron_exon_pos (dict): A dictionary where each transcript name is associated with another dictionary,
        specifying whether each genomic position is part of an exon (E) or an intron (I).
      - dict_transcript_chrom (dict): A dictionary associating each transcript name with its corresponding chromosome.
    """
    dict_intron_exon_pos = {}  # Dictionary to store intron and exon positions for each transcript
    dict_transcript_chrom = {}  # Dictionary to store the chromosome for each transcript

    coord_all_exons = dict_gene_exon_pos[transcript_name][2]

    # Extract intron and exon positions and store in a dictionary
    filing_dict = extract_intron_exon_pos(coord_all_exons)
    dict_intron_exon_pos[transcript_name] = filing_dict

    # Store the chromosome associated with each transcript
    dict_transcript_chrom[transcript_name] = dict_gene_exon_pos[transcript_name][0]

    return dict_intron_exon_pos, dict_transcript_chrom


def build_unspliced_fasta_seq_lowered_intron(dict_pos_unspliced_ORFs_plus_stop, genome_file, dict_gene_exon_pos):
    """
    This function generates unspliced FASTA sequences with lowered introns based on provided genomic information,
    unspliced ORF coordinates, chromosome details, and the genome file.

    Parameters:
    - dict_intron_exon_pos (dict): A dictionary where each transcript name is associated with intron and exon positions.
    - dict_pos_unspliced_ORFs_plus_stop (dict): A dictionary associating each unspliced ORF name with its start, end, and stop codon positions.
    - dict_transcript_chrom (dict): A dictionary associating each transcript name with its corresponding chromosome.
    - genome_file (str): The path to the file containing genomic information in FASTA format.

    Returns:
    dict: A dictionary where each unspliced ORF name is associated with its corresponding unspliced FASTA sequence
          with lowered introns.
    """
    dict_cord_unspliced_orfs = {}
    dict_unspliced_orfs_fasta_lowered_introns = {}
    dict_chromosomes_fasta = build_dict_chrom(genome_file)
    # Build a dictionary associating each unspliced ORF with its transcript, chromosome, start, and end positions
    for orf_name in dict_pos_unspliced_ORFs_plus_stop.keys():
        dict_intron_exon_pos, dict_transcript_chrom = build_dico_intron_exon_specified(dict_gene_exon_pos, orf_name.split("_")[0])
        #print (dict_intron_exon_pos.keys())
        dict_cord_unspliced_orfs[orf_name] = [orf_name.split("_")[0], dict_transcript_chrom[orf_name.split("_")[0]],
                                              dict_pos_unspliced_ORFs_plus_stop[orf_name][0],
                                              dict_pos_unspliced_ORFs_plus_stop[orf_name][1]]

        transcript_id = dict_cord_unspliced_orfs[orf_name][0]
        chromosome = dict_cord_unspliced_orfs[orf_name][1]
        start_unspliced_orf = dict_cord_unspliced_orfs[orf_name][2]
        end_unspliced_orf = dict_cord_unspliced_orfs[orf_name][3]
        my_transcript_splicing = dict_intron_exon_pos[transcript_id]

        if end_unspliced_orf > start_unspliced_orf:
            my_unspliced_seq_with_lowered_intron = ""

            # Build unspliced sequence with lowered introns
            for coord in range(start_unspliced_orf, end_unspliced_orf + 1):
                my_nucleotide = dict_chromosomes_fasta[chromosome][coord]
                if coord in my_transcript_splicing.keys():
                    if my_transcript_splicing[coord] == "E":
                        my_nucleotide_final = my_nucleotide.upper()
                    else:
                        my_nucleotide_final = my_nucleotide.lower()
                else:
                    my_nucleotide_final = my_nucleotide.upper()
                my_unspliced_seq_with_lowered_intron += my_nucleotide_final
            #if orf == "STRG.17426.1_1_15_215": ###### add:
                #print ("forward") ###### add:
                #print (dict_chromosomes_fasta[chromosome][start_unspliced_orf:end_unspliced_orf + 1]) ###### add:
                #print ("*****************") ###### add:
                #print (my_unspliced_seq_with_lowered_intron) ###### add:
            dict_unspliced_orfs_fasta_lowered_introns[orf_name] = my_unspliced_seq_with_lowered_intron
        else:
            my_unspliced_seq_with_lowered_intron = ""

            # Build reversed unspliced sequence with lowered introns
            for coord in range(end_unspliced_orf, start_unspliced_orf + 1):
                my_nucleotide = dict_chromosomes_fasta[chromosome][coord]
                if coord in my_transcript_splicing.keys():
                    if my_transcript_splicing[coord] == "E":
                        my_nucleotide_final = my_nucleotide.upper()
                    else:
                        my_nucleotide_final = my_nucleotide.lower()
                else:
                    my_nucleotide_final = my_nucleotide.upper()
                my_unspliced_seq_with_lowered_intron += my_nucleotide_final

            Final_unspliced_seq_with_lowered_intron = reverse_sequence(my_unspliced_seq_with_lowered_intron)
            #if orf == "STRG.17426.1_1_15_215": ###### add:
                #print ("reverse") ###### add:
                #lala = "" ###### add:
                #for coord in range(end_unspliced_orf, start_unspliced_orf + 1): ###### add:
                    #lala += dict_chromosomes_fasta[chromosome][coord] ###### add:
                #lala = reverse_sequence(lala) ###### add:
                #print (lala) ###### add:
                #print ("*****************") ###### add:
                #print (Final_unspliced_seq_with_lowered_intron) ###### add:
            dict_unspliced_orfs_fasta_lowered_introns[orf_name] = Final_unspliced_seq_with_lowered_intron

        #c += 1

    return dict_unspliced_orfs_fasta_lowered_introns
