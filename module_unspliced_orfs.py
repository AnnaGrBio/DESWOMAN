import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_orfs import * 



def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


### Functions part 4

def build_dict_unspliced_seqs(opened_gtf_file):
    """
    This function parses a GTF file containing transcriptome assembly information.
    It creates a dictionary where each transcript name is associated with a list containing
    chromosome, direction of transcription, and a sublist with exon start and stop positions.

    Parameters:
    - opened_gtf_file: The opened GTF file for transcriptome assembly.

    Returns:
    A dictionary where keys are transcript names, and values are lists containing chromosome,
    direction of transcription, and a sublist with exon start and stop positions.
    """
    
    # Initialize an empty dictionary to store transcript information
    dict_gene_exon_pos = {}
    transcript_name = ""

    # Iterate through each line in the opened GTF file
    for line in opened_gtf_file:
        # Exclude comment lines
        if line[0] != "#":
            # Split the line into elements
            list_elts_line = line.split("	")

            # Check if the line represents a transcript entry
            if list_elts_line[2] == "transcript":
                # If a previous transcript was being processed, store its information in the dictionary
                if transcript_name != "":
                    transcript_info = [chrom, direction, list_exons]
                    dict_gene_exon_pos[transcript_name] = transcript_info

                # Reset variables for the new transcript
                list_exons = []
                chrom = list_elts_line[0]
                direction = list_elts_line[6]
                
                # Extract transcript name from the line
                subLigne = list_elts_line[8].split(";")
                subsubLigne = subLigne[1].split(" ")
                transcript_name = subsubLigne[2][1:len(subsubLigne[2])-1]

            else:
                # Extract start and end positions for exons
                start_exon = int(list_elts_line[3])
                end_exon = int(list_elts_line[4])
                subList_exon = [start_exon, end_exon]
                list_exons.append(subList_exon)

    # Store the information for the last transcript in the dictionary
    transcript_info = [chrom, direction, list_exons]
    dict_gene_exon_pos[transcript_name] = transcript_info

    return dict_gene_exon_pos


def build_dict_chrom(genome_fasta): 
    
    """
    This function takes the fasta file containing the genome and builds a dictionary with all chromosomes, which is returned
    """
    
    dict_chrom = {}
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        dict_chrom[str(seq_record.id)] = str(seq_record.seq)
    return dict_chrom


def reverse_sequence(transcript): 
    
    """
    This funtion takes a fasta nucleotide sequence as an input are reverse transcribe it. I assume biopython would also have done that... returns the reversed sequence
    """
    
    reverse = ""
    Dico = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N","a":"t", "t":"a", "g":"c", "c":"g", "n":"n"}
    L = list(transcript.strip())
    for i in reversed(L):
        new_nucl = Dico[i]
        reverse+=new_nucl
    return reverse


def build_new_spliced_transcript(list_transcript_elts, dict_chromosomes_fasta):
    """
    This function takes as input the properties of a transcript and the genome sequence.
    It retrieves the unspliced transcript, slightly extended in both directions.
    
    Parameters:
    - list_transcript_elts: A list containing information about the transcript
        0: Name of the chromosome of the transcript
        1: Direction of transcription
        2: List of exons with their start and stop in the genome.
    - dict_chromosomes_fasta: A dictionary containing the genome in fasta format, chromosome by chromosome.

    Returns:
    As an output, it retrieve the unsplcied transcript, a bit longer in both direction. It returns such a transcript and the position of the splicing events in the genome.
    """
    
    chrom_name = list_transcript_elts[0]
    direction_transcription = list_transcript_elts[1]
    list_exons_pos = list_transcript_elts[2]
    chrom_sequence = dict_chromosomes_fasta[chrom_name]
    
    c = 3  # Counter for position in the spliced transcript
    dict_pos_spliced_transcript = {1: list_exons_pos[0][0] - 3, 2: list_exons_pos[0][0] - 2}

    if direction_transcription == "+":
        new_spliced_transcript = ""
        
        # Extend in the 5' direction
        new_spliced_transcript += chrom_sequence[list_exons_pos[0][0] - 3 : list_exons_pos[0][0] - 1]

        # Add exonic sequences
        for sublist_exon in list_exons_pos:
            my_exon_seq = chrom_sequence[sublist_exon[0] - 1 : sublist_exon[1]]
            for j in range(sublist_exon[0] - 1, sublist_exon[1]):
                dict_pos_spliced_transcript[c] = j
                c += 1
            new_spliced_transcript += my_exon_seq

        # Extend in the 3' direction
        if len(chrom_sequence) > list_exons_pos[len(list_exons_pos) - 1][1]:
            new_spliced_transcript += chrom_sequence[list_exons_pos[len(list_exons_pos) - 1][1]]
        dict_pos_spliced_transcript[c] = list_exons_pos[len(list_exons_pos) - 1][1]

    elif direction_transcription == "-" or direction_transcription == ".":  # Decouverte donc, les "." sont des reverse
        new_spliced_transcript = ""
        
        # Extend in the 5' direction
        new_spliced_transcript += chrom_sequence[list_exons_pos[0][0] - 3 : list_exons_pos[0][0] - 1]

        # Add exonic sequences
        for sublist_exon in list_exons_pos:
            my_exon_seq = chrom_sequence[sublist_exon[0] - 1 : sublist_exon[1]]
            for j in range(sublist_exon[0] - 1, sublist_exon[1]):
                dict_pos_spliced_transcript[c] = j
                c += 1
            new_spliced_transcript += my_exon_seq

        # Extend in the 3' direction
        if len(chrom_sequence) > list_exons_pos[len(list_exons_pos) - 1][1]:
            new_spliced_transcript += chrom_sequence[list_exons_pos[len(list_exons_pos) - 1][1]]
        dict_pos_spliced_transcript[c] = list_exons_pos[len(list_exons_pos) - 1][1]

        # Reverse the transcript if direction is "-" or "."
        new_spliced_transcript = reverse_sequence(new_spliced_transcript)

    return new_spliced_transcript, dict_pos_spliced_transcript


def max_dict(my_dict): 
    
    """
    This function takes the max of a dictionary
    """
    
    max_detected =0
    for i in my_dict.keys():
        if i>max_detected:
            max_detected = i
    return max_detected


def rebuild_spliced_orfs(genome_file, dict_gene_structure, dict_orfs):
    """
    This function reconstructs spliced transcripts based on GTF information, extracts spliced ORFs based on the GTF information,
    ensures the reconstructed ORFs correspond perfectly to the reconstructed ones. If validated, it retrieves the stop codon of
    the spliced ORF in the genome and the coordinates of the unspliced ORF. It returns a dictionary with all ORFs that had
    a stop codon validated by certain criteria, and takes the name of correct ORFs as an item. The list contains the coordinates
    in the genome of the unspliced ORF start, unspliced ORF stop, and the stop codon.

    Parameters:
    - genome_file: The file containing the genome information.
    - dict_gene_structure: A dictionary containing information about gene structure.
    - dict_orfs: A dictionary containing information about ORFs.

    Returns:
    A dictionary with information about validated ORFs and their coordinates in the genome.
    """
    dict_orf_unspliced_start_end_stop = {}
    dict_chromosomes_fasta = build_dict_chrom(genome_file)
    dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_orfs)
    stop_codon_list = ["TAG", "TAA", "TGA"]
    compteur_problematic_seq = 0
    dic_nucl_translation = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}

    for transcript_name in dict_gene_structure.keys():
        if transcript_name in dict_transcrit_all_orfs.keys():  # remove transcripts that do not have an ORF
            list_transcript_elts = dict_gene_structure[transcript_name]
            reconstructed_spliced_transcript, dict_pos_spliced_transcript = build_new_spliced_transcript(list_transcript_elts, dict_chromosomes_fasta)
            sequence_to_be_taken = True

            for orf_name in dict_transcrit_all_orfs[transcript_name]:
                if orf_name in dict_orfs.keys():
                    list_elts_orf_name = orf_name.split("_")
                    start_orf = int(list_elts_orf_name[2])
                    end_orf = int(list_elts_orf_name[3])
                    spliced_orf = dict_orfs[orf_name]

                    if list_transcript_elts[1] == "+":
                        reconstructed_spliced_orf = reconstructed_spliced_transcript[start_orf + 1 : end_orf + 2]
                        chromosome_fasta = dict_chromosomes_fasta[list_transcript_elts[0]]
                        start_orf_in_genome = dict_pos_spliced_transcript[start_orf + 2]
                        end_orf_in_genome = dict_pos_spliced_transcript[end_orf + 2]
                        end_orf_in_genome_with_stop = dict_pos_spliced_transcript[end_orf + 2] + 3
                        orf_stop_codon = chromosome_fasta[end_orf_in_genome + 1 : end_orf_in_genome + 4].upper()
                        #if orf_name == "STRG.17426.1_1_15_215": ### added
                            #print (spliced_orf)### added
                            #print ("*********")### added
                            #print (reconstructed_spliced_orf)### added

                        if orf_stop_codon not in stop_codon_list:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False

                        if spliced_orf != reconstructed_spliced_orf:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False


                    elif list_transcript_elts[1] == "-" or list_transcript_elts[1] == ".":  # Decouverte donc, les "." sont des reverse
                        reconstructed_spliced_orf = reconstructed_spliced_transcript[start_orf:end_orf + 1]
                        chromosome_fasta = dict_chromosomes_fasta[list_transcript_elts[0]]
                        start_orf_in_genome = dict_pos_spliced_transcript[max_dict(dict_pos_spliced_transcript) - start_orf]
                        end_orf_in_genome = dict_pos_spliced_transcript[max_dict(dict_pos_spliced_transcript) - end_orf]
                        end_orf_in_genome_with_stop = end_orf_in_genome - 3
                        orf_stop_codon = reverse_sequence(chromosome_fasta[end_orf_in_genome_with_stop:end_orf_in_genome]).upper()

                        if spliced_orf != reconstructed_spliced_orf:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False

                        if orf_stop_codon not in stop_codon_list:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False

                    if sequence_to_be_taken is True:
                        dict_orf_unspliced_start_end_stop[orf_name] = [start_orf_in_genome, end_orf_in_genome_with_stop, orf_stop_codon]

    return dict_orf_unspliced_start_end_stop


def add_stop_codon_to_spliced_orfs(dict_pos_unspliced_ORFs_plus_stop, dict_all_ORFs): 
    
    """
    this function creates a dictionary of all orfs that passed the stop codon and identity test, and and the corresponding stop codon to their fasta sequence (returns a dictionary with orf name as key and orf seq with the stop as an item).
    """
    
    dict_all_orfs_with_stop = {}
    for orf_name in dict_pos_unspliced_ORFs_plus_stop.keys():
        orf_seq = dict_all_ORFs[orf_name]
        orf_seq_with_stop = orf_seq + dict_pos_unspliced_ORFs_plus_stop[orf_name][2]
        dict_all_orfs_with_stop[orf_name] = orf_seq_with_stop
        #if orf_name == "STRG.17426.1_1_15_215": ### add
            #print (orf_seq_with_stop) ### add

    return dict_all_orfs_with_stop


