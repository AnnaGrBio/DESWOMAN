import os
from Bio import SeqIO


def get_list_name_query_rep(path_rep):
    list_query_name = []
    elts_dir = os.listdir(path_rep)
    for names in elts_dir:
        header = names.split(".")[0]
        if header not in list_query_name:
            list_query_name.append(header)
    return list_query_name
    

def get_path_gtf(path_rep,query_name):
    searched_name = query_name + "." + "gtf"
    elts_dir = os.listdir(path_rep)
    new_path = path_rep + "/" + searched_name
    return new_path

def get_path_to_fasta(path_rep,query_name):
    searched_name = ""
    elts_dir = os.listdir(path_rep)
    for file in elts_dir:
        lst = file.split(".")
        if lst[0] == query_name:
            if lst[1] == "fna" or lst[1] == "fa" or lst[1] == "fasta":
                searched_name = file
                break
    new_path = path_rep + "/" + searched_name
    return new_path

def get_path_gff(path_rep,query_name):
    searched_name = query_name + "." + "gff"
    elts_dir = os.listdir(path_rep)
    new_path = path_rep + "/" + searched_name
    return new_path



def assessTE(fasta_file):
    presenceTE = False
    list_nucl = ["a", "t", "c", "g", "n"]
    for nucl in list_nucl:
        if nucl in fasta_file:
            presenceTE = True
            break
    return presenceTE

def associate_gene_and_transcripts(opened_gtf_file, min_tpm_asked,  dict_transcript_fasta_all, te_exclusion): 

    """
    This function uses the gtf file of transcriptome assembly. It creates a dictionary with genes as keys and a list with all their transcripts as items.
    It returns the dictionary.
    """

    list_transcripts_to_exclude = []
    if te_exclusion == "True":
        for transcript_name in dict_transcript_fasta_all.keys():
            presenceTE = assessTE(dict_transcript_fasta_all[transcript_name])
            if presenceTE == True:
                list_transcripts_to_exclude.append(transcript_name)
    # Initialize an empty dictionary to store gene-transcript associations
    dico_genes = {}
    # Iterate through each line in the opened gtf file
    for line in opened_gtf_file:
        # Check if the line is not a comment
        if line[0] != "#":
            # Split the line into elements
            elts_line = line.split()

            if len(elts_line) > 9:
                # Check if the element is a transcript
                if elts_line[2] == "transcript":   # Note: Handle other applications if needed
                    # Extract gene_id, transcript_id, and TPM from the line
                    gene_id = elts_line[9][1:len(elts_line[9])-2]  # Note: Ensure this works for all cases
                    transcript_id = elts_line[11][1:len(elts_line[11])-2]  # Note: Ensure this works for all cases
                    TPM = elts_line[17][1:len(elts_line[17])-2]  # Note: Ensure this works for all cases

                    # Check if the TPM value is greater than or equal to the specified minimum TPM
                    if float(TPM) >= min_tpm_asked and transcript_id not in list_transcripts_to_exclude:
                        # Check if the gene_id is not in the dictionary, then add it
                        if gene_id not in dico_genes.keys():
                            dico_genes[gene_id] = [transcript_id]
                        else:
                            # If the gene_id is already in the dictionary, append the transcript_id to its list
                            dico_genes[gene_id].append(transcript_id)

    # Return the dictionary of gene-transcript associations
    return dico_genes


def transcripts_properties(opened_gtf_file):
    """
    This function uses the gtf file of transcriptome assembly. It creates a dictionary with transcripts as keys and a list with all their properties as items.
    The first key and item lists are already inserted in the dictionary and correspond to the properties we are looking at. It returns the dictionary.
    """

    # Initialize a dictionary with the first key and item lists
    dico_transcripts = {"key_transcript_name": ["chromosome", "direction", "start", "stop", "FPKM", "TPM", "coverage", "nb_exons"]}

    # Initialize an empty dictionary for inter-exon information
    dico_inter_exons = {}

    # Iterate through each line in the opened gtf file
    for line in opened_gtf_file:
        # Check if the line is not a comment
        if line[0] != "#":
            # Split the line into elements
            elts_line = line.split()

            # Check if the element is a transcript
            if elts_line[2] == "transcript":   # Note: Handle other applications if needed
                # Extract properties from the line for the transcript
                chromosome = elts_line[0]
                direction = elts_line[6]
                start = elts_line[3]
                stop = elts_line[4]
                transcript_id = elts_line[11][1:len(elts_line[11])-2]  # Note: Ensure this works for all cases
                FPKM = elts_line[15][1:len(elts_line[15])-2]  # Note: Ensure this works for all cases
                TPM = elts_line[17][1:len(elts_line[17])-2]  # Note: Ensure this works for all cases
                coverage = elts_line[13][1:len(elts_line[13])-2]  # Note: Ensure this works for all cases

                # Add transcript information to the dictionary
                dico_transcripts[transcript_id] = [chromosome, direction, start, stop, FPKM, TPM, coverage]

            # Check if the element is an exon
            elif elts_line[2] == "exon":    # Note: Handle other applications if needed
                # Extract properties from the line for the exon
                transcript_id = elts_line[11][1:len(elts_line[11])-2]  # Note: Ensure this works for all cases
                exon_nb = int(elts_line[13][1:len(elts_line[13])-2])  # Note: Ensure this works for all cases

                # Add or update the number of exons for the transcript in the inter-exon dictionary
                if transcript_id not in dico_inter_exons.keys():
                    dico_inter_exons[transcript_id] = exon_nb
                else:
                    if exon_nb > dico_inter_exons[transcript_id]:
                        dico_inter_exons[transcript_id] = exon_nb

    # Merge the inter-exon information into the main dictionary
    for transcripts in dico_transcripts.keys():
        if transcripts in dico_inter_exons.keys():
            dico_transcripts[transcripts].append(dico_inter_exons[transcripts])

    # Return the dictionary of transcript properties
    return dico_transcripts


def assess_te(fasta_seq):
    presence = False
    list_elts = ["a", "t", "c", "g", "n"]
    for nucl in fasta_seq:
        if nucl in list_elts:
            presence = True
            break
    return presence

def store_transcriptome(opened_transcriptome_file): 
    
    """
    This function opens the fasta file of the transcriptome, and store each transcripts in a dictionary with the name as a key and the sequence as an item. It returns the dictionary
    """
    
    dico_transcripts = {}
    for seq_record in SeqIO.parse(opened_transcriptome_file, "fasta"):
        dico_transcripts[str(seq_record.id)] = str(seq_record.seq)
    return dico_transcripts

