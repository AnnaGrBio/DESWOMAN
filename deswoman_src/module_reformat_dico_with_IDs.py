import os
from Bio import SeqIO


__author__ = "Anna Grandchamp"
__contributor__=""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def implement_dict_transcripts_properties(dict_transcripts_properties : dict, dict_transcripts_properties_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the keys of a dictionary of transcript properties and stores them in a new dictionary.

    This function takes a dictionary of transcript properties (`dict_transcripts_properties`), adds a prefix (the 
    transcriptome ID) to each transcript name, and stores the result in a new dictionary (`dict_transcripts_properties_total`).
    The prefix helps to uniquely identify the transcripts from different transcriptomes.

    Parameters:
    -----------
    dict_transcripts_properties : dict
        A dictionary where keys are transcript names and values are lists or tuples containing the properties of each transcript. 

    dict_transcripts_properties_total : dict
        The dictionary where the new entries (with the prefixed transcript names) will be stored.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each transcript name in the 
        `dict_transcripts_properties` to generate unique keys in the `dict_transcripts_properties_total`.

    Returns:
    --------
    None
        The function modifies `dict_transcripts_properties_total` in place by adding the new transcript names and their properties.
    """
    for transcript_name in dict_transcripts_properties:
        new_transcript_name = ID_transcriptome+"_" + transcript_name
        dict_transcripts_properties_total[new_transcript_name] = dict_transcripts_properties[transcript_name]


def implement_dict_gene_status(dict_gene_status : dict, dict_gene_status_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the keys of a dictionary of gene statuses and stores them in a new dictionary.

    This function takes a dictionary of gene statuses (`dict_gene_status`), adds a prefix (the 
    transcriptome ID) to each gene name, and stores the result in a new dictionary (`dict_gene_status_total`).
    The prefix helps to uniquely identify the genes from different transcriptomes.

    Parameters:
    -----------
    dict_gene_status : dict
        A dictionary where keys are gene names and values are the status of each gene (e.g., "genic" or "non-genic").

    dict_gene_status_total : dict
        The dictionary where the new entries (with the prefixed gene names) will be stored.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each gene name in the 
        `dict_gene_status` to generate unique keys in the `dict_gene_status_total`.

    Returns:
    --------
    None
        The function modifies `dict_gene_status_total` in place by adding the new gene names and their statuses.
    """
    for gene_name in dict_gene_status:
        new_gene_name = ID_transcriptome + "_" + gene_name
        dict_gene_status_total[new_gene_name] = dict_gene_status[gene_name]


def implement_dict_gene_transcripts_denovo(dict_gene_transcripts_denovo : dict, dict_gene_transcripts_denovo_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the gene and transcript names in a dictionary of gene-transcript associations.

    This function takes a dictionary (`dict_gene_transcripts_denovo`) where keys are gene names and values are lists 
    of associated transcript names. It adds a prefix (the transcriptome ID) to both gene and transcript names, 
    then stores the result in a new dictionary (`dict_gene_transcripts_denovo_total`).

    Parameters:
    -----------
    dict_gene_transcripts_denovo : dict
        A dictionary where keys are gene names and values are lists of associated transcript names. These represent 
        the relationships between genes and their corresponding de novo transcripts.

    dict_gene_transcripts_denovo_total : dict
        The dictionary where the new entries (with the prefixed gene and transcript names) will be stored.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each gene and transcript name 
        in `dict_gene_transcripts_denovo` to generate unique keys in the `dict_gene_transcripts_denovo_total`.

    Returns:
    --------
    None
        The function modifies `dict_gene_transcripts_denovo_total` in place by adding the new gene-transcript relationships.
    """
    for gene_name in dict_gene_transcripts_denovo:
        new_gene_name = ID_transcriptome + "_" + gene_name
        new_list_transcs = []
        for transcript in dict_gene_transcripts_denovo[gene_name]:
            new_transcript_name = ID_transcriptome + "_" + transcript
            new_list_transcs.append(new_transcript_name)
        dict_gene_transcripts_denovo_total[new_gene_name] = new_list_transcs


def implement_dict_transcrit_filtered_orfs(dict_transcrit_filtered_orfs : dict, dict_transcrit_filtered_orfs_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the transcript and ORF names in a dictionary of transcript-ORF associations.

    This function takes a dictionary (`dict_transcrit_filtered_orfs`) where keys are transcript names and values 
    are lists of associated ORF names. It adds a prefix (the transcriptome ID) to both transcript and ORF names, 
    then stores the result in a new dictionary (`dict_transcrit_filtered_orfs_total`).

    Parameters:
    -----------
    dict_transcrit_filtered_orfs : dict
        A dictionary where keys are transcript names and values are lists of associated ORF names. These represent 
        the relationships between filtered transcripts and their associated ORFs.

    dict_transcrit_filtered_orfs_total : dict
        The dictionary where the new entries (with the prefixed transcript and ORF names) will be stored.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each transcript and ORF name 
        in `dict_transcrit_filtered_orfs` to generate unique keys in the `dict_transcrit_filtered_orfs_total`.

    Returns:
    --------
    None
        The function modifies `dict_transcrit_filtered_orfs_total` in place by adding the new transcript-ORF relationships.
    """
    for transcript_name in dict_transcrit_filtered_orfs:
        new_transcript_name = ID_transcriptome + "_" + transcript_name
        list_new_orfs_names = []
        for orf in dict_transcrit_filtered_orfs[transcript_name]:
            new_orf_name = ID_transcriptome + "_" + orf
            list_new_orfs_names.append(new_orf_name)
        dict_transcrit_filtered_orfs_total[new_transcript_name] = list_new_orfs_names


def implement_dict_pos_unspliced_ORFs_plus_stop(dict_pos_unspliced_ORFs_plus_stop : dict, dict_pos_unspliced_ORFs_plus_stop_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the ORF names and stores their positions in a new dictionary.

    This function iterates over a dictionary of unspliced ORFs with stop codons, adding a prefix (the transcriptome ID) 
    to each ORF name. It then stores the ORF names along with their corresponding positions in a new dictionary.

    Parameters:
    -----------
    dict_pos_unspliced_ORFs_plus_stop : dict
        A dictionary where keys are unspliced ORF names (without the transcriptome ID prefix) and values are lists or tuples
        representing the genomic positions of these ORFs (e.g., start and stop positions).
    
    dict_pos_unspliced_ORFs_plus_stop_total : dict
        The dictionary where the new ORF names (with the prefix) and their associated positions will be stored.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each ORF name to ensure unique
        naming within the new dictionary.

    Returns:
    --------
    None
        The function modifies `dict_pos_unspliced_ORFs_plus_stop_total` in place by adding the new ORF names with their positions.
    """
    for orf_name in dict_pos_unspliced_ORFs_plus_stop:
        new_orf_name = ID_transcriptome + "_" + orf_name
        dict_pos_unspliced_ORFs_plus_stop_total[new_orf_name] = dict_pos_unspliced_ORFs_plus_stop[orf_name]


def implement_sublist_transcripts_overlap(sublist_transcripts_overlap : list, sublist_transcripts_overlap_total : list, ID_transcriptome : str) -> list:
    """
    Adds a prefix (name of the considered transcriptome) to the transcript names in the overlap sublists and updates the total list.

    This function takes a list of sublists (representing different categories of overlapping transcripts) and adds a prefix 
    (the transcriptome ID) to each transcript name. It then appends the updated transcript names to a corresponding 
    sublist in a total list, ensuring that the overlap categories are maintained.

    Parameters:
    -----------
    sublist_transcripts_overlap : list of lists
        A list of 4 sublists, where each sublist contains transcript names that belong to specific overlap categories:
        - Index 0: Intergenic Transcripts
        - Index 1: Intronic Transcripts
        - Index 2: Antisense Transcripts
        - Index 3: Genic Transcripts
    
    sublist_transcripts_overlap_total : list of lists
        The total list of sublists where the updated transcript names (with the prefix) will be stored. Each sublist in 
        this list corresponds to a category from `sublist_transcripts_overlap`.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each transcript name.

    Returns:
    --------
    sublist_transcripts_overlap_total : list of lists
        The updated list of sublists with the transcript names prefixed with the transcriptome ID.
    """
    if len(sublist_transcripts_overlap_total) == 0:
        sublist_transcripts_overlap_total = [[], [], [], []]
    for transcript in sublist_transcripts_overlap[0]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[0].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[1]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[1].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[2]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[2].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[3]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[3].append(new_transcript_name)
    return sublist_transcripts_overlap_total


def implement_dict_all_ORFs_filtered_with_stop(dict_all_ORFs_filtered_with_stop : dict, dict_all_ORFs_filtered_with_stop_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the ORF names and updates the total dictionary with the new names.

    This function adds a unique prefix (the transcriptome ID) to the names of all ORFs in the 
    provided dictionary and stores the results in the total dictionary. This allows the ORF names 
    to be uniquely identified by their corresponding transcriptome.

    Parameters:
    -----------
    dict_all_ORFs_filtered_with_stop : dict
        A dictionary where the keys are ORF names (strings) and the values are the corresponding ORF sequences or information.
    
    dict_all_ORFs_filtered_with_stop_total : dict
        A dictionary where the new ORF names (with the added prefix) will be stored as keys, 
        with their corresponding ORF sequences or information as values.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each ORF name.

    Returns:
    --------
    dict_all_ORFs_filtered_with_stop_total : dict
        The updated dictionary where ORF names are prefixed with the transcriptome ID.
    """
    for orf_name in dict_all_ORFs_filtered_with_stop:
        new_orf_name = ID_transcriptome + "_" + orf_name
        dict_all_ORFs_filtered_with_stop_total[new_orf_name] = dict_all_ORFs_filtered_with_stop[orf_name]


def implement_dict_transcript_fasta_denovo(dict_transcript_fasta_denovo : dict, dict_transcript_fasta_denovo_total : dict, ID_transcriptome : str) -> None:
    """
    Adds a prefix (name of the considered transcriptome) to the transcript names and updates the total dictionary with the new names and their corresponding sequences.

    This function adds a unique prefix (the transcriptome ID) to the names of all transcripts in the 
    provided dictionary and stores the results in the total dictionary. This allows the transcript names 
    to be uniquely identified by their corresponding transcriptome.

    Parameters:
    -----------
    dict_transcript_fasta_denovo : dict
        A dictionary where the keys are transcript names (strings) and the values are the corresponding transcript sequences.
    
    dict_transcript_fasta_denovo_total : dict
        A dictionary where the new transcript names (with the added prefix) will be stored as keys, 
        with their corresponding transcript sequences as values.

    ID_transcriptome : str
        The unique identifier for the transcriptome. This ID will be added as a prefix to each transcript name.

    Returns:
    --------
    dict_transcript_fasta_denovo_total : dict
        The updated dictionary where transcript names are prefixed with the transcriptome ID.
    """
    for transcript_name in dict_transcript_fasta_denovo:
        new_transcript_name = ID_transcriptome + "_" + transcript_name
        dict_transcript_fasta_denovo_total[new_transcript_name] = dict_transcript_fasta_denovo[transcript_name]
