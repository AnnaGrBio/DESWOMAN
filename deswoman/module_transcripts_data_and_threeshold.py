import os
from Bio import SeqIO
from deswoman.module_colors import *


__author__ = "Anna Grandchamp"
__contributor__ = "Marie Lebherz"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def get_list_name_query_rep(path_rep: str) -> list:
    """
    This function retrieves unique query names from filenames in the specified directory.
    The query name is taken as the part of the filename before the first dot (".").

    Parameters:
    -----------
    path_rep : str
        The path to the directory containing the files. The filenames will be processed
        to extract the query names.

    Returns:
    --------
    list_query_name : list
        A list of unique query names derived from the filenames in the specified directory.
        Each query name corresponds to the part before the first dot (".") in the filename.
    """
    list_query_name = []
    elts_dir = os.listdir(path_rep)
    for names in elts_dir:
        header = names.split(".")[0]
        if header not in list_query_name:
            list_query_name.append(header)
    return list_query_name


def get_path_gtf(path_rep: str, query_name: str) -> str:
    """
    Returns the path to the GTF file corresponding to the given query name.

    This function constructs the file path for a transcriptome annotation file
    (in GTF format) based on the provided query name and the directory where
    the files are stored. It assumes that the GTF file is named using the format
    "<query_name>.gtf" and is located within the specified directory.

    Parameters:
    -----------
    path_rep : str
        The directory where the GTF file is located.

    query_name : str
        The base name of the query (without file extension) used to search for
        the corresponding GTF file.

    Returns:
    --------
    new_path : str
        The complete file path to the GTF file corresponding to the provided
        query name.
    """
    searched_name = query_name + "." + "gtf"
    elts_dir = os.listdir(path_rep)
    # this does not include a condition because the presence of the query file in the folder was tested in the validation step
    new_path = path_rep + "/" + searched_name
    return new_path


def get_path_to_fasta(path_rep: str, query_name: str) -> str:
    """
    Returns the path to a FASTA file corresponding to the given query name.

    This function searches the specified directory for a file that starts with
    the given query name and ends with one of the recognized FASTA file extensions
    (".fna", ".fa", or ".fasta"). It returns the full path to the first matching
    file found.

    Parameters:
    -----------
    path_rep : str
        The directory where the FASTA file is located.

    query_name : str
        The base name of the query (without file extension) used to search for
        the corresponding FASTA file.

    Returns:
    --------
    new_path : str
        The complete file path to the FASTA file corresponding to the provided
        query name.
    """
    searched_name = ""
    elts_dir = os.listdir(path_rep)
    for file in elts_dir:
        lst = file.split(".")
        if lst[0] == query_name:
            # searches for the 3 output file termination
            if lst[1] == "fna" or lst[1] == "fa" or lst[1] == "fasta":
                searched_name = file
                break
    new_path = path_rep + "/" + searched_name
    return new_path


def get_path_gff(path_rep: str, query_name: str) -> str:
    """
    Returns the path to a GFF file corresponding to the given query name.

    This function searches the specified directory for a file that starts with
    the given query name and ends with one of the recognized GFF file extensions
    (".gff" or ".gff3"). It returns the full path to the first matching
    file found.

    Parameters:
    -----------
    path_rep : str
        The directory where the GFF file is located.

    query_name : str
        The base name of the query (without file extension) used to search for
        the corresponding GFF file.

    Returns:
    --------
    new_path : str
        The complete file path to the GFF file corresponding to the provided
        query name.
    """
    searched_name = query_name + "." + "gff"
    elts_dir = os.listdir(path_rep)
    for file in elts_dir:
        lst = file.split(".")
        if lst[0] == query_name:
            # searches for the 2 output file termination
            if lst[1] == "gff" or lst[1] == "gff3":
                searched_name = file
                break
    new_path = path_rep + "/" + searched_name
    return new_path


def assessTE(fasta_seq: str) -> bool:
    """
    Assesses whether there are lowercase letters in the given FASTA sequence.

    This function checks if the input FASTA sequence contains any lowercase
    nucleotides, which may indicate the presence of transposable elements (TEs).
    If any lowercase nucleotide is found, the function returns True.
    Otherwise, it returns False.

    Parameters:
    -----------
    fasta_seq : str
        A string representing a nucleotide sequence in FASTA format.

    Returns:
    --------
    presenceTE : bool
        True if any lowercase nucleotide is found in the sequence, indicating
        the presence of transposable elements (TEs); False otherwise.
    """
    presenceTE = False
    for nucl in fasta_seq:
        if nucl.islower():
            presenceTE = True
            break
    return presenceTE


def exclude_transcripts(dict_transcript_fasta_all: dict, te_exclusion: str) -> list:
    """
    Excludes transcripts that overlap with transposable elements (TEs) based on lowercase letters in the sequence.

    If the option to remove transcripts overlapping with TEs is activated (`te_exclusion == "True"`), this function checks each transcript in the
    provided dictionary (`dict_transcript_fasta_all`). It identifies transcripts that contain lowercase letters, which are assumed to
    correspond to annotated TEs, and stores these transcript names in a list. This list is then returned.

    Parameters:
    -----------
    dict_transcript_fasta_all : dict
        A dictionary where keys are transcript names and values are corresponding nucleotide sequences in FASTA format.

    te_exclusion : str
        A string that indicates whether to exclude transcripts overlapping with TEs. If set to `"True"`, transcripts containing lowercase
        letters (indicating TEs) will be excluded.

    Returns:
    --------
    list_transcripts_to_exclude : list
        A list of transcript names that overlap with TEs (contain lowercase letters in their sequence).
    """
    list_transcripts_to_exclude = []
    if te_exclusion == "True":
        for transcript_name in dict_transcript_fasta_all:
            presenceTE = assessTE(dict_transcript_fasta_all[transcript_name])
            if presenceTE == True:
                list_transcripts_to_exclude.append(transcript_name)
    return list_transcripts_to_exclude


def detect_tpm_in_file(elts_line: list) -> int | None:
    """
    Detects the column in a GTF file that contains TPM (Transcripts Per Million) values.

    This function scans through a list of elements from a line in a GTF file to identify
    the column that holds the TPM values. It returns the column number where TPM values
    are found, or `None` if no such column exists.

    Parameters:
    -----------
    elts_line : list
        A list of strings representing the elements (columns) of a single line in a GTF file.

    Returns:
    --------
    num_column_tpm : int or None
        The column number (1-indexed) where the TPM values are found, or `None` if no TPM column is present.
    """
    num_column_tpm = None
    nb_column = 0
    for elt in elts_line:
        if elt == "TPM":
            num_column_tpm = nb_column + 1
            break
        nb_column += 1
    return num_column_tpm


def detect_fpkm_in_file(elts_line: list) -> int | None:
    """
    Detects the column in a GTF file where the FPKM (Fragments Per Kilobase of transcript per Million mapped reads) values are located.

    Args:
        elts_line (list of str): A list representing a single line from a GTF file, where each element corresponds to a column in the file.

    Returns:
        int or None: The 1-based index of the column containing FPKM values, or None if no such column is found.
    """
    num_column_fpkm = None
    nb_column = 0
    for elt in elts_line:
        if elt == "FPKM":
            num_column_fpkm = nb_column + 1
            break
        nb_column += 1
    return num_column_fpkm


def detect_cov_in_file(elts_line: list) -> int | None:
    """
    Detects the column in a GTF file where the coverage values are located.

    Args:
        elts_line (list of str): A list representing a single line from a GTF file, where each element corresponds to a column in the file.

    Returns:
        int or None: The 1-based index of the column containing coverage values, or None if no such column is found.
    """
    num_column_cov = None
    nb_column = 0
    for elt in elts_line:
        if elt == "cov":
            num_column_cov = nb_column + 1
            break
        nb_column += 1
    return num_column_cov


def associate_gene_and_transcripts(
    opened_gtf_file: list,
    min_tpm_asked: int,
    dict_transcript_fasta_all: dict,
    te_exclusion: str,
    unoriented_transc_exclusion: str,
) -> dict:
    """
    Associates genes with their corresponding transcripts based on a GTF file of transcriptome assembly.

    This function processes a GTF file to create a dictionary where genes are the keys and the values are lists of their corresponding transcript IDs.
    It also filters transcripts based on TPM values and exclusions such as transcripts overlapping with TEs or those with unspecified orientation.

    Args:
        opened_gtf_file (iterable): An iterable object (e.g., a file object) representing the opened GTF file.
        min_tpm_asked (float): The minimum TPM value to filter the transcripts. Only transcripts with TPM greater than or equal to this value will be included.
        dict_transcript_fasta_all (dict): A dictionary of transcript sequences (used for TE exclusion).
        te_exclusion (str): Whether to exclude transcripts overlapping with TEs ('True' or 'False').
        unoriented_transc_exclusion (str): Whether to exclude transcripts with unspecified orientation ('True' or 'False').

    Returns:
        dict: A dictionary with gene IDs as keys and lists of transcript IDs as values for genes that meet the filtering criteria.
    """
    # if the user wants to filters transcripts overlapping with annotated TEs (in lower letter in genome), it happens just below
    list_transcripts_to_exclude = exclude_transcripts(
        dict_transcript_fasta_all, te_exclusion
    )
    dico_genes = {}

    for line in opened_gtf_file:
        # Check if the line is not a comment
        if line[0] != "#":
            # Split the line into elements
            elts_line = line.split()
            if len(elts_line) > 9:
                # Check if the element is a transcript
                if (
                    elts_line[2] == "transcript"
                ):  # Note: Handle other applications if needed
                    # Extract gene_id, transcript_id, and TPM from the line
                    gene_id = elts_line[9][
                        1 : len(elts_line[9]) - 2
                    ]  # Note: Ensure this works for all cases
                    transcript_id = elts_line[11][
                        1 : len(elts_line[11]) - 2
                    ]  # Note: Ensure this works for all cases
                    num_column_tpm = detect_tpm_in_file(
                        elts_line
                    )  # access the column where the TPM values are
                    TPM = elts_line[num_column_tpm][
                        1 : len(elts_line[num_column_tpm]) - 2
                    ]  # Note: Ensure this works for all cases
                    transcript_orientation = elts_line[6]
                    # exclude transcripts with "." as orientation is the user asked for it
                    if (
                        unoriented_transc_exclusion == "True"
                        and transcript_orientation == "."
                    ):
                        list_transcripts_to_exclude.append(transcript_id)

                    # Check if the TPM value is greater than or equal to the specified minimum TPM
                    if (
                        float(TPM) >= float(min_tpm_asked)
                        and transcript_id not in list_transcripts_to_exclude
                    ):
                        # Check if the gene_id is not in the dictionary, then add it
                        if gene_id not in dico_genes:
                            dico_genes[gene_id] = [transcript_id]
                        else:
                            # If the gene_id is already in the dictionary, append the transcript_id to its list
                            dico_genes[gene_id].append(transcript_id)

    # Return the dictionary of gene-transcript associations
    return dico_genes


def extract_transcripts_properties_from_gff(
    dico_transcripts: dict, elts_line: list
) -> None:
    """
    Extracts properties of a transcript from a line in a GFF file and stores the information in a dictionary.

    This function processes a line from a GFF file containing transcript information, extracting properties such as the chromosome, direction,
    start and stop positions, TPM, FPKM, and coverage. The extracted information is stored in the provided dictionary `dico_transcripts`
    using the transcript ID as the key.

    Args:
        dico_transcripts (dict): A dictionary to store the transcript information. The transcript ID will be used as the key,
                                  and the value will be a list containing the transcript's properties (chromosome, direction,
                                  start, stop, FPKM, TPM, coverage).
        elts_line (list): A list of elements representing a line from a GFF file. This line is assumed to contain information about a transcript.

    Returns:
        None: The function modifies the `dico_transcripts` dictionary in place by adding a new entry for the transcript.
    """
    chromosome = elts_line[0]
    direction = elts_line[6]
    start = elts_line[3]
    stop = elts_line[4]
    transcript_id = elts_line[11][
        1 : len(elts_line[11]) - 2
    ]  # Note: Ensure this works for all cases
    nb_coulumn_tpm = detect_tpm_in_file(elts_line)  # extract pos of the TPM
    nb_coulumn_fpkm = detect_fpkm_in_file(elts_line)  # extract pos of the FPKM
    nb_coulumn_cov = detect_cov_in_file(elts_line)  # extract pos of the coverage
    FPKM = elts_line[nb_coulumn_fpkm][
        1 : len(elts_line[nb_coulumn_fpkm]) - 2
    ]  # Note: Ensure this works for all cases
    TPM = elts_line[nb_coulumn_tpm][
        1 : len(elts_line[nb_coulumn_tpm]) - 2
    ]  # Note: Ensure this works for all cases
    coverage = elts_line[nb_coulumn_cov][
        1 : len(elts_line[nb_coulumn_cov]) - 2
    ]  # Note: Ensure this works for all cases

    # Add transcript information to the dictionary
    dico_transcripts[transcript_id] = [
        chromosome,
        direction,
        start,
        stop,
        FPKM,
        TPM,
        coverage,
    ]


def get_transcripts_properties(opened_gtf_file: list) -> dict:
    """
    Processes a GTF file to extract and organize transcript properties into a dictionary.

    This function reads through each line of the provided GTF file, extracting properties for both transcripts and exons.
    It stores these properties in a dictionary where each transcript ID is the key, and the value is a list containing the transcript's
    properties such as chromosome, direction, start and stop positions, FPKM, TPM, coverage, and the number of exons.

    The function also creates a dictionary for inter-exon information to track the maximum number of exons for each transcript, which
    is then added to the main dictionary.

    Args:
        opened_gtf_file (file object): The opened GTF file containing transcript and exon data.

    Returns:
        dict: A dictionary where each key is a transcript ID and each value is a list of transcript properties, including the
              number of exons.
    """
    # Initialize a dictionary with the first key and item lists
    dico_transcripts = {
        "key_transcript_name": [
            "chromosome",
            "direction",
            "start",
            "stop",
            "FPKM",
            "TPM",
            "coverage",
            "nb_exons",
        ]
    }

    # Initialize an empty dictionary for inter-exon information
    dico_inter_exons = {}

    # Iterate through each line in the opened gtf file
    for line in opened_gtf_file:
        # Check if the line is not a comment
        if line[0] != "#":
            # Split the line into elements
            elts_line = line.split()

            # Check if the element is a transcript
            if (
                elts_line[2] == "transcript"
            ):  # Note: Handle other applications if needed
                extract_transcripts_properties_from_gff(dico_transcripts, elts_line)

            # Check if the element is an exon
            elif elts_line[2] == "exon":  # Note: Handle other applications if needed
                # Extract properties from the line for the exon
                transcript_id = elts_line[11][
                    1 : len(elts_line[11]) - 2
                ]  # Note: Ensure this works for all cases
                exon_nb = int(
                    elts_line[13][1 : len(elts_line[13]) - 2]
                )  # Note: Ensure this works for all cases

                # Add or update the number of exons for the transcript in the inter-exon dictionary
                if transcript_id not in dico_inter_exons:
                    dico_inter_exons[transcript_id] = exon_nb
                else:
                    if exon_nb > dico_inter_exons[transcript_id]:
                        dico_inter_exons[transcript_id] = exon_nb

    # implement the number of exons the main dictionary
    for transcripts in dico_transcripts:
        if transcripts in dico_inter_exons:
            dico_transcripts[transcripts].append(dico_inter_exons[transcripts])

    # Return the dictionary of transcript properties
    return dico_transcripts


def store_transcriptome_seq(opened_transcriptome_file: list) -> dict:
    """
    Parses a transcriptome FASTA file and stores each transcript in a dictionary.

    This function reads a FASTA file containing transcript sequences and stores them in a dictionary. The dictionary has transcript
    names as keys and their corresponding sequences as values.

    Args:
        opened_transcriptome_file (file object): The opened FASTA file containing the transcriptome sequences.

    Returns:
        dict: A dictionary where each key is the transcript name (ID) and each value is the transcript's nucleotide sequence.
    """
    dico_transcripts = {}
    for seq_record in SeqIO.parse(opened_transcriptome_file, "fasta"):
        dico_transcripts[str(seq_record.id)] = str(seq_record.seq)
    return dico_transcripts
