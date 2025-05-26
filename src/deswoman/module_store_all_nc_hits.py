import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from deswoman.module_colors import openFile


__author__ = "Anna Grandchamp"
__contributor__ = ""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def make_dico_len_unspliced_denovo(name_output_directory: str) -> dict:
    """
    Build a dictionary mapping the names of de novo sequences to the lengths of their sequences
    (including introns if present).

    This function reads a FASTA file containing de novo sequences with introns and computes
    the length of each sequence. It stores these lengths in a dictionary where the key is the
    sequence name (identifier) and the value is the length of the corresponding sequence.

    Parameters:
    -----------
    name_output_directory : str
        The directory path where the de novo sequences file is located. This function expects
        the presence of a file named "denovo_unspliced_lowered_introns.fa" in the directory,
        which contains the sequences with introns included.

    Returns:
    --------
    dico_name_size : dict
        A dictionary where the keys are the names (IDs) of the de novo sequences, and the values
        are the lengths of these sequences, including any introns.
    """
    my_denovo_file = name_output_directory + "/denovo_unspliced_lowered_introns.fa"
    dico_name_size = {}
    for seq_record in SeqIO.parse(my_denovo_file, "fasta"):
        dico_name_size[str(seq_record.id)] = len(str(seq_record.seq))
    return dico_name_size


def extract_all_blast_outputs(
    sp_pop_ind_name: str, name_intermediate_directory: str
) -> dict:
    """
    This function processes the BLAST output for denovo genes and organizes the hits into a dictionary
    where each de novo gene is associated with its corresponding BLAST hits.

    The BLAST hits are formatted as a range (start-end) on the hit sequence, along with the hit's
    target gene name and strand orientation (forward or reverse).

    Parameters:
    -----------
    sp_pop_ind_name : str
        The name of the species/population/individual combination used to identify the specific
        BLAST output file. It is used to form the filename for the BLAST result.

    name_intermediate_directory : str
        The directory path where the intermediate BLAST output files are stored.

    Returns:
    --------
    dict_my_blast_selected_outputs : dict
        A dictionary where the keys are denovo gene names (from the query in the BLAST results),
        and the values are lists of corresponding BLAST hit names. Each hit name is formatted as
        "start-end-target_gene-strand", with the strand indicated by either "-f" (forward) or "-r" (reverse).
    """
    dict_my_blast_selected_outputs = {}
    my_file_name = (
        name_intermediate_directory
        + "/blast_denovo_to_target_genome/"
        + sp_pop_ind_name
        + "_all_bast_output.txt"
    )
    my_blast_results = openFile(my_file_name)
    for line in my_blast_results:
        if line[0] != "#":
            list_elts = line.split()
            denovo_name = list_elts[0]
            if int(list_elts[5]) < int(list_elts[6]):
                my_hit_name = (
                    list_elts[5] + "-" + list_elts[6] + "-" + list_elts[1] + "-f"
                )
            else:
                my_hit_name = (
                    list_elts[6] + "-" + list_elts[5] + "-" + list_elts[1] + "-r"
                )
            if denovo_name not in dict_my_blast_selected_outputs:
                dict_my_blast_selected_outputs[denovo_name] = [my_hit_name]
            else:
                dict_my_blast_selected_outputs[denovo_name].append(my_hit_name)
    return dict_my_blast_selected_outputs


def extract_best_blast_outputs(
    sp_pop_ind_name: str, name_intermediate_directory: str
) -> dict:
    """
    This function processes BLAST results to associate the best homologous hit to each de novo gene
    identified in the target genome. The function assumes that the user has opted not to use synteny
    for gene detection, and instead, it selects the best hit based on the BLAST results.

    It works by selecting the first hit for each de novo gene, which is considered the best hit
    based on the BLAST results (the first one encountered in the file).

    Parameters:
    -----------
    sp_pop_ind_name : str
        The name of the species/population/individual combination used to identify the specific
        BLAST output file. It is used to form the filename for the BLAST result.

    name_intermediate_directory : str
        The directory path where the intermediate BLAST output files are stored.

    Returns:
    --------
    dict_my_blast_selected_outputs : dict
        A dictionary where the keys are denovo gene names (from the query in the BLAST results),
        and the values are the best BLAST hit for each denovo gene, formatted as
        "start-end-target_gene-strand" with the strand indicated as either "-f" (forward) or "-r" (reverse).
    """
    dict_my_blast_selected_outputs = {}
    my_file_name = (
        name_intermediate_directory
        + "/blast_denovo_to_target_genome/"
        + sp_pop_ind_name
        + "_all_bast_output.txt"
    )
    my_blast_results = openFile(my_file_name)
    for line in my_blast_results:
        if line[0] != "#":
            list_elts = line.split()
            denovo_name = list_elts[0]
            if int(list_elts[5]) < int(list_elts[6]):
                # detect forward direction in the genome
                my_hit_name = (
                    list_elts[5] + "-" + list_elts[6] + "-" + list_elts[1] + "-f"
                )
            else:
                # detect reverse direction in the genome
                my_hit_name = (
                    list_elts[6] + "-" + list_elts[5] + "-" + list_elts[1] + "-r"
                )
            # take only the first denovo met, which is the best hit, in case of several occurences
            if denovo_name not in dict_my_blast_selected_outputs:
                dict_my_blast_selected_outputs[denovo_name] = my_hit_name
    return dict_my_blast_selected_outputs


def store_de_novo_informations(name_output_directory: str) -> dict:
    """
    This function extracts key information about de novo genes from the `information_file.txt`
    output file and stores it in a dictionary.

    The `information_file.txt` contains multiple columns with details about each de novo gene,
    and this function specifically extracts the following information:
    - Gene name (from column 9)
    - Genomic location (from column 12)
    - Chromosome (from column 1)
    - Start and stop positions (from columns 6 and 7 respectively)

    Parameters:
    -----------
    name_output_directory : str
        The directory path where the `information_file.txt` is stored.

    Returns:
    --------
    dico_denovo_info : dict
        A dictionary where the keys are de novo gene names (from column 9),
        and the values are lists containing the following information:
        - Genomic location (column 12)
        - Chromosome (column 1)
        - Start position (column 6)
        - Stop position (column 7)
    """
    dico_denovo_info = {}
    file_name = name_output_directory + "/information_file.txt"
    my_file = openFile(file_name)
    for line in my_file[1:]:
        line = line.split("\n")[0]
        name_denovo = line.split(",")[9]
        genomic_location = line.split(",")[12]
        chromosome = line.split(",")[1]
        start = line.split(",")[6]
        stop = line.split(",")[7]
        dico_denovo_info[name_denovo] = [genomic_location, chromosome, start, stop]
    return dico_denovo_info
