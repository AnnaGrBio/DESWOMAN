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


def get_file_size(my_file: str) -> dict:
    """
    Get the length of sequences in a FASTA file.

    This function parses a FASTA file and returns a dictionary with the sequence IDs as keys
    and the corresponding sequence lengths (in base pairs) as values.

    Parameters:
    -----------
    my_file : str or file-like object
        The path to the FASTA file or a file-like object containing sequences in FASTA format.

    Returns:
    --------
    dico_name_size : dict
        A dictionary where the keys are sequence IDs (strings) and the values are the lengths
        of the corresponding sequences (integers, in base pairs).
    """
    dico_name_size = {}
    for seq_record in SeqIO.parse(my_file, "fasta"):
        dico_name_size[str(seq_record.id)] = len(str(seq_record.seq))
    return dico_name_size


def retrieve_name_hit_nucl_blast(name_intermediate_directory: str) -> list:
    """
    Extract unique query names from a nucleotide BLAST output file.

    This function processes the output of a nucleotide BLAST (typically in text format),
    extracts the names of queries that have hits, and returns a list of unique query names.

    Parameters:
    -----------
    name_intermediate_directory : str
        The directory where the BLAST output file is stored. The expected file is named
        "denovo_blast_output_nucl.txt" and contains the BLAST search results.

    Returns:
    --------
    list_correct_hits : list of str
        A list of unique query names (strings) that have at least one hit, extracted from the
        BLAST output file. Duplicates are removed.
    """
    list_correct_hits = []
    name_output = name_intermediate_directory + "/denovo_blast_output_nucl.txt"
    my_file = openFile(name_output)
    if len(my_file) > 0:
        for line in my_file:
            elts_line = line.split()
            if elts_line[0] != "#":
                list_correct_hits.append(elts_line[0])
    list_correct_hits = list(set(list_correct_hits))
    return list_correct_hits


def retrieve_name_hit_prot_blast(name_intermediate_directory: str) -> list:
    """
    Extract unique query names from a protein BLAST (DIAMOND) output file.

    This function processes the output of a protein DIAMOND BLAST (typically in text format),
    extracts the names of queries that have hits, and returns a list of unique query names.

    Parameters:
    -----------
    name_intermediate_directory : str
        The directory where the BLAST output file is stored. The expected file is named
        "denovo_blast_output_prot.txt" and contains the DIAMOND BLAST search results.

    Returns:
    --------
    list_correct_hits : list of str
        A list of unique query names (strings) that have at least one hit, extracted from the
        BLAST output file. Duplicates are removed.
    """
    list_correct_hits = []
    name_output = name_intermediate_directory + "/denovo_blast_output_prot.txt"
    my_file = openFile(name_output)
    if len(my_file) > 0:
        for line in my_file:
            elts_line = line.split()
            if elts_line[0] not in list_correct_hits:
                list_correct_hits.append(elts_line[0])
    list_correct_hits = list(set(list_correct_hits))
    return list_correct_hits


def merge_hit_lists(list1: list, list2: list) -> list:
    """
    Merge two lists of hit names and return a list of unique values.

    This function combines two input lists (`list1` and `list2`) and removes any duplicate
    entries, ensuring that only unique hit names remain in the final list.

    Parameters:
    -----------
    list1 : list of str
        The first list of hit names to be merged. Each element in the list should be a string
        representing a hit name.

    list2 : list of str
        The second list of hit names to be merged. Each element in the list should be a string
        representing a hit name.

    Returns:
    --------
    list_final : list of str
        A list of unique hit names resulting from the union of `list1` and `list2`. Any duplicates
        from the two input lists will be removed.
    """
    list_final = list(set(list1 + list2))
    return list_final


def reshufe_files_in_denovo(
    name_output_directory: str, list_denovo_to_delete: list
) -> None:
    """
    Recreate de novo generated output files by removing entries that have a BLAST hit.

    This function processes multiple output files generated during a de novo analysis pipeline.
    It removes sequences (from FASTA files and a tabular information file) that correspond to entries
    in the list of de novo sequences that have associated BLAST hits. The remaining sequences are then
    written back into new files, preserving the original structure but excluding those with hits.

    Parameters:
    -----------
    name_output_directory : str
        The directory where the output files from the de novo analysis are stored. These files include
        the denovo nucleic acid and protein sequences, as well as an information file with related metadata.

    list_denovo_to_delete : list of str
        A list containing the names (IDs) of de novo sequences that should be removed, based on having
        a corresponding BLAST hit. Sequences in the list will not be included in the recreated output files.

    Modifies:
    --------
    The following files in the `name_output_directory` will be modified:
    - "denovo_nucl.fa" (nucleotide sequences)
    - "denovo_unspliced_lowered_introns.fa" (if it exists, nucleotide sequences with introns)
    - "denovo_protein.fa" (protein sequences)
    - "information_file.txt" (metadata associated with sequences)
    """
    # denovo nucl
    dico_new_denovo_nucl = {}
    denovo_nucl_file = name_output_directory + "/denovo_nucl.fa"
    # store all denovo with no hit
    for seq_record in SeqIO.parse(denovo_nucl_file, "fasta"):
        if str(seq_record.id) not in list_denovo_to_delete:
            dico_new_denovo_nucl[">" + str(seq_record.id) + "\n"] = (
                str(seq_record.seq) + "\n"
            )
    # destroy and recreate denovo file
    os.system("rm " + denovo_nucl_file)
    my_file = open(denovo_nucl_file, "w")
    for name in dico_new_denovo_nucl:
        my_file.write(name)
        my_file.write(dico_new_denovo_nucl[name])
    my_file.close()

    # denovo nucl lowered intron
    denovo_unspliced_lowered_introns_file = (
        name_output_directory + "/denovo_unspliced_lowered_introns.fa"
    )
    # make sure an intron file exist (which is only the case for strat1)
    if os.path.isfile(denovo_unspliced_lowered_introns_file):
        dico_new_denovo_nucl_lowered_intron = {}
        for seq_record in SeqIO.parse(denovo_unspliced_lowered_introns_file, "fasta"):
            if str(seq_record.id) not in list_denovo_to_delete:
                dico_new_denovo_nucl_lowered_intron[">" + str(seq_record.id) + "\n"] = (
                    str(seq_record.seq) + "\n"
                )
        # destroy and recreate denovo file
        os.system("rm " + denovo_unspliced_lowered_introns_file)
        my_file = open(denovo_unspliced_lowered_introns_file, "w")
        for name in dico_new_denovo_nucl_lowered_intron:
            my_file.write(name)
            my_file.write(dico_new_denovo_nucl_lowered_intron[name])
        my_file.close()

    # denovo prot
    dico_new_denovo_prot = {}
    denovo_protein_file = name_output_directory + "/denovo_protein.fa"
    for seq_record in SeqIO.parse(denovo_protein_file, "fasta"):
        if str(seq_record.id) not in list_denovo_to_delete:
            dico_new_denovo_prot[">" + str(seq_record.id) + "\n"] = (
                str(seq_record.seq) + "\n"
            )
    os.system("rm " + denovo_protein_file)
    my_file = open(denovo_protein_file, "w")
    for name in dico_new_denovo_prot:
        my_file.write(name)
        my_file.write(dico_new_denovo_prot[name])
    my_file.close()

    # info_file
    list_info_file = []
    information_file = name_output_directory + "/information_file.txt"
    study_file = openFile(information_file)
    for line in study_file:
        name_orf = line.split(",")[9]
        if name_orf not in list_denovo_to_delete:
            list_info_file.append(line)
    os.system("rm " + information_file)
    my_file = open(information_file, "w")
    for i in list_info_file:
        my_file.write(i)
    my_file.close()
