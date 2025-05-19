import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_graphical_interface_strat1 import my_graphical_interface_strategy1
from module_handle_config_file import my_config_file_extract_parameters
from module_colors import *


__author__ = "Anna Grandchamp"
__contributor__="Marie Lebherz"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"
    

def validate_presence_of_mandatory_parameters_strat1(dico_variables : dict, RUN_PYTHON : bool) -> bool:
    """
    Validates the presence of mandatory parameters in the provided dictionary.

    This function checks if essential keys in `dico_variables` contain valid (non-empty) values.
    If any mandatory parameter is missing or empty, an error message is printed, and the `RUN_PYTHON`
    flag is set to False.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing the necessary parameters. Expected keys:
        - "query": Name of the query.
        - "path_to_genome_repository": Path to the genome repository.
        - "path_to_transcriptome_repository": Path to the transcriptome repository.

    RUN_PYTHON : bool
        A flag indicating whether the script should continue execution.

    Returns:
    --------
    bool
        Updated `RUN_PYTHON` flag (False if any mandatory parameter is missing).
    """
    if dico_variables["query"] == "":
        print (BRIGHT_RED + "ERROR... NO QUERY NAME PROVIDED ... " + RESET)
        RUN_PYTHON = False
    if dico_variables["path_to_genome_repository"] == "":
        print (BRIGHT_RED + "ERROR... NO GENOME FOLDER PROVIDED ... " + RESET)
        RUN_PYTHON = False
    if dico_variables["path_to_transcriptome_repository"] == "":
        print (BRIGHT_RED + "ERROR... NO TRANSCRIPTOME FOLDER PROVIDED ... " + RESET)
        RUN_PYTHON = False
    return RUN_PYTHON


def validate_presence_query_genome_strat1(query_name : str, path_rep : str) -> bool:
    """
    Checks whether the specified query FASTA genome file is present in the given genome folder.

    This function scans the directory at `path_rep` to verify the existence of a FASTA file 
    corresponding to `query_name`. It considers common FASTA file extensions: "fa", "fasta", and "fna".

    Parameters:
    -----------
    query_name : str
        The expected name of the query genome file (without the file extension).

    path_rep : str
        The path to the genome folder where the FASTA file should be located.

    Returns:
    --------
    bool
        True if a FASTA file matching `query_name` is found in `path_rep`, otherwise False.
    """
    elts_dir = os.listdir(path_rep)
    list_fasta_naming = ["fa", "fasta", "fna"]
    presence_query = False
    for names in elts_dir:
        if (".") in names:
            header = names.split(".")[0]
            format = names.split(".")[1]
            if header == query_name and format in list_fasta_naming:
                presence_query = True
                break
    return presence_query


def validate_presence_query_gff_strat1(query_name : str, path_rep : str) -> bool:
    """
    Checks whether the specified query GFF genome file is present in the given genome folder.

    This function scans the directory at `path_rep` to verify the existence of a GFF file 
    corresponding to `query_name`. It considers common GFF file extensions: "gff" and "gff3".

    Parameters:
    -----------
    query_name : str
        The expected name of the query genome file (without the file extension).

    path_rep : str
        The path to the genome folder where the GFF file should be located.

    Returns:
    --------
    bool
        True if a GFF file matching `query_name` is found in `path_rep`, otherwise False.
    """
    elts_dir = os.listdir(path_rep)
    list_gff_naming = ["gff", "gff3"]
    presence_query = False
    for names in elts_dir:
        if (".") in names:
            header = names.split(".")[0]
            format = names.split(".")[1]
            if header == query_name and format in list_gff_naming:
                presence_query = True
                break
    return presence_query


def validate_presence_query_transcriptome_strat1(query_name : str, path_rep : str) -> bool:
    """
    Checks whether the specified query FASTA transcriptome file is present in the given transcriptome folder.

    This function scans the directory at `path_rep` to verify the existence of a FASTA file 
    corresponding to `query_name`. It considers common FASTA file extensions: "fa", "fasta", and "fna".

    Parameters:
    -----------
    query_name : str
        The expected name of the query transcriptome file (without the file extension).

    path_rep : str
        The path to the transcriptome folder where the FASTA file should be located.

    Returns:
    --------
    bool
        True if a FASTA file matching `query_name` is found in `path_rep`, otherwise False.
    """
    elts_dir = os.listdir(path_rep)
    list_fasta_naming = ["fa", "fasta", "fna"]
    presence_query = False
    for names in elts_dir:
        if (".") in names:
            header = names.split(".")[0]
            format = names.split(".")[1]
            if header == query_name and format in list_fasta_naming:
                presence_query = True
                break
    return presence_query


def validate_presence_query_gtf_strat1(query_name : str, path_rep : str) -> bool:
    """
    Checks whether the specified query GTF genome file is present in the given transcriptome folder.

    This function scans the directory at `path_rep` to verify the existence of a GTF file 
    corresponding to `query_name`. It ensures that the file has a `.gtf` extension.

    Parameters:
    -----------
    query_name : str
        The expected name of the query genome file (without the file extension).

    path_rep : str
        The path to the transcriptome folder where the GTF file should be located.

    Returns:
    --------
    bool
        True if a GTF file matching `query_name` is found in `path_rep`, otherwise False.
    """
    elts_dir = os.listdir(path_rep)
    presence_query = False
    for names in elts_dir:
        if (".") in names:
            header = names.split(".")[0]
            format = names.split(".")[1]
            if header == query_name and format == "gtf":
                presence_query = True
                break
    return presence_query
    

def validate_presence_query_name_strat1(dico_variables : dict, RUN_PYTHON : bool) -> bool:
    """
    Validates the presence of query genome and transcriptome data in the specified directories.

    This function checks whether the required genome and transcriptome files exist in their respective 
    folders. It verifies the presence of:
    - A genome FASTA file (.fasta, .fa, or .fna) in the genome repository.
    - A genome annotation file (.gff or .gff3) in the genome repository.
    - A transcriptome FASTA file (.fasta, .fa, or .fna) in the transcriptome repository.
    - A transcriptome annotation file (.gtf) in the transcriptome repository.

    If any of the required files are missing, an error message is displayed, and `RUN_PYTHON` is set to `False`.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing the necessary paths and query name. Expected keys:
        - "query": The name of the query sequence.
        - "path_to_genome_repository": Path to the genome repository.
        - "path_to_transcriptome_repository": Path to the transcriptome repository.

    RUN_PYTHON : bool
        A flag indicating whether the script should continue execution.

    Returns:
    --------
    bool
        Updated `RUN_PYTHON` flag (False if any required file is missing).
    """
    query_name = dico_variables["query"]
    path_to_genome = dico_variables["path_to_genome_repository"]
    path_to_transcriptome = dico_variables["path_to_transcriptome_repository"]

    # validate that a genome in fasta and  format exit in the genome folder with the query name
    presence_fasta_genome = validate_presence_query_genome_strat1(query_name, path_to_genome) 
    presence_gff = validate_presence_query_gff_strat1(query_name, path_to_genome)

    # validate that a transcriptome in fasta and  format exit in the genome folder with the query name
    presence_fasta_transcriptome = validate_presence_query_transcriptome_strat1(query_name, path_to_transcriptome)
    presence_fasta_gtf = validate_presence_query_gtf_strat1(query_name, path_to_transcriptome)

    # Present corresponding error message if one of the validation steps is not satistied
    if presence_fasta_genome == False:
        RUN_PYTHON = False
        name_fasta_genome = query_name + (".fasta")
        print (BRIGHT_RED + "ERROR : Query genome " + name_fasta_genome + " (or .fa|.fna) not present in genome folder" + RESET)
    if presence_gff == False:
        RUN_PYTHON = False
        name_gff_genome = query_name + (".gff")
        print (BRIGHT_RED + "ERROR : Query genome annotation " + name_gff_genome + " (or .gff3) not present in genome folder" + RESET)
    if presence_fasta_transcriptome == False:
        RUN_PYTHON = False
        name_fasta_transcriptome = query_name + (".fasta")
        print (BRIGHT_RED + "ERROR : Query transcriptome " + name_fasta_transcriptome + " (or .fa|.fna) not present in transcriptome folder" + RESET)
    if presence_fasta_gtf == False:
        RUN_PYTHON = False
        name_gtf_transcriptome = query_name + (".gtf")
        print (BRIGHT_RED + "ERROR : Query transcriptome annotation " + name_gtf_transcriptome + " not present in transcriptome folder" + RESET)
    return RUN_PYTHON


def check_genomic_overlap_strat1(dico_variables : dict, RUN_PYTHON : bool) -> bool:
    """
    Checks if genomic locations of the denovo selected by the user (intergnic, etc) are provided.
    This because by mistake the ser could unselect all options on the graphical interface.

    This function verifies whether the "transcript_overlap" key in `dico_variables` contains any data. 
    If it is empty, an error message is displayed, and `RUN_PYTHON` is set to `False`.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing various input parameters. Expected key:
        - "transcript_overlap": A list or other iterable storing transcript location data.

    RUN_PYTHON : bool
        A flag indicating whether the script should continue execution.

    Returns:
    --------
    bool
        Updated `RUN_PYTHON` flag (False if no transcript locations are provided).
    """
    if len(dico_variables["transcript_overlap"]) == 0:
        RUN_PYTHON = False
        print (BRIGHT_RED + "ERROR : No transcript location provided" + RESET)
    return RUN_PYTHON


def is_fasta_strat1(filename : str) -> bool:
    """
    Checks whether a file is in FASTA format.

    This function attempts to parse the given file using the Biopython `SeqIO.parse()` function 
    with the "fasta" format. If the file contains at least one valid FASTA record, the function 
    returns `True`; otherwise, it returns `False`.

    Parameters:
    -----------
    filename : str
        Path to the file to be checked.

    Returns:
    --------
    bool
        True if the file is a valid FASTA file (contains at least one FASTA record), otherwise False.
    """
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
    
    
def is_nucl_strat1(filename : str) -> bool:
    """
    Checks whether a FASTA file contains valid nucleotide sequences.

    This function parses a FASTA file and verifies if the sequence consists only of valid 
    nucleotide characters (including standard bases and ambiguous IUPAC nucleotide codes).

    Parameters:
    -----------
    filename : str
        Path to the FASTA file to be checked.

    Returns:
    --------
    bool
        True if the first sequence in the file consists only of valid nucleotide characters, 
        otherwise False.
    """
    list_nucl = ["A", "C", "G", "T", "U", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N"]
    is_n = True
    dico_seq = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        dico_seq[str(seq_record.id)] = str(seq_record.seq)
    for name in dico_seq:
        my_seq = dico_seq[name]
        for nucl in my_seq:
            new_nucl = nucl.upper()
            if new_nucl not in list_nucl:
                is_n = False
        break
    return is_n


def is_prot_strat1(filename : str) -> bool:
    """
    Checks whether a FASTA file contains valid protein sequences.

    This function parses a FASTA file and verifies if the sequence consists only of valid 
    amino acid characters (including standard amino acids and the stop codon "*").

    Parameters:
    -----------
    filename : str
        Path to the FASTA file to be checked.

    Returns:
    --------
    bool
        True if the first sequence in the file consists only of valid amino acid characters, 
        otherwise False.
    """
    list_aa = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "*"]
    is_aa = True
    dico_seq = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        dico_seq[str(seq_record.id)] = str(seq_record.seq)
    for name in dico_seq:
        my_seq = dico_seq[name]
        for aa in my_seq:
            new_aa = aa.upper()
            if new_aa not in list_aa:
                is_aa = False
        break
    return is_aa
    

def is_gff_strat1(filename : str) -> bool:
    """
    Checks if a file is in GFF (General Feature Format) format.

    This function performs a basic check to see if the file follows the GFF format by verifying 
    that the first non-comment line contains exactly 9 columns. It doesn't perform a full validation 
    of the GFF file, but it checks for basic structure.

    Parameters:
    -----------
    filename : str
        Path to the file to be checked.

    Returns:
    --------
    bool
        True if the file contains at least one valid GFF-like line (with 9 columns), 
        otherwise False.
    """
    my_file = openFile(filename)
    correct_gff = False
    for line in my_file:
        if line[0] != "#":
            elts_file = line.split("\t")
            if len(elts_file) == 9:
                correct_gff = True
            break
    return correct_gff


def is_gtf_strat1(filename : str) -> bool:
    """
    Checks if a file is in GTF (Gene Transfer Format) format.

    This function performs a basic check to see if the file follows the GTF format by verifying 
    that the first non-comment line contains exactly 9 columns, using tab characters as delimiters. 
    It doesn't perform a full validation of the GTF file, but checks for basic structure.

    Parameters:
    -----------
    filename : str
        Path to the file to be checked.

    Returns:
    --------
    bool
        True if the file contains at least one valid GTF-like line (with 9 columns), 
        otherwise False.
    """
    my_file = openFile(filename)
    correct_gtf = False
    for line in my_file:
        if line[0] != "#":
            elts_file = line.split("\t")
            if len(elts_file) == 9:
                correct_gtf = True
            break
    return correct_gtf


def search_liste_name_target_genome_strat1(query_name : str, path_genome_folder : str) -> (list, list, list):
    """
    Extracts the names of target genomes and their annotations from a genome folder.

    This function scans the directory at `path_genome_folder` and extracts the names of genome 
    and annotation files that do not match the specified `query_name`. It searches for files 
    in common genome and annotation formats (FASTA and GFF).

    Parameters:
    -----------
    query_name : str
        The name of the query genome (used to exclude files matching this name).

    path_genome_folder : str
        Path to the genome folder to be scanned for target genome and annotation files.

    Returns:
    --------
    tuple
        A tuple containing three lists:
        - list_target_name: The names of the target genomes (excluding the query genome).
        - liste_target_genome: The list of genome file names in FASTA format (excluding the query genome).
        - liste_target_genome_annotation: The list of genome annotation file names (in GFF format).
    """
    list_fasta_format = ["fa", "fna", "fasta"]
    liste_gff_format = ["gff", "gff3"]
    elts_dir = os.listdir(path_genome_folder)
    liste_target_genome = []
    liste_target_genome_annotation = []
    list_target_name = []
    for name_candidate in elts_dir:
        if (".") in name_candidate:
            header = name_candidate.split(".")[0]
            format = name_candidate.split(".")[1]
            if header != query_name and format in list_fasta_format:
                if header not in list_target_name:
                    list_target_name.append(header)
                liste_target_genome.append(name_candidate)
            elif header != query_name and format in liste_gff_format:
                if header not in list_target_name:
                    list_target_name.append(header)
                liste_target_genome_annotation.append(name_candidate)
    return list_target_name, liste_target_genome, liste_target_genome_annotation


def validate_target_genomes_strat1(dico_variables : dict, RUN_PYTHON : bool) -> (bool, list):
    """
    Validates the presence and format of genome and annotation files for target genomes.

    This function checks whether each target genome, identified by a name in `dico_variables`, 
    has both a corresponding genome file in FASTA format and a corresponding annotation file in GFF format. 
    It also verifies the validity of the file formats for both the genome (FASTA) and the annotation (GFF).

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing various input parameters, including the path to the genome repository 
        (key: "path_to_genome_repository") and the query name (key: "query").

    RUN_PYTHON : bool
        A flag indicating whether the script should continue execution. The flag is set to False if 
        any validation errors occur.

    Returns:
    --------
    tuple
        A tuple containing:
        - RUN_PYTHON (bool): Updated flag indicating whether all validations passed.
        - list_target_name (list): List of target genome names that were validated.
    """
    path_to_genome = dico_variables["path_to_genome_repository"]
    list_target_name, liste_target_genome, liste_target_genome_annotation = search_liste_name_target_genome_strat1(dico_variables["query"], path_to_genome)
    # loop in the list of target names
    for name in list_target_name:
        genome_p = False
        genome_valid = True
        annotation_p = False
        annotation_valid = True
        # check if fasta genomes is present and correct
        for file in liste_target_genome:
            header = file.split(".")[0]
            if name == header:
                genome_p = True
                link_to_genome = path_to_genome + "/" + file
                genome_valid = is_fasta_strat1(link_to_genome)
                break
        # check if gff genomes is present and correct
        for file in liste_target_genome_annotation:
            header = file.split(".")[0]
            if name == header:
                annotation_p = True
                link_to_annotation = path_to_genome + "/" + file
                annotation_valid = is_gff_strat1(link_to_annotation)
                break
        if genome_p == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " has genome annotation file but no genome" + RESET)
        if annotation_p == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " has genome file but no annotation" + RESET)
        if genome_valid == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " genome seems to NOT be in FASTA format" + RESET)
        if annotation_valid == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " genome annotation format (GFF/GFF3) seems to be invalid" + RESET)
            print (" ")
    return RUN_PYTHON, list_target_name


def search_liste_name_target_transcriptome_strat1(query_name : str, path_transcriptome_folder : str) -> (list, list, list):
    """
    Extracts the names of target transcriptomes and their annotations from a transcriptome folder.

    This function scans the directory at `path_transcriptome_folder` and extracts the names of transcriptome 
    and annotation files that do not match the specified `query_name`. It searches for files in common transcriptome 
    and annotation formats (FASTA and GTF).

    Parameters:
    -----------
    query_name : str
        The name of the query transcriptome (used to exclude files matching this name).

    path_transcriptome_folder : str
        Path to the transcriptome folder to be scanned for target transcriptome and annotation files.

    Returns:
    --------
    tuple
        A tuple containing three lists:
        - list_target_name: The names of the target transcriptomes (excluding the query transcriptome).
        - liste_target_transcriptome: The list of transcriptome file names in FASTA format (excluding the query transcriptome).
        - liste_target_transcriptome_annotation: The list of transcriptome annotation file names (in GTF format).
    """
    list_fasta_format = ["fa", "fna", "fasta"]
    liste_gtf_format = ["gtf"]
    elts_dir = os.listdir(path_transcriptome_folder)
    liste_target_transcriptome = []
    liste_target_transcriptome_annotation = []
    list_target_name = []
    for name_candidate in elts_dir:
        if (".") in name_candidate:
            header = name_candidate.split(".")[0]
            format = name_candidate.split(".")[1]
            if header != query_name and format in list_fasta_format:
                if header not in list_target_name:
                    list_target_name.append(header)
                liste_target_transcriptome.append(name_candidate)
            elif header != query_name and format in liste_gtf_format:
                if header not in list_target_name:
                    list_target_name.append(header)
                liste_target_transcriptome_annotation.append(name_candidate)
    return list_target_name, liste_target_transcriptome, liste_target_transcriptome_annotation


def validate_target_transcriptomes_strat1(dico_variables : dict, RUN_PYTHON : bool, list_genome_name : list) -> bool:
    """
    Validates the presence, format, and correspondence of outgroup transcriptomes.

    This function checks the presence of transcriptomes for target genomes and their corresponding annotation files. 
    If a transcriptome is present, it validates the file formats (FASTA for the transcriptome and GTF for the annotation).
    Additionally, the function checks if the transcriptome corresponds to any genome. If the transcriptome does not 
    correspond to any genome, a warning is issued. If any transcriptome or its annotation file is missing or has an invalid format, 
    an error message is raised.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing input variables, including the path to the transcriptome repository 
        (key: "path_to_transcriptome_repository") and the query name (key: "query").

    RUN_PYTHON : bool
        A flag indicating whether the script should continue execution. The flag is set to False if any validation errors occur.

    list_genome_name : list
        A list of genome names that are considered valid for the analysis.

    Returns:
    --------
    bool
        The updated value of `RUN_PYTHON` indicating whether all validations passed (True) or failed (False).
    """
    # extract transcriptomes and check whether one of them does not correspond to any genome
    path_to_transcriptome = dico_variables["path_to_transcriptome_repository"]
    list_transcriptome_name, liste_target_transcriptome, liste_target_transcriptome_annotation = search_liste_name_target_transcriptome_strat1(dico_variables["query"], path_to_transcriptome)
    for name_transcriptome in list_transcriptome_name:
        if name_transcriptome not in list_genome_name:
            print (ORANGE + "WARNING : " + name_transcriptome + " transcriptome does not correspond to a genome" + RESET)
            print (ORANGE + "---> " + name_transcriptome + " is not included in the analysis" + RESET)
            print (" ")
    # iterate through genomes names (transcriptome are required to have the same name that the genome they correspond to)
    for name in list_genome_name:
        name_is_present = False
        transcriptome_p = False
        transcriptome_valid = True
        annotation_p = False
        annotation_valid = True
        # validate presence and fasta format
        for file in liste_target_transcriptome:
            header = file.split(".")[0]
            if name == header:
                name_is_present = True
                transcriptome_p = True
                link_to_transcriptome = path_to_transcriptome + "/" + file
                transcriptome_valid = is_fasta_strat1(link_to_transcriptome)
                break
        # validate presence and annotation format
        for file in liste_target_transcriptome_annotation:
            header = file.split(".")[0]
            if name == header:
                name_is_present = True
                annotation_p = True
                link_to_annotation = path_to_transcriptome + "/" + file
                annotation_valid = is_gtf_strat1(link_to_annotation)
                break
        # write warning messages if needed
        if name_is_present == False:
            print (ORANGE + "WARNING : " + name + " does not have a transcriptome" + RESET)
            print (ORANGE + "---> " + name + " outgroup homologs will not be analysed for transcription" + RESET)
        if transcriptome_p == False and name_is_present == True:
            print (ORANGE + "WARNING : target " + name + " has transcriptome annotation file but no transcriptome" + RESET)
            print (ORANGE + "---> " + name + " outgroup homologs will not be analysed for transcription" + RESET)
        if annotation_p == False and name_is_present == True:
            print (ORANGE + "WARNING : target " + name + " has transcriptome but no corresponding annotation" + RESET)
            print (ORANGE + "---> " + name + " outgroup homologs will not be analysed for transcription" + RESET)
        if transcriptome_valid == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " transcriptome seems to NOT be in FASTA format" + RESET)
        if annotation_valid == False:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : target " + name + " transcriptome annotation format (GTF) seems to be invalid" + RESET)
            print (" ")
    return RUN_PYTHON


def search_name_in_folder_strat1(path_folder : str, name_query : str, list_format : list) -> str:
    """
    Searches for a file with the specified name and format in a given folder.

    This function scans a specified folder (`path_folder`) for a file that matches the `name_query` and has one of 
    the formats specified in the `list_format`. The function returns the full path to the matching file if found, 
    or an empty string if no match is found.

    Parameters:
    -----------
    path_folder : str
        The directory path where the function should search for the file.

    name_query : str
        The name of the file to search for (without extension).

    list_format : list
        A list of acceptable file extensions (formats) to check for. Only files with extensions in this list will be matched.

    Returns:
    --------
    str
        The full path to the file if a match is found; otherwise, an empty string.
    """
    elts_dir = os.listdir(path_folder)
    link_query = ""
    for names in elts_dir:
        if (".") in names:
            header = names.split(".")[0]
            format = names.split(".")[1]
            if header == name_query and format in list_format:
                link_query = path_folder + "/" + names
                break
    return link_query


def validate_file_content_query_strat1(dico_variables : dict, RUN_PYTHON : bool) -> bool:
    """
    Validates the content of query files (genome and transcriptome) to ensure they are in the correct format.

    This function checks if the query files (genome and transcriptome) specified in the `dico_variables` dictionary 
    are present in the corresponding folders and verifies that the content of the files adheres to the correct format 
    for FASTA, GFF, and GTF files. If any of the files are in the wrong format, the function prints an error message 
    and sets `RUN_PYTHON` to `False`.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing paths to the genome and transcriptome folders, as well as the query name under the key "query".
    
    RUN_PYTHON : bool
        A flag indicating whether the script should continue running. It will be set to `False` if any validation fails.

    Returns:
    --------
    bool
        The updated value of `RUN_PYTHON` after performing the validation checks. If any file is in the wrong format, 
        `RUN_PYTHON` will be set to `False`.
    """
    # extract query files
    fasta_genome = search_name_in_folder_strat1(dico_variables["path_to_genome_repository"], dico_variables["query"], ["fasta", "fna", "fa"])
    gff_genome = search_name_in_folder_strat1(dico_variables["path_to_genome_repository"], dico_variables["query"], ["gff", "gff3"])
    fasta_transcriptome = search_name_in_folder_strat1(dico_variables["path_to_transcriptome_repository"], dico_variables["query"], ["fasta", "fna", "fa"])
    gtf_transcriptome = search_name_in_folder_strat1(dico_variables["path_to_transcriptome_repository"], dico_variables["query"], ["gtf"])

    # validate query files
    fasta_genome_correct = is_fasta_strat1(fasta_genome)
    gff_genome_correct = is_gff_strat1(gff_genome)
    fasta_transcriptome_correct = is_fasta_strat1(fasta_transcriptome)
    gtf_transcriptome_correct = is_gtf_strat1(gtf_transcriptome)

    # error messages if needed
    if fasta_genome_correct == False:
        RUN_PYTHON = False
        print (BRIGHT_RED + "ERROR : The query genome does not seem to be in FASTA format" + RESET)

    if gff_genome_correct == False:
        RUN_PYTHON = False
        print (BRIGHT_RED + "ERROR : The query genome annotation does not seem to be in correct GFF format" + RESET)

    if fasta_transcriptome_correct == False:
        RUN_PYTHON = False
        print (BRIGHT_RED + "ERROR : The query transcriptome does not seem to be in FASTA format" + RESET)

    if gtf_transcriptome_correct == False:
        RUN_PYTHON = False
        print (BRIGHT_RED + "ERROR : The query transcriptome annotation does not seem to be in correct GTF format" + RESET)
    return RUN_PYTHON


def display_parameters_strat1(dico_variables : dict) -> None:
    """
    Display all parameters entered by the user for the strategy.

    This function prints out the values of various parameters that have been provided by the user 
    in the `dico_variables` dictionary. These parameters include information about the strategy, 
    query, file paths, filtering options, and any databases that are linked for validation of homology.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing all the parameters and values entered by the user. The dictionary is expected 
        to contain keys such as "strategy", "query", "path_to_genome_repository", "path_to_transcriptome_repository", 
        "TPM_threeshold", "transcript_overlap", "ORFs_choice", "filter_genic", "filter_TE", "rm_undir_transc", 
        "rec_best_hit", "synteny_window", "premature_stop", and keys for linked databases ("link_database_outgroup_prot" 
        and "link_database_outgroup_nucl").

    Returns:
    --------
    None
        The function does not return any value. It simply prints out the parameter values to the console.
    
    Notes:
    ------
    - If "ORFs_choice" contains options for UTR size, the corresponding sizes for 5'UTR and 3'UTR will be printed 
      if they are non-zero.
    - Warnings will be displayed if no protein or DNA homology databases are provided for validation.
    - The function makes use of constants `BRIGHT_GREEN` and `RESET` for formatting the output (likely used for color output).
    """
    print (BRIGHT_GREEN + "Strategy : " + RESET + str(dico_variables["strategy"]))
    print (BRIGHT_GREEN + "Query name : " + RESET + str(dico_variables["query"]))
    print (BRIGHT_GREEN + "Genomes folder : " + RESET + str(dico_variables["path_to_genome_repository"]))
    print (BRIGHT_GREEN + "Transcriptome(s) folder : " + RESET + str(dico_variables["path_to_transcriptome_repository"]))
    print (BRIGHT_GREEN + "Minimal transcripts TPM threshold : " + RESET + str(dico_variables["TPM_threeshold"]))
    print (BRIGHT_GREEN + "Genomic position of candidate : " + RESET + str(dico_variables["transcript_overlap"]))
    for option in dico_variables["ORFs_choice"]:
        if option[0] == "utr_size":
            if option[1] != 0:
                print (BRIGHT_GREEN + "Minimum size of 5'UTR (nb nucl): " + RESET + str(option[1]))
            if option[2] != 0:
                print (BRIGHT_GREEN + "Minimum size of 3'UTR (nb nucl): " + RESET + str(option[2]))
        elif option[0] != "utr_size" and option[0] != "duplicate_handle":
            print (BRIGHT_GREEN + "ORFs in transcripts : " + RESET + str(option[0]))
    print (BRIGHT_GREEN + "Remove spliced events : " + RESET + str(dico_variables["filter_genic"]))
    print (BRIGHT_GREEN + "Filter TEs : " + RESET + str(dico_variables["filter_TE"]))
    print (BRIGHT_GREEN + "Remove unoriented transcripts : " + RESET + str(dico_variables["rm_undir_transc"]))
    print (BRIGHT_GREEN + "Reciprocal BEST hits : " + RESET + str(dico_variables["rec_best_hit"]))
    print (BRIGHT_GREEN + "Synteny window : " + RESET + str(dico_variables["synteny_window"]))
    print (BRIGHT_GREEN + "Search premature stop in (sequence percentage) : " + RESET + str(dico_variables["premature_stop"]) + " %")
    print (" ")
    if dico_variables["link_database_outgroup_prot"] != "":
        print (BRIGHT_GREEN + "Validation of NO PROTEIN homology to neORF will be performed against : " + RESET)
        print (str(dico_variables["link_database_outgroup_prot"]))
        print (BRIGHT_GREEN + "With mode : " + RESET + str(dico_variables["parameters_database_prot"]["mode"]))
    else:
        print (ORANGE + "WARNING : Validation of NO PROTEIN homology NOT performed (no database provided)" + RESET)

    print (" ")

    if dico_variables["link_database_outgroup_nucl"] != "":
        print (BRIGHT_GREEN + "Validation of NO DNA homology to neORF will be performed against : " + RESET)
        print (str(dico_variables["link_database_outgroup_nucl"]))
        print (BRIGHT_GREEN + "With e-value : " + RESET + str(dico_variables["parameters_database_nucl"]["e_value"]))
    else:
        print (ORANGE + "WARNING : Validation of NO DNA homology NOT performed (no database provided)" + RESET)

    print (" ")


def assess_blast_datasets_strat1(dico_variables : dict, RUN_PYTHON : bool) -> bool:
    """
    Assess the format and content of the datasets for BLAST search against outgroup (step 2).

    This function validates the input datasets for protein and nucleotide/rna homology searches, 
    as indicated by the user in the `dico_variables`. It checks if the datasets are in the correct 
    format (FASTA) and whether they contain the expected content (proteins for protein datasets, 
    nucleotides or RNA for nucleotide datasets). It also prints error messages if the datasets are 
    not correctly formatted or do not contain the expected data.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing the paths to the protein and nucleotide/rna homology search datasets. 
        Expected keys include:
        - "link_database_outgroup_prot" : Path to the protein dataset (string).
        - "link_database_outgroup_nucl" : Path to the nucleotide/rna dataset (string).
    
    RUN_PYTHON : bool
        A flag indicating whether the program should continue running. If any validation fails, 
        this flag is set to `False`.

    Returns:
    --------
    bool
        The updated value of `RUN_PYTHON`, which will be `False` if any dataset is in an invalid 
        format or does not contain the expected data.

    Notes:
    ------
    - The function first checks the protein dataset (if provided) for being in FASTA format and 
      containing proteins. If both checks pass, it continues.
    - If the nucleotide/rna dataset is provided, the function checks for FASTA format and verifies 
      that the dataset contains DNA or RNA sequences.
    - If any of the checks fail, an error message is printed and `RUN_PYTHON` is set to `False`.
    """
    link_prot_dataset = dico_variables["link_database_outgroup_prot"]
    link_nucl_dataset = dico_variables["link_database_outgroup_nucl"]
    # run only for proteins if a protein file was indicated
    if link_prot_dataset != "":
        # assess whether protein file is in FASTA
        RUN_PYTHON = is_fasta_strat1(link_prot_dataset)
        if RUN_PYTHON == False:
            print (BRIGHT_RED + "ERROR : dataset for protein homology search seems to not be in FASTA format" + RESET)
        else:
            # Assess whether protein file contains proteins
            RUN_PYTHON = is_prot_strat1(link_prot_dataset)
            if RUN_PYTHON == False:
                print (BRIGHT_RED + "ERROR : dataset for protein homology search does not seem to contain proteins" + RESET)
    # run only for nucleotides/rna if nucleotide/rna file was indicated
    if RUN_PYTHON == True and link_nucl_dataset != "":
        # assess whether dna/rna file is in FASTA
        RUN_PYTHON = is_fasta_strat1(link_nucl_dataset)
        if RUN_PYTHON == False:
            print (BRIGHT_RED + "ERROR : dataset for nucl homology search seems to not be in FASTA format" + RESET)
        else:
            # Assess whether dna/rna file contain dna/rna
            RUN_PYTHON = is_nucl_strat1(link_nucl_dataset)
            if RUN_PYTHON == False:
                print (BRIGHT_RED + "ERROR : dataset for nucl homology search does not seem to contain DNA/RNA" + RESET)
    return RUN_PYTHON


def display_welcome() -> None:
    ## display DESwoMAN welcome message ##
    print (BRIGHT_GREEN + " " + RESET)
    print (YELLOW + "Welcome to DESwoMAN !!" + RESET)
    print (BRIGHT_GREEN + " " + RESET)

    desman = r"""
        🌸   
    (\_ _/)
    ( o.o )  
     (>o<)
   (   "   )    
    |   
    |   
    |
    """
    print(desman)

    print (BRIGHT_GREEN + "- - - - - - - - - - - - -" + RESET)
    print (BRIGHT_GREEN + "| PARAMETERS ASSESSMENT |" + RESET)
    print (BRIGHT_GREEN + "- - - - - - - - - - - - -" + RESET)
    print (BRIGHT_GREEN + " " + RESET)


def assess_parameters_strat1(link_config : str) -> (bool, dict):
    """
    Main function to assess the user-provided parameters and file contents in the DESwoMAN pipeline.

    This function is responsible for performing a series of checks on the parameters provided by the user
    through the graphical interface. It validates the presence of required folders and query files (genomes, 
    transcriptomes), ensures that the formats of these files are correct, and confirms that all necessary datasets 
    (for BLAST searches) are available and properly formatted. It also checks that at least one target genome is 
    present and that the genomic overlap settings are valid.

    The function ensures that all necessary parameters and files are ready for the next steps of the DESwoMAN pipeline. 
    If any validation step fails, the function will stop execution and provide an error message.

    Parameters:
    -----------
    link_config (bool|str): The link to the config file, or FALSE which indicates DESwoMAN runs with the graphical interface

    Returns:
    --------
    RUN_PYTHON : bool
        A flag indicating whether the program should continue (`True` if all checks pass, `False` otherwise).
    
    dico_variables : dict
        The updated dictionary of variables containing the user-provided parameters, including paths to input data 
        and settings for the pipeline. If the function completes successfully, additional output paths are also included.
    """



    RUN_PYTHON = True
    if link_config == False:
        dico_variables = my_graphical_interface_strategy1() # dico_variables is retrived from the parameters chosen by the user in the graphical interface
        display_welcome()
    else:
        display_welcome()
        dico_variables, RUN_PYTHON = my_config_file_extract_parameters(link_config, 1)
    if RUN_PYTHON == True:
        # display the user choices entered in the graphical interface so that the person can realise if some parameter has to be modified
        display_parameters_strat1(dico_variables) 
        # assess the presence of mandatory folders path and query name
        RUN_PYTHON = validate_presence_of_mandatory_parameters_strat1(dico_variables, RUN_PYTHON) 
    if RUN_PYTHON == True:
        # assess that query genome and transcriptomes are in the folders
        RUN_PYTHON = validate_presence_query_name_strat1(dico_variables, RUN_PYTHON) 
    if RUN_PYTHON == True:
        # assess that the user did not unselect all denovo genomic locations (intergenic, intronic, genic, antisense)
        RUN_PYTHON = check_genomic_overlap_strat1(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # try to make sure all query files (fasta, gff, gtf) are not corrupted or in a wrong format
        RUN_PYTHON = validate_file_content_query_strat1(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # verify that all target have the genomes files (fasta and gff) and that the format are correct
        RUN_PYTHON, list_target_name = validate_target_genomes_strat1(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # make sure there is at least 1 target genome
        if len(list_target_name) == 0:
            RUN_PYTHON = False
            print (BRIGHT_RED + "ERROR : The genome folder contains 0 outgroup genomes (min required : 1)" + RESET)
    if RUN_PYTHON == True:
        # make sure the dataset for BLAST contain DNA and proteins
        RUN_PYTHON = assess_blast_datasets_strat1(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # verify that all target have the transcriptome files (fasta and gff) and that the format are correct
        RUN_PYTHON = validate_target_transcriptomes_strat1(dico_variables, RUN_PYTHON, list_target_name)
        # implement the arbitraty path to the folders where to store results
        dico_variables["path_output"] = dico_variables["path_to_genome_repository"] + "/DESwoMAN_denovo_output"
        dico_variables["path_output_intermediate"] = dico_variables["path_output"] + "/Intermediate_output"
    return RUN_PYTHON, dico_variables



##### test ######

#dico_variables = my_graphical_interface_strategy1()
#dico_variables = {'strategy': '1', 
                  #'query': 'AK5', 
                  #'path_to_genome_repository': '/home/anna/Bureau/Allemagne_recherche/DESMAN/desman_feb16_with_strat_2/genomes_strat1', 
                  #'path_to_transcriptome_repository': '/home/anna/Bureau/Allemagne_recherche/DESMAN/desman_feb16_with_strat_2/transcriptomes_strat1', 
                  #'link_database_outgroup_prot': '/home/anna/Bureau/Allemagne_recherche/DESMAN/blast_database/All_Ants.fasta', 
                  #'link_database_outgroup_nucl': '', 
                  #'TPM_threeshold': 0.5, 
                  #'transcript_overlap': ['intergenic'], 
                  #'ORFs_choice': [['longest'], ['duplicate_handle'], ['utr_size', '10', '10']], 
                  #'filter_genic': False, 
                  #'filter_TE': 'False', 
                  #'rm_undir_transc': 'False', 
                  #'parameters_database_prot': {'type': 'blastp', 'mode': '--more-sensitive'}, 
                  #'parameters_database_nucl': {'type': 'blastn', 'e_value': '0.01', 'coverage': None, 'strand': None}, 
                  #'rec_best_hit': 'False', 
                  #'synteny_window': 2, 
                  #'premature_stop': 50}
#print (dico_variables)
#/home/anna/Bureau/Allemagne_recherche/DESMAN/blast_database/AK5-families.fa
