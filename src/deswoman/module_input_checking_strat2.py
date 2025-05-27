import os
from Bio import SeqIO
from deswoman.module_graphical_interface_strat2 import my_graphical_interface_strategy2
from deswoman.module_handle_config_file import my_config_file_extract_parameters
from deswoman.module_colors import *


def validate_presence_of_mandatory_parameters_strat2(
    dico_variables: dict, RUN_PYTHON: bool
) -> bool:
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
        print(BRIGHT_RED + "ERROR... NO QUERY NAME PROVIDED ... " + RESET)
        RUN_PYTHON = False
    if dico_variables["path_to_genome_repository"] == "":
        print(BRIGHT_RED + "ERROR... NO GENOME FOLDER PROVIDED ... " + RESET)
        RUN_PYTHON = False
    if dico_variables["path_to_transcriptome_repository"] == "":
        print(BRIGHT_RED + "ERROR... NO TRANSCRIPTOME FOLDER PROVIDED ... " + RESET)
        RUN_PYTHON = False
    return RUN_PYTHON


def validate_presence_query_genome_strat2(query_name: str, path_rep: str) -> bool:
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


def validate_presence_query_gff_strat2(query_name: str, path_rep: str) -> bool:
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


def validate_presence_query_transcriptome_strat2(
    query_name: str, path_rep: str
) -> bool:
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


def validate_presence_query_gtf_strat2(query_name: str, path_rep: str) -> bool:
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


def validate_presence_query_name_strat2(dico_variables: dict, RUN_PYTHON: bool) -> bool:
    """
    Validate the presence of the query genome files in the genome folder.

    This function checks if the genome files corresponding to the query genome (in FASTA and GFF formats)
    are present in the specified genome folder. It ensures that both the genome sequence file (FASTA format)
    and the genome annotation file (GFF format) exist before proceeding further with the analysis.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing the user-provided parameters. It should include:
        - "query": the name of the query genome to be searched.
        - "path_to_genome_repository": the directory path where genome files are stored.

    RUN_PYTHON : bool
        A flag indicating whether the program should continue (`True` if all checks pass, `False` if any check fails).

    Returns:
    --------
    RUN_PYTHON : bool
        A flag indicating whether the program should continue after validating the presence of the query genome files.
        If any of the required files are missing, `RUN_PYTHON` is set to `False`, and error messages are printed.
    """
    query_name = dico_variables["query"]
    path_to_genome = dico_variables["path_to_genome_repository"]
    # path_to_transcriptome = dico_variables["path_to_transcriptome_repository"]
    presence_fasta_genome = validate_presence_query_genome_strat2(
        query_name, path_to_genome
    )
    presence_gff = validate_presence_query_gff_strat2(query_name, path_to_genome)
    # make sure a reference genome is present in fasta and gff format
    if presence_fasta_genome == False:
        RUN_PYTHON = False
        name_fasta_genome = query_name + (".fasta")
        print(
            BRIGHT_RED
            + "ERROR : Query genome "
            + name_fasta_genome
            + " (or .fa|.fna) not present in genome folder"
            + RESET
        )
    if presence_gff == False:
        RUN_PYTHON = False
        name_gff_genome = query_name + (".gff")
        print(
            BRIGHT_RED
            + "ERROR : Query genome annotation "
            + name_gff_genome
            + " (or .gff3) not present in genome folder"
            + RESET
        )
    return RUN_PYTHON


def check_genomic_overlap_strat2(dico_variables: dict, RUN_PYTHON: bool) -> bool:
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
        print(BRIGHT_RED + "ERROR : No transcript location provided" + RESET)
    return RUN_PYTHON


def is_fasta_strat2(filename: str) -> bool:
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


def is_nucl_strat2(filename: str) -> bool:
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
    list_nucl = [
        "A",
        "C",
        "G",
        "T",
        "U",
        "R",
        "Y",
        "K",
        "M",
        "S",
        "W",
        "B",
        "D",
        "H",
        "V",
        "N",
    ]
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


def is_prot_strat2(filename: str) -> bool:
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
    list_aa = [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z",
        "*",
    ]
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


def is_gff_strat2(filename: str) -> bool:
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


def is_gtf_strat2(filename: str) -> bool:
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
            elts_file = line.split("	")
            if len(elts_file) == 9:
                correct_gtf = True
            break
    return correct_gtf


def search_liste_name_target_transcriptome_strat2(
    query_name: str, path_transcriptome_folder: str
) -> (list, list, list):
    """
    Extracts the names of target transcriptomes from a transcriptome folder.

    This function scans the given folder containing transcriptome files and identifies
    transcriptome files that belong to targets (excluding the query genome). It organizes
    the files into two categories: transcriptome files (e.g., FASTA format) and their corresponding
    annotation files (e.g., GTF format). The function returns the list of target names, transcriptome files,
    and annotation files for further processing.

    Parameters:
    -----------
    query_name : str
        The name of the query genome. This genome will be excluded from the search.

    path_transcriptome_folder : str
        The directory path where the transcriptome files are located.

    Returns:
    --------
    list_target_name : list
        A list of names of target genomes (excluding the query genome) that have corresponding transcriptome files.

    liste_target_transcriptome : list
        A list of names of transcriptome files in FASTA format corresponding to the target genomes.

    liste_target_transcriptome_annotation : list
        A list of names of transcriptome annotation files in GTF format corresponding to the target genomes.
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
    return (
        list_target_name,
        liste_target_transcriptome,
        liste_target_transcriptome_annotation,
    )


def validate_target_transcriptomes_strat2(
    dico_variables: dict, RUN_PYTHON: bool
) -> (bool, list):
    """
    Assesses the presence and validity of target transcriptomes and their corresponding GTF annotations.

    This function verifies that each target genome has both a corresponding transcriptome file
    (in FASTA format) and an annotation file (in GTF format). It checks the format validity of
    these files and removes any target genomes that lack either file or have invalid formats.
    Warnings are issued for missing files, and errors are raised for invalid formats.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing parameters and file paths, including the query genome name and the path
        to the transcriptome repository.

    RUN_PYTHON : bool
        A flag indicating whether to continue the process. It will be set to False if errors are encountered
        with the format or presence of transcriptome and annotation files.

    Returns:
    --------
    RUN_PYTHON : bool
        The updated flag indicating whether the process should continue after validation.

    list_transcriptome_name : list
        A list of target transcriptome names (genomes) that passed the validation for the presence and format
        of transcriptome and annotation files.
    """
    path_to_transcriptome = dico_variables["path_to_transcriptome_repository"]
    (
        list_transcriptome_name,
        liste_target_transcriptome,
        liste_target_transcriptome_annotation,
    ) = search_liste_name_target_transcriptome_strat2(
        dico_variables["query"], path_to_transcriptome
    )
    liste_name_to_remove = []
    # iterate in query transcriptomes
    for name in list_transcriptome_name:
        transcriptome_p = False
        transcriptome_valid = True
        annotation_p = False
        annotation_valid = True
        # validate FASTA transcriptome
        for file in liste_target_transcriptome:
            header = file.split(".")[0]
            if name == header:
                transcriptome_p = True
                link_to_transcriptome = path_to_transcriptome + "/" + file
                transcriptome_valid = is_fasta_strat2(link_to_transcriptome)
                break
        # validate gtf transcriptome
        for file in liste_target_transcriptome_annotation:
            header = file.split(".")[0]
            if name == header:
                annotation_p = True
                link_to_annotation = path_to_transcriptome + "/" + file
                annotation_valid = is_gtf_strat2(link_to_annotation)
                break
        # write warning messages if needed
        if transcriptome_p == False:
            print(
                ORANGE
                + "WARNING : query "
                + name
                + " has transcriptome annotation file but no transcriptome"
                + RESET
            )
            print(
                ORANGE
                + "---> "
                + name
                + " will not be analysed for transcription"
                + RESET
            )
            liste_name_to_remove.append(name)
        if annotation_p == False:
            print(
                ORANGE
                + "WARNING : query "
                + name
                + " has transcriptome file but no annotation"
                + RESET
            )
            print(
                ORANGE
                + "---> "
                + name
                + " will not be analysed for transcription"
                + RESET
            )
            liste_name_to_remove.append(name)
        if transcriptome_valid == False:
            RUN_PYTHON = False
            print(
                BRIGHT_RED
                + "ERROR : target "
                + name
                + " transcriptome seems to NOT be in FASTA format"
                + RESET
            )
        if annotation_valid == False:
            RUN_PYTHON = False
            print(
                BRIGHT_RED
                + "ERROR : target "
                + name
                + " transcriptome annotation format (GTF) seems to be invalid"
                + RESET
            )
            print(" ")
    for i in liste_name_to_remove:
        list_transcriptome_name.remove(i)
    return RUN_PYTHON, list_transcriptome_name


def search_name_in_folder_strat2(
    path_folder: str, name_query: str, list_format: list
) -> str:
    """
    Search for a file with a specific name and format in a given folder.

    This function looks for a file in the specified folder (`path_folder`) that has a name matching the
    provided query (`name_query`) and is of one of the acceptable formats listed in `list_format`.
    It returns the full path to the file if found, otherwise returns an empty string.

    Parameters:
    -----------
    path_folder : str
        The path to the folder where the search will be performed.

    name_query : str
        The query name (without the file extension) that the function will search for in the folder.

    list_format : list
        A list of acceptable file formats (extensions) to look for. Only files with a format from this list
        will be considered.

    Returns:
    --------
    link_query : str
        The full path to the file if a matching file is found, otherwise an empty string.
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


def validate_file_content_query_strat2(dico_variables: dict, RUN_PYTHON: bool) -> bool:
    """
    Assess whether the query genome files are in the correct format.

    This function validates the content of the query genome and its associated genome annotation file.
    It checks if the genome file is in the correct FASTA format and if the genome annotation file is
    in a valid GFF format. If any of the files are in an incorrect format, the function will print
    error messages and set `RUN_PYTHON` to `False` to prevent further execution.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing configuration variables such as paths to the genome repository and the query genome name.

    RUN_PYTHON : bool
        A flag indicating whether the process should continue. If any files are in an incorrect format, this will be set to `False`.

    Returns:
    --------
    RUN_PYTHON : bool
        Returns the `RUN_PYTHON` flag, which is set to `False` if any file is in an invalid format; otherwise, it remains `True`.
    """
    fasta_genome = search_name_in_folder_strat2(
        dico_variables["path_to_genome_repository"],
        dico_variables["query"],
        ["fasta", "fna", "fa"],
    )
    gff_genome = search_name_in_folder_strat2(
        dico_variables["path_to_genome_repository"],
        dico_variables["query"],
        ["gff", "gff3"],
    )
    fasta_genome_correct = is_fasta_strat2(fasta_genome)
    gff_genome_correct = is_gff_strat2(gff_genome)

    if fasta_genome_correct == False:
        RUN_PYTHON = False
        print(
            BRIGHT_RED
            + "ERROR : The query genome does not seem to be in FASTA format"
            + RESET
        )

    if gff_genome_correct == False:
        RUN_PYTHON = False
        print(
            BRIGHT_RED
            + "ERROR : The query genome annotation does not seem to be in correct GFF format"
            + RESET
        )
    return RUN_PYTHON


def display_parameters_strat2(dico_variables: dict) -> None:
    """
    Display all parameters entered by the user for Strategy 2.

    This function prints out the key parameters that the user has entered or selected in the graphical interface
    for Strategy 2. These parameters include details about the strategy type, query genome, transcriptome,
    minimum thresholds, genomic positions, and various choices related to Open Reading Frames (ORFs), filters,
    and database validation. It provides a summary of the user's configuration to ensure that the settings are
    correctly chosen and to help the user identify any necessary modifications.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing the user-defined parameters for Strategy 2. The dictionary should include:
            - "strategy": The selected strategy (e.g., "strategy 2").
            - "query": The name of the query genome.
            - "path_to_genome_repository": Path to the genome repository.
            - "path_to_transcriptome_repository": Path to the transcriptome repository.
            - "TPM_threeshold": Minimal transcripts TPM threshold.
            - "transcript_overlap": Genomic position of candidate.
            - "ORFs_choice": A list of choices related to ORFs in the transcripts.
            - "filter_genic": Boolean for removing spliced events.
            - "filter_TE": Boolean for filtering TEs (transposable elements).
            - "rm_undir_transc": Boolean for removing unoriented transcripts.
            - "link_database_outgroup_prot": Path to a protein database for homology validation.
            - "link_database_outgroup_nucl": Path to a nucleotide database for homology validation.
            - "parameters_database_prot": Parameters for protein database (e.g., mode).
            - "parameters_database_nucl": Parameters for nucleotide database (e.g., e-value).
    """
    print(BRIGHT_GREEN + "Strategy : " + RESET + str(dico_variables["strategy"]))
    print(BRIGHT_GREEN + "Query name : " + RESET + str(dico_variables["query"]))
    print(
        BRIGHT_GREEN
        + "Genomes folder : "
        + RESET
        + str(dico_variables["path_to_genome_repository"])
    )
    print(
        BRIGHT_GREEN
        + "Transcriptome(s) folder : "
        + RESET
        + str(dico_variables["path_to_transcriptome_repository"])
    )
    print(
        BRIGHT_GREEN
        + "Minimal transcripts TPM threshold : "
        + RESET
        + str(dico_variables["TPM_threeshold"])
    )
    print(
        BRIGHT_GREEN
        + "Genomic position of candidate : "
        + RESET
        + str(dico_variables["transcript_overlap"])
    )
    for option in dico_variables["ORFs_choice"]:
        if option[0] == "utr_size":
            if option[1] != 0:
                print(
                    BRIGHT_GREEN
                    + "Minimum size of 5'UTR (nb nucl): "
                    + RESET
                    + str(option[1])
                )
            if option[2] != 0:
                print(
                    BRIGHT_GREEN
                    + "Minimum size of 3'UTR (nb nucl): "
                    + RESET
                    + str(option[2])
                )
        elif option[0] != "utr_size" and option[0] != "duplicate_handle":
            print(BRIGHT_GREEN + "ORFs in transcripts : " + RESET + str(option[0]))
    print(
        BRIGHT_GREEN
        + "Remove spliced events : "
        + RESET
        + str(dico_variables["filter_genic"])
    )
    print(BRIGHT_GREEN + "Filter TEs : " + RESET + str(dico_variables["filter_TE"]))
    print(
        BRIGHT_GREEN
        + "Remove unoriented transcripts :  : "
        + RESET
        + str(dico_variables["rm_undir_transc"])
    )
    print(" ")
    if dico_variables["link_database_outgroup_prot"] != "":
        print(
            BRIGHT_GREEN
            + "Validation of NO PROTEIN homology will be performed against : "
            + RESET
        )
        print(str(dico_variables["link_database_outgroup_prot"]))
        print(
            BRIGHT_GREEN
            + "With mode : "
            + RESET
            + str(dico_variables["parameters_database_prot"]["mode"])
        )
    else:
        print(
            ORANGE
            + "WARNING : Validation of NO PROTEIN homology NOT performed (no database provided)"
            + RESET
        )

    print(" ")

    if dico_variables["link_database_outgroup_nucl"] != "":
        print(
            BRIGHT_GREEN
            + "Validation of NO DNA homology will be performed against : "
            + RESET
        )
        print(str(dico_variables["link_database_outgroup_nucl"]))
        print(
            BRIGHT_GREEN
            + "With e-value : "
            + RESET
            + str(dico_variables["parameters_database_nucl"]["e_value"])
        )
    else:
        print(
            ORANGE
            + "WARNING : Validation of NO DNA homology NOT performed (no database provided)"
            + RESET
        )

    print(" ")


def assess_blast_datasets_strat2(dico_variables: dict, RUN_PYTHON: bool) -> bool:
    """
    Assess the format and content of the datasets for BLAST search against outgroup in Strategy 2.

    This function evaluates whether the provided datasets for protein and nucleotide homology searches
    (for outgroup comparison) are correctly formatted in FASTA format, and whether they contain valid protein
    or nucleotide sequences. If any dataset fails these checks, the function will print an error message and
    set `RUN_PYTHON` to False, signaling that the process should be halted.

    Parameters:
    -----------
    dico_variables : dict
        A dictionary containing user-defined variables, including the paths to the protein and nucleotide
        datasets (for homology searches).

        - "link_database_outgroup_prot" : Path to the protein homology database.
        - "link_database_outgroup_nucl" : Path to the nucleotide homology database.

    RUN_PYTHON : bool
        A boolean that controls whether the program should continue running. If any dataset fails the checks,
        this value is set to False.

    Returns:
    --------
    bool
        The updated `RUN_PYTHON` status, which will be False if any dataset fails the format/content validation.
    """
    link_prot_dataset = dico_variables["link_database_outgroup_prot"]
    link_nucl_dataset = dico_variables["link_database_outgroup_nucl"]

    if link_prot_dataset != "":
        RUN_PYTHON = is_fasta_strat2(link_prot_dataset)
        if RUN_PYTHON == False:
            print(
                BRIGHT_RED
                + "ERROR : dataset for protein homology search seems to not be in FASTA format"
                + RESET
            )
        else:
            RUN_PYTHON = is_prot_strat2(link_prot_dataset)
            if RUN_PYTHON == False:
                print(
                    BRIGHT_RED
                    + "ERROR : dataset for protein homology search does not seem to contain proteins"
                    + RESET
                )
    if RUN_PYTHON == True and link_nucl_dataset != "":
        RUN_PYTHON = is_fasta_strat2(link_nucl_dataset)
        if RUN_PYTHON == False:
            print(
                BRIGHT_RED
                + "ERROR : dataset for nucl homology search seems to not be in FASTA format"
                + RESET
            )
        else:
            RUN_PYTHON = is_nucl_strat2(link_nucl_dataset)
            if RUN_PYTHON == False:
                print(
                    BRIGHT_RED
                    + "ERROR : dataset for nucl homology search does not seem to contain DNA/RNA"
                    + RESET
                )
    return RUN_PYTHON


def display_welcome() -> None:
    ## display DESwoMAN welcome message ##
    print(BRIGHT_GREEN + " " + RESET)
    print(YELLOW + "Welcome to DESwoMAN !!" + RESET)
    print(BRIGHT_GREEN + " " + RESET)

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

    print(BRIGHT_GREEN + "- - - - - - - - - - - - -" + RESET)
    print(BRIGHT_GREEN + "| PARAMETERS ASSESSMENT |" + RESET)
    print(BRIGHT_GREEN + "- - - - - - - - - - - - -" + RESET)
    print(BRIGHT_GREEN + " " + RESET)


def assess_parameters_strat2(link_config: str) -> (bool, dict):
    """
    Main function for Strategy 2 parameter assessment in the DESwoMAN workflow.

    This function handles the user input parameters, checks the necessary folders and files,
    and validates the dataset formats required for Strategy 2. The purpose is to ensure that
    all necessary data and configurations are present and correct before proceeding with the analysis.

    It performs the following checks:
    1. Validates the presence and paths of required folders (e.g., genome and transcriptome folders).
    2. Ensures the query genome and transcriptomes are available in the appropriate folders.
    3. Verifies that at least two valid transcriptomes are present.
    4. Confirms that the query genome files (FASTA and GFF) are correctly formatted and not corrupted.
    5. Checks for the validity and presence of BLAST datasets for protein and nucleotide homology searches.
    6. Generates output directories for storing results.

    Parameters:
    -----------
    link_config (bool|str): The link to the config file, or FALSE which indicates DESwoMAN runs with the graphical interface

    Returns:
    --------
    bool
        A boolean value `RUN_PYTHON` indicating whether the checks passed (`True`) or failed (`False`).

    dict
        The updated `dico_variables` dictionary containing user-provided parameters and paths for the workflow.
    """

    RUN_PYTHON = True
    if link_config == False:
        dico_variables = my_graphical_interface_strategy2()  # dico_variables is retrived from the parameters chosen by the user in the graphical interface
        display_welcome()
    else:
        display_welcome()
        dico_variables, RUN_PYTHON = my_config_file_extract_parameters(link_config, 2)
    if RUN_PYTHON == True:
        # display the user choices entered in the graphical interface so that the person can realise if some parameter has to be modified
        display_parameters_strat2(dico_variables)
        # assess the presence of mandatory folders path and query name
        RUN_PYTHON = validate_presence_of_mandatory_parameters_strat2(
            dico_variables, RUN_PYTHON
        )
    if RUN_PYTHON == True:
        # assess that query genome is in the folders
        RUN_PYTHON = validate_presence_query_name_strat2(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # assess that the user did not unselect all transcripts genomic locations
        RUN_PYTHON = check_genomic_overlap_strat2(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # try to make sure all query files (fasta, gff) are not corrupted or in a wrong format
        RUN_PYTHON = validate_file_content_query_strat2(dico_variables, RUN_PYTHON)
    if RUN_PYTHON == True:
        # withdraw query transcriptomes names and validate formats for transcriptoems
        RUN_PYTHON, list_transcriptome_name = validate_target_transcriptomes_strat2(
            dico_variables, RUN_PYTHON
        )
    if RUN_PYTHON == True:
        # a minimum of 2 query transcriptomes are required for strategy 2
        if len(list_transcriptome_name) < 2:
            RUN_PYTHON = False
            print(
                BRIGHT_RED
                + "ERROR : The transcriptome folder contains 0 or 1 valid transcriptome (min required : 2)"
                + RESET
            )
    if RUN_PYTHON == True:
        # make sure the dataset for BLAST contain DNA and proteins
        RUN_PYTHON = assess_blast_datasets_strat2(dico_variables, RUN_PYTHON)
        # implement the arbitraty path to the folders where to store results
        dico_variables["path_output"] = (
            dico_variables["path_to_genome_repository"] + "/DESwoMAN_denovo_output"
        )
        dico_variables["path_output_intermediate"] = (
            dico_variables["path_output"] + "/Intermediate_output"
        )
    return RUN_PYTHON, dico_variables


### test ###

# dico_variables = my_graphical_interface_strategy1()
# dico_variables = {'strategy': '2',
#'query': 'Ref',
#'path_to_genome_repository': '/home/anna/Bureau/Allemagne_recherche/DESMAN/desman_feb16_with_strat_2/genomes_strat2',
#'path_to_transcriptome_repository': '/home/anna/Bureau/Allemagne_recherche/DESMAN/desman_feb16_with_strat_2/transcriptomes_strat2',
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
# print (dico_variables)
# /home/anna/Bureau/Allemagne_recherche/DESMAN/blast_database/AK5-families.fa
