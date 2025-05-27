import os
from deswoman.module_colors import openFile


def perform_blast_nucl_homology_filter(
    name_output_directory: str,
    name_intermediate_directory: str,
    path_to_db: str,
    parameters_db: dict,
) -> None:
    """
    Performs a nucleotide BLAST search against a target DNA/RNA dataset for query denovo candidates.

    This function constructs and runs the necessary BLAST commands to:
    1. Create a BLAST nucleotide database from the provided dataset.
    2. Execute a BLAST search (blastn) of query denovo candidates against the target nucleotide database.
    3. Store the output of the BLAST search in a specified output file.
    4. Clean up by removing temporary BLAST database files.

    Args:
        name_output_directory (str): Path to the directory containing denovo nucleotide sequences.
        name_intermediate_directory (str): Path to the intermediate directory where output files will be stored.
        path_to_db (str): Path to the nucleotide sequence database to be used for the BLAST search.
        parameters_db (dict): A dictionary containing parameters for the BLAST search, specifically the e-value ('e_value').

    Returns:
        None: This function does not return anything but generates output files in the intermediate directory.

    Side Effects:
        - Creates a BLAST database using the `makeblastdb` command.
        - Executes a BLAST search using `blastn`.
        - Generates output files (e.g., `denovo_blast_output_nucl.txt`, `blast_n_makebd.txt`, `blast_n_output.txt`).
        - Removes temporary BLAST database files after the search.
    """
    # Constructing the command to create a BLAST database
    command1 = "makeblastdb -in " + path_to_db + " -dbtype nucl"

    # Executing the command to create the BLAST database and redirecting output to a file
    redirect_place1 = name_intermediate_directory + "/blast_n_makebd.txt"
    os.system(command1 + ">> " + redirect_place1)

    # Setting the name for the output file of the BLAST search
    name_output = name_intermediate_directory + "/denovo_blast_output_nucl.txt"

    # Constructing the BLAST command based on evalue provided
    location_denovo_nucl = name_output_directory + "/denovo_nucl.fa"
    command2 = (
        "blastn -evalue "
        + parameters_db["e_value"]
        + " -query "
        + location_denovo_nucl
        + " -db "
        + path_to_db
        + " -out "
        + name_output
        + ' -outfmt "7 qacc sacc evalue Identities qstart qend sstart send qcovs"'
    )

    # Executing the BLAST search command and redirecting output to a file
    redirect_place2 = name_intermediate_directory + "/blast_n_output.txt"
    os.system(command2 + ">> " + redirect_place2)
    os.system("rm " + path_to_db + ".*")


def perform_blast_prot_transcripts(
    name_intermediate_directory: str, name_output_directory: str
) -> None:
    """
    Performs a protein BLAST search using the DIAMOND tool between query and target protein sequences.

    This function:
    1. Creates a DIAMOND database for the target protein sequences.
    2. Executes a protein-to-protein BLAST (blastp) search of the query protein sequences against the target protein database.
    3. Stores the BLAST search results in a specified output file.
    4. Cleans up temporary DIAMOND database files after the search.

    Args:
        name_intermediate_directory (str): The directory where intermediate files will be stored, including the target and query protein sequences.
        name_output_directory (str): The directory containing the original denovo sequences (not directly used here but part of the system).

    Returns:
        None: This function does not return anything but generates output files in the intermediate directory.

    Side Effects:
        - Creates a DIAMOND database using the `diamond makedb` command.
        - Executes the DIAMOND `blastp` command to search for homologous proteins.
        - Generates BLAST output files (e.g., `diamond_transc_prot.out`, `blast_diamond_ORFs_makebd.txt`, `blast_diamond_ORFs.txt`).
        - Removes temporary DIAMOND database files (`*.dmnd`).
    """
    path_target_prot = (
        name_intermediate_directory + "/Intermediate_prot_BLAST/target_prot.fa"
    )
    command1 = "diamond makedb --in " + path_target_prot + " -d target_prot"
    path_output_1 = name_intermediate_directory + "/blast_diamond_ORFs_makebd.txt"
    os.system(command1 + ">> " + path_output_1)
    path_query = name_intermediate_directory + "/Intermediate_prot_BLAST/query_prot.fa"
    link_to_blast_output = name_intermediate_directory + "/diamond_transc_prot.out"
    command2 = (
        "diamond blastp -d target_prot -q "
        + path_query
        + " --more-sensitive -o "
        + link_to_blast_output
    )
    path_output_2 = name_intermediate_directory + "/blast_diamond_ORFs.txt"
    os.system(command2 + ">> " + path_output_2)
    final_rm = name_intermediate_directory + "/*.dmnd"


def perform_blast_prot_homology_filter(
    name_output_directory: str,
    name_intermediate_directory: str,
    path_to_db: str,
    parameters_db: dict,
) -> None:
    """
    Performs a protein BLAST search using the DIAMOND tool between query and target protein sequences.

    This function:
    1. Creates a DIAMOND database for the target protein sequences using the provided `path_to_db`.
    2. Executes a protein-to-protein BLAST (blastp) search of the query protein sequences (`denovo_protein.fa`) against the target protein database.
    3. Stores the BLAST search results in a specified output file.
    4. Cleans up temporary DIAMOND database files after the search.

    Args:
        name_output_directory (str): Directory containing the original denovo protein sequences (`denovo_protein.fa`).
        name_intermediate_directory (str): Directory where intermediate files (such as BLAST output and temporary databases) will be stored.
        path_to_db (str): Path to the target protein database used for the BLAST search.
        parameters_db (dict): Dictionary containing parameters for the BLAST search, including:
            - `"mode"`: Mode for the DIAMOND blastp command (e.g., `--more-sensitive`).

    Returns:
        None: This function does not return anything but generates output files in the intermediate directory.

    Side Effects:
        - Creates a DIAMOND database using the `diamond makedb` command.
        - Executes the DIAMOND `blastp` command to search for homologous proteins.
        - Generates BLAST output files (e.g., `denovo_blast_output_prot.txt`).
        - Removes temporary DIAMOND database files (`*.dmnd`).
    """
    # Extracting the name of the database from the path
    list_path = path_to_db.split("/")
    name_database = list_path[len(list_path) - 1].split(".")[0]
    name_database = name_intermediate_directory + "/" + name_database

    # Constructing the command to create a Diamond database
    command1 = "diamond makedb --in " + path_to_db + " -d " + name_database

    # Executing the command to create the Diamond database
    redirect_place1 = name_intermediate_directory + "/output_step2_prot_makeblastdb.txt"
    os.system(command1 + ">> " + redirect_place1)

    # Setting the name for the output file of the Diamond search
    name_output = name_intermediate_directory + "/denovo_blast_output_prot.txt"

    # Constructing the Diamond command based on parameters provided
    location_denovo_prot = name_output_directory + "/denovo_protein.fa"
    command2 = (
        "diamond blastp -d "
        + name_database
        + " -q "
        + location_denovo_prot
        + " "
        + parameters_db["mode"]
        + " -o "
        + name_output
    )

    # Executing the Diamond search command
    redirect_place2 = name_intermediate_directory + "/output_step2_prot_diamond.txt"
    os.system(command2 + ">> " + redirect_place2)

    # Removing temporary Diamond database files
    os.system("rm " + name_intermediate_directory + "/*.dmnd")


def perform_blast_to_genome(
    ID_genome: str, link_to_my_target_genome: str, name_output_directory: str
) -> int:
    """
    Performs a nucleotide BLAST search (BLASTn) of denovo genes candidates against a target genome.

    This function:
    1. Creates a BLAST database from the target genome.
    2. Executes a BLAST search from the query denovo gene sequences (`denovo_unspliced_lowered_introns.fa`) against the target genome.
    3. Stores the results in an output file.
    4. Cleans up temporary BLAST database files after the search.

    Args:
        ID_genome (str): The identifier for the genome being processed. Used to name the output file.
        link_to_my_target_genome (str): Path to the target genome file (in FASTA format) that will be used for BLAST search.
        name_output_directory (str): Directory where intermediate and output files will be stored.

    Returns:
        int: The size of the BLAST output file (in bytes). A value of `0` indicates that the output file is empty.

    Side Effects:
        - Creates a BLAST database from the target genome using the `makeblastdb` command.
        - Executes the `blastn` command to perform the search between denovo genes and the target genome.
        - Generates an output file with the BLAST search results (e.g., `*_all_blast_output.txt`).
        - Removes the temporary BLAST database files after the search.
    """
    name_intermediate_directory = (
        name_output_directory + "/Intermediate_output/blast_denovo_to_target_genome"
    )
    # Defining parameters for BLAST to genome
    parameters_blast_genome = {"type": "nucl", "e_value": "0.01", "coverage": None}

    # Setting the name for the output file of the BLAST search
    name_output = name_intermediate_directory + "/" + ID_genome + "_all_bast_output.txt"

    command_blast_db = "makeblastdb -in " + link_to_my_target_genome + " -dbtype nucl"
    # Creating a BLAST database from the target genome
    link_db_output = (
        name_output_directory
        + "/Intermediate_output/output_last_blast_to_genome_makedb.txt"
    )
    os.system(command_blast_db + ">> " + link_db_output)

    # Constructing the BLAST command based on parameters provided
    name_my_query = name_output_directory + "/denovo_unspliced_lowered_introns.fa"
    command = (
        "blastn -evalue "
        + parameters_blast_genome["e_value"]
        + " -query "
        + name_my_query
        + " -db "
        + link_to_my_target_genome
        + " -out "
        + name_output
        + ' -outfmt "7 qacc sacc evalue Identities qstart qend sstart send qcovs"'
    )

    # Executing the BLAST search command
    link_output_blast = (
        name_output_directory + "/Intermediate_output/output_last_blast_to_genome.txt"
    )
    os.system(command + ">> " + link_output_blast)
    check_file = os.stat(name_output).st_size

    # Removing temporary BLAST database files
    list_path_to_folder = link_to_my_target_genome.split("/")
    path_to_folder = ""
    for i in list_path_to_folder[0 : len(list_path_to_folder) - 1]:
        path_to_folder += i + "/"
    os.system("rm " + path_to_folder + "*.ndb")
    os.system("rm " + path_to_folder + "*.nin")
    os.system("rm " + path_to_folder + "*.not")
    os.system("rm " + path_to_folder + "*.nsq")
    os.system("rm " + path_to_folder + "*.ntf")
    os.system("rm " + path_to_folder + "*.nto")
    os.system("rm " + path_to_folder + "*.nhr")
    return check_file


def get_transcription_hits(
    dico_target_species_lines: dict, pop_species_name: str, name_output_directory: str
) -> dict:
    """
    Blasts all de novo candidates (nucleotide sequences without introns) against the target species' transcriptome.
    Matches are stored in a dictionary, categorized as 'complete' or 'partial' based on coverage.

    This function performs the following steps:
    1. Checks if the target species has a transcriptome GTF file and a corresponding transcriptome FASTA file.
    2. Creates a BLAST database from the transcriptome FASTA file.
    3. Performs a BLAST search of the query denovo nucleotide sequences against the target transcriptome.
    4. Categorizes the BLAST hits as 'complete' or 'partial' based on the coverage percentage.
    5. Stores the results in a dictionary where the key is the query denovo sequence and the value is another dictionary
       with target transcripts and their corresponding coverage status.

    Args:
        dico_target_species_lines (dict): A dictionary containing species-specific data, including transcriptome information.
        pop_species_name (str): The species name for which the BLAST search is performed.
        name_output_directory (str): Directory to store intermediate and output files during the process.

    Returns:
        dict: A dictionary containing BLAST hit information. The keys are the de novo gene IDs, and the values are dictionaries
              where keys are the transcript IDs and values are 'complete' or 'partial' based on the coverage.
    """
    dico_transcript_hits = {}  # Dictionary to store transcript hits

    # Checking if transcriptome GTF file exists in the dictionary for the given species
    if "transcriptome_gtf" in dico_target_species_lines[pop_species_name]:
        # Obtaining the path to the transcriptome FASTA file
        link_transcripts_fasta = dico_target_species_lines[pop_species_name][
            "transcriptome_fasta"
        ]

        # Creating a BLAST database from the transcriptome FASTA file
        command1 = "makeblastdb -in " + link_transcripts_fasta + " -dbtype nucl"
        link_transc_db = (
            name_output_directory
            + "/Intermediate_output/output_last_blast_denovo_to_homolog_db.txt"
        )
        os.system(command1 + ">> " + link_transc_db)

        # Performing BLAST of the query denovo genes against the target transcriptome
        link_denovo_nucl = name_output_directory + "/denovo_nucl.fa"
        link_blast_output = (
            name_output_directory + "/Intermediate_output/transcript_blast_output.txt"
        )
        command2 = (
            "blastn -evalue 0.01 -query "
            + link_denovo_nucl
            + " -db "
            + link_transcripts_fasta
            + " -out "
            + link_blast_output
            + ' -outfmt "7 qacc sacc evalue Identities qstart qend sstart send qcovs"'
        )
        link_transc_blast_out = (
            name_output_directory
            + "/Intermediate_output/output_last_blast_denovo_to_homolog.txt"
        )
        os.system(command2 + ">> " + link_transc_blast_out)

        # Removing temporary BLAST database files
        list_path_to_folder = link_transcripts_fasta.split("/")
        path_to_folder = ""
        for i in list_path_to_folder[0 : len(list_path_to_folder) - 1]:
            path_to_folder += i + "/"
        os.system("rm " + path_to_folder + "*.ndb")
        os.system("rm " + path_to_folder + "*.nin")
        os.system("rm " + path_to_folder + "*.not")
        os.system("rm " + path_to_folder + "*.nsq")
        os.system("rm " + path_to_folder + "*.ntf")
        os.system("rm " + path_to_folder + "*.nto")
        os.system("rm " + path_to_folder + "*.nhr")

        # Processing BLAST output to extract hits
        my_blast_hits = openFile(link_blast_output)
        for line in my_blast_hits:
            if line[0] != "#":
                elts_line = line.split()
                query_denovo = elts_line[0]
                target_transcript = elts_line[1]
                cov = int(elts_line[7])
                value_cov = ""

                # Determining coverage type
                if cov == 100:
                    value_cov = "complete"
                else:
                    value_cov = "partial"

                # Storing hits in the dictionary
                if query_denovo not in dico_transcript_hits:
                    dico_transcript_hits[query_denovo] = {target_transcript: value_cov}
                else:
                    if target_transcript not in dico_transcript_hits[query_denovo]:
                        dico_transcript_hits[query_denovo][target_transcript] = (
                            value_cov
                        )
        # Removing temporary BLAST output file
        os.system("rm " + link_blast_output)

    return dico_transcript_hits
