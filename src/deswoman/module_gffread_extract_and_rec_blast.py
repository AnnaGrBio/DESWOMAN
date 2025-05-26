import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from deswoman.module_transcripts_data_and_threeshold import get_path_to_fasta, get_path_gff
from deswoman.module_colors import openFile

__author__ = "Anna Grandchamp"
__contributor__ = ""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def reshuffle_prot_file(link_prot1):
    """
    This function does not serve, i could just have renamed the link to prot 1 in new file. It was written with a version of gffread where
    a modif had to be perfored in the output file.
    Still, I leave it that way also for if in future gffread evolve and modifications have to be done again.
    """
    # Opening the protein file for reading
    my_file1 = openFile(link_prot1)
    # Creating a new file name for the reshuffled protein file
    name_prot1_reshuffled = link_prot1.split(".")[0] + "_new.fa"

    # Opening the new file for writing
    my_new_file1 = open(name_prot1_reshuffled, "w")

    # Iterating through each line of the protein file
    for line in my_file1:
        my_new_file1.write(line)  # new

    my_new_file1.close()

    # Returning the name of the reshuffled protein file
    return name_prot1_reshuffled


def extract_all_prots_genome(
    pop_ind_name: str, path_to_genome: str, name_intermediate_directory: str
) -> None:
    """
    Extracts all protein sequences from annotated genes in a genome based on the provided GFF annotation.

    This function checks if the protein sequences have already been extracted for a given population/individual
    (specified by `pop_ind_name`). If they haven't been extracted yet, it performs the extraction using the `gffread` tool,
    which takes the genome's GFF and FASTA files as input. The protein sequences are then reshuffled and saved into a
    designated directory.

    The extracted protein file is only generated once for each genome. If the file already exists,
    the function skips the extraction process.

    Args:
        pop_ind_name (str): The name of the population or individual for which protein sequences should be extracted.
                             This name is used to identify the appropriate genome and annotation files.
        path_to_genome (str): The directory path containing the genome files (FASTA and GFF) for the species.
        name_intermediate_directory (str): The directory where intermediate files (e.g., extracted protein files) will be stored.

    Returns:
        None: This function does not return any values. It performs file operations (creating directories, extracting proteins, etc.).
    """
    file_to_assess = (
        name_intermediate_directory + "/old_prot_blast/" + pop_ind_name + "_prot_new.fa"
    )

    # Checking if the file exists. It is convenient because for each pops, the query genes are required. It would be unusefull to create them each time, only once is enough.
    if os.path.isfile(file_to_assess) == False:
        # Creating a directory if it doesn't exist
        path_to_create_old_prot = name_intermediate_directory + "/old_prot_blast"
        if os.path.isdir(path_to_create_old_prot) == False:
            os.system("mkdir " + path_to_create_old_prot)

        # Constructing paths to required files
        link_output = path_to_genome + "/" + pop_ind_name + "_prot.fa"
        link_genome = get_path_to_fasta(path_to_genome, pop_ind_name)
        link_gff = get_path_gff(path_to_genome, pop_ind_name)
        # Constructing the command to extract proteins from GFF annotation
        command = "gffread -y " + link_output + " -g " + link_genome + " " + link_gff
        # Executing the command
        os.system(command)

        # Deleting temporary files
        file_to_delete = link_genome + ".fai"
        os.system("rm " + file_to_delete)
        # print ("gffread extraction done")

        # Reshuffling the protein file
        link_final_output = reshuffle_prot_file(link_output)
        # print ("filtering of established proteins done")
        os.system("rm " + link_output)

        # Moving the final output file to the appropriate directory
        os.system("mv " + link_final_output + " " + path_to_create_old_prot)


def filter_prot_file(link_to_file: str, link_to_new_file: str) -> None:
    """
    Filters out stop codons represented by dots ('.') from a protein sequence file.

    In some protein sequences, stop codons are represented by dots (with an old version of gffread), which can be problematic in downstream steps of
    the analysis. This function removes these dots to ensure that the sequences can be processed correctly by DESMAN.

    The primary goal of DESMAN is not to determine whether a protein is complete or broken, but to ensure that
    sequences are clean for further analysis. The dots are removed to avoid any issues with downstream steps.

    Args:
        link_to_file (str): The path to the input protein sequence file that contains stop codons represented by dots.
        link_to_new_file (str): The path to the new file where the filtered protein sequences will be saved, with dots removed.

    Returns:
        None: The function writes the filtered sequences directly to the specified output file.
    """
    old_file = openFile(link_to_file)
    new_file = open(link_to_new_file, "w")
    for line in old_file:
        if line[0] == ">":
            new_file.write(line)
        else:
            for i in line:
                if i != ".":
                    new_file.write(i)
    new_file.close()


def make_reciprocal_prot_blasts(
    query_name: str, target_name: str, name_intermediate_directory: str
) -> None:
    """
    Performs reciprocal protein BLAST searches between query and target protein sets.

    This function prepares protein sequences by filtering out stop codons ('.'), then creates Diamond
    databases for both query and target proteins. Subsequently, it performs reciprocal BLAST searches
    between the two protein sets to identify homologous sequences.

    Args:
        query_name (str): Name of the query species for the BLAST search.
        target_name (str): Name of the target species for the BLAST search.
        name_intermediate_directory (str): Path to the intermediate directory where protein files and results will be stored.

    Returns:
        None: The function performs the necessary steps and writes results to output files.
    """
    # Constructing paths to query and target protein files
    ## Here the links to old and new represent the same file, but we use 2 names to rewrite them.
    link_old_prot_query = (
        name_intermediate_directory + "/old_prot_blast/" + query_name + "_prot_new.fa"
    )
    link_old_prot_target = (
        name_intermediate_directory + "/old_prot_blast/" + target_name + "_prot_new.fa"
    )

    link_prot_query = (
        name_intermediate_directory + "/old_prot_blast/" + query_name + "_prot_new.fa"
    )
    link_prot_target = (
        name_intermediate_directory + "/old_prot_blast/" + target_name + "_prot_new.fa"
    )

    filter_prot_file(link_old_prot_query, link_prot_query)
    filter_prot_file(link_old_prot_target, link_prot_target)

    # Creating Diamond databases for query and target proteins
    output_file_query = name_intermediate_directory + "/" + query_name
    output_file_target = name_intermediate_directory + "/" + target_name
    command1 = (
        "./diamond makedb --in " + link_prot_query + " -d " + output_file_query
    )  # query_name
    command2 = (
        "./diamond makedb --in " + link_prot_target + " -d " + output_file_target
    )  # target_name
    redirect_comm1 = (
        name_intermediate_directory + "/output_last_diamond_gene_rec_blast_makedb_1.txt"
    )
    redirect_comm2 = (
        name_intermediate_directory + "/output_last_diamond_gene_rec_blast_makedb_2.txt"
    )
    os.system(command1 + ">> " + redirect_comm1)
    os.system(command2 + ">> " + redirect_comm2)
    # Performing reciprocal protein BLAST searches
    name_output_query = (
        name_intermediate_directory + "/old_prot_blast/prot_query_blast_out.txt"
    )
    name_output_target = (
        name_intermediate_directory + "/old_prot_blast/prot_target_blast_out.txt"
    )
    command_blast1 = (
        "./diamond blastp -d "
        + output_file_target
        + " -q "
        + link_prot_query
        + " -o "
        + name_output_query
    )  # -d was target_name
    command_blast2 = (
        "./diamond blastp -d "
        + output_file_query
        + " -q "
        + link_prot_target
        + " -o "
        + name_output_target
    )  # -d was query_name
    redirect_blast1 = (
        name_intermediate_directory + "/output_last_diamond_gene_rec_blast_1.txt"
    )
    redirect_blast2 = (
        name_intermediate_directory + "/output_last_diamond_gene_rec_blast_2.txt"
    )
    os.system(command_blast1 + ">> " + redirect_blast1)
    os.system(command_blast2 + ">> " + redirect_blast2)


def make_dico_blast_for_reciprocal_hits(
    rec_best_hit: str, name_intermediate_directory: str
) -> (dict, dict):
    """
    Creates two dictionaries from reciprocal BLAST output files:
    one mapping query genes to target genes, and another mapping target genes to query genes.

    If the reciprocal best hit option is set to "True", only the best first hit is included.
    Otherwise, all hits are considered.

    Args:
        rec_best_hit (str): A string that is either "True" or "False" indicating whether to keep only the best reciprocal hit.
        name_intermediate_directory (str): Path to the intermediate directory containing BLAST output files.

    Returns:
        tuple: Two dictionaries:
            - dico_blast_query_to_target: Maps query genes to target genes.
            - dico_blast_target_to_query: Maps target genes to query genes.
    """
    dico_blast_query_to_target = {}
    # Opening the first protein BLAST output file for reading
    file_name_query = (
        name_intermediate_directory + "/old_prot_blast/prot_query_blast_out.txt"
    )
    file_blast1 = openFile(file_name_query)
    for line in file_blast1:
        if line[0] != "#":  # Skipping comment lines
            elts_line = line.split()
            gene_query = elts_line[0]
            gene_target = elts_line[1]
            # Adding target gene to dictionary under query gene key
            if gene_query not in dico_blast_query_to_target:
                dico_blast_query_to_target[gene_query] = [gene_target]
            else:
                # will implement all other less good hits only if the user does not want the best reciprocal hit
                if rec_best_hit != "True":
                    dico_blast_query_to_target[gene_query].append(gene_target)

    # Creating an empty dictionary for the second protein BLAST results
    dico_blast_target_to_query = {}
    file_name_target = (
        name_intermediate_directory + "/old_prot_blast/prot_target_blast_out.txt"
    )
    file_blast2 = openFile(file_name_target)
    # Iterating through each line in the second BLAST output file
    for line in file_blast2:
        if line[0] != "#":  # Skipping comment lines
            elts_line = line.split()
            gene_target = elts_line[
                0
            ]  # Extracting target genes that where here query of the blast
            gene_query = elts_line[
                1
            ]  # Extracting query genes that were here target of the blast
            # Adding target gene to dictionary under query gene key
            if gene_target not in dico_blast_target_to_query:
                dico_blast_target_to_query[gene_target] = [gene_query]
            else:
                # will implement all other less good hits only if the user does not want the best reciprocal hit
                if rec_best_hit != "True":
                    dico_blast_target_to_query[gene_target].append(gene_query)
    # Removing temporary BLAST output files
    link1 = name_intermediate_directory + "/old_prot_blast/prot_query_blast_out.txt"
    link2 = name_intermediate_directory + "/old_prot_blast/prot_target_blast_out.txt"
    os.system("rm " + link1)
    os.system("rm " + link2)
    # Returning both dictionaries
    return dico_blast_query_to_target, dico_blast_target_to_query
