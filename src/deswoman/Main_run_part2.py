from Bio.SeqRecord import SeqRecord
import os
from deswoman.module_blast_and_diamond import (
    perform_blast_prot_homology_filter,
    perform_blast_nucl_homology_filter,
)
from deswoman.module_select_orphans import (
    retrieve_name_hit_prot_blast,
    retrieve_name_hit_nucl_blast,
    merge_hit_lists,
    get_file_size,
    reshufe_files_in_denovo,
)
from deswoman.module_colors import *


__author__ = "Anna Grandchamp"
__contributor__ = ""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def run_part_2(dico_variables: dict) -> None:
    """
    Performs homology search of candidate de novo sequences against external protein and nucleotide databases.

    This function runs two separate homology searches: one against a protein database and another
    against a nucleotide (DNA/RNA) database. It uses BLAST to find matching sequences and removes
    any de novo sequences that show significant homology. After homology searching, the function
    updates the output files to reflect only the candidate sequences that have no homology hits.

    Args:
        dico_variables (dict): A dictionary containing user-defined parameters, including:
            - "path_output": Directory where output files are saved.
            - "path_output_intermediate": Directory for intermediate output files.
            - "link_database_outgroup_prot": Path to the protein database for BLAST search (optional).
            - "link_database_outgroup_nucl": Path to the nucleotide database for BLAST search (optional).
            - "parameters_database_prot": Parameters for protein BLAST search (optional).
            - "parameters_database_nucl": Parameters for nucleotide BLAST search (optional).

    Returns:
        None: This function modifies the output files directly and does not return a value.

    Workflow:
        1. **Protein Homology Search**: BLAST candidate de novo protein sequences against a protein database.
        2. **Nucleotide Homology Search**: BLAST candidate de novo nucleotide sequences against a nucleotide database.
        3. **Merge Results**: Combines the results of both BLAST searches and determines which sequences to remove.
        4. **Update Output Files**: Removes de novo sequences with significant homology from the output files.

    Notes:
        - If no homology database is provided, the corresponding search step will be skipped.
        - A final count of the remaining de novo candidates is displayed.
    """
    # lists ready to store names of candidate ORFs that show homology and will be discarded by the end of the step.
    list_name_to_remove_nucl = []
    list_name_to_remove_prot = []
    list_name_to_remove_total = []
    name_output_directory = dico_variables["path_output"]
    name_intermediate_directory = dico_variables["path_output_intermediate"]
    # conduct homology search against protein dataset only if a protein file is given
    if dico_variables["link_database_outgroup_prot"] != "":
        print(BRIGHT_BLUE + "Performing protein homology search ..." + RESET)
        # perform diamond BLAST of candidate de novo prot against protein dataset and stocks the results in intermediate directory
        perform_blast_prot_homology_filter(
            name_output_directory,
            name_intermediate_directory,
            dico_variables["link_database_outgroup_prot"],
            dico_variables["parameters_database_prot"],
        )
        # if hits were detected, the name of the query denovo with hits are stoded in the list "list_name_to_remove_prot"
        list_name_to_remove_prot = retrieve_name_hit_prot_blast(
            name_intermediate_directory
        )

    # conduct homology search against protein dataset only if a nucleotide/RNA file is given
    if dico_variables["link_database_outgroup_nucl"] != "":
        print(BRIGHT_BLUE + "Performing DNA/RNA homology search ..." + RESET)
        # perform nucl BLAST of candidate de novo nucl against protein dataset and stocks the results in intermediate directory
        perform_blast_nucl_homology_filter(
            name_output_directory,
            name_intermediate_directory,
            dico_variables["link_database_outgroup_nucl"],
            dico_variables["parameters_database_nucl"],
        )
        # if hits were detected, the name of the query denovo with hits are stoded in the list "list_name_to_remove_nucl"
        list_name_to_remove_nucl = retrieve_name_hit_nucl_blast(
            name_intermediate_directory
        )
    # merge the 2 lists (uniq name is common to both)
    list_name_to_remove_total = merge_hit_lists(
        list_name_to_remove_prot, list_name_to_remove_nucl
    )

    file_denovo_nucl = name_output_directory + "/denovo_nucl.fa"
    file_denovo_prot = name_output_directory + "/denovo_protein.fa"
    # get size of denovo_nucl (and prot) before removing ORFs that had homology
    dico_size_denovo_nucl = get_file_size(file_denovo_nucl)
    dico_size_denovo_prot = get_file_size(file_denovo_prot)

    print(
        BRIGHT_BLUE
        + "Number of denovo candidate with a hit : "
        + RESET
        + str(len(list_name_to_remove_total))
        + " out of "
        + str(len(dico_size_denovo_prot))
    )
    new_nb = len(dico_size_denovo_prot) - len(list_name_to_remove_total)
    print(
        BRIGHT_BLUE
        + "New number of denovo candidate (without homology) : "
        + RESET
        + str(new_nb)
    )
    # here are re made all output files from step 1 with only denovo genes validated in here step 2.
    reshufe_files_in_denovo(name_output_directory, list_name_to_remove_total)
