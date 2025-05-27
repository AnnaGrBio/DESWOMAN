import os
from Bio import SeqIO
from deswoman.module_colors import openFile


def build_target_and_query_prot_list(
    name_intermediate_directory: str, name_output_directory: str
) -> None:
    """
    Creates two identical FASTA files containing all denovo proteins extracted from previous steps: one for the query and one for the target.
    These files are stored in the 'Intermediate_prot_BLAST' directory.

    Parameters:
        name_intermediate_directory (str): Path to the intermediate directory where the "Intermediate_prot_BLAST" folder will be created.
        name_output_directory (str): Path to the directory containing the "denovo_protein.fa" file with denovo proteins.

    Returns:
        None: The function writes the query and target protein files into the 'Intermediate_prot_BLAST' folder.
    """
    dico_all_seq = {}
    new_path = name_intermediate_directory + "/Intermediate_prot_BLAST"
    # if the folder does not exist as such, it is created
    if os.path.exists(new_path) is False:
        os.system("mkdir " + new_path)
    # path to all de novo prots extracted before
    path_denovo_prot = name_output_directory + "/denovo_protein.fa"
    for seq_record in SeqIO.parse(path_denovo_prot, "fasta"):
        dico_all_seq[str(seq_record.id)] = str(seq_record.seq)
    # all proteins are stored in a query and target file
    path_to_query_prot = new_path + "/query_prot.fa"
    path_to_target_prot = new_path + "/target_prot.fa"
    query_prot_file = open(path_to_query_prot, "w")
    target_prot_file = open(path_to_target_prot, "w")
    # fill the 2 files exactly in the same way.
    for prot_name in dico_all_seq:
        query_prot_file.write(">" + prot_name + "\n")
        target_prot_file.write(">" + prot_name + "\n")
        query_prot_file.write(dico_all_seq[prot_name] + "\n")
        target_prot_file.write(dico_all_seq[prot_name] + "\n")
    query_prot_file.close()
    target_prot_file.close()


def extract_coord(link_file: str) -> dict:
    """
    Parses the information file generated in Step 1 and 2, and extracts the unspliced position and chromosome of all denovo sequences.

    Parameters:
        link_file (str): Path to the information file containing denovo sequence data.

    Returns:
        dict: A dictionary where the keys are the sequence names and the values are lists containing chromosome, start position, and stop position.
    """
    dico = {}
    info_file = openFile(link_file)
    for line in info_file[1:]:
        elts_line = line.split(",")
        chrom = elts_line[1]
        name = elts_line[9]
        start = int(elts_line[10])
        stop = int(elts_line[11])
        dico[name] = [chrom, start, stop]
    return dico


def assess_overlap(l1_old: list, l2_old: list) -> bool:
    """
    This function checks if two given intervals (representing suspected orthologs) overlap on the same chromosome.

    Parameters:
        l1_old (list): A list of three elements representing the first interval. The format is [chromosome, start, stop].
        l2_old (list): A list of three elements representing the second interval. The format is [chromosome, start, stop].

    Returns:
        bool: True if the two intervals overlap, otherwise False.
    """
    # first, range the position from the smaller to the highest (ex if start is 50 and stop 40 it reverse them).
    if l1_old[1] > l1_old[2]:
        l1 = [l1_old[0], l1_old[2], l1_old[1]]
    else:
        l1 = [l1_old[0], l1_old[1], l1_old[2]]
    if l2_old[1] > l2_old[2]:
        l2 = [l2_old[0], l2_old[2], l2_old[1]]
    else:
        l2 = [l2_old[0], l2_old[1], l2_old[2]]
    overlap = False
    if l1[0] == l2[0]:
        start1 = l1[1]
        start2 = l2[1]
        stop1 = l1[2]
        stop2 = l2[2]
        if start1 <= start2 and stop1 >= start2:
            overlap = True
        elif start1 <= start2 and stop1 >= stop2:
            overlap = True
        elif start1 >= start2 and stop1 <= stop2:
            overlap = True
        elif stop1 >= stop2 and start1 <= stop2:
            overlap = True
        elif start1 == start2 and stop1 == stop2:
            overlap = True
    return overlap


def reduce_by_location(list_orthologs: list, query: str, dico_all_coord: dict) -> dict:
    """
    This function filters a list of suspected homologs for a given query based on genomic location,
    keeping only those that overlap with the query's genomic coordinates.

    Parameters:
        list_orthologs (list): A list of suspected homologs (target genes).
        query (str): The name of the query gene.
        dico_all_coord (dict): A dictionary containing the genomic coordinates (chromosome, start, stop)
                               for each gene, with gene names as keys.

    Returns:
        dict: A dictionary where the query gene is the key, and the value is a list of homologs
              that overlap with the query's genomic coordinates.
    """
    correct_dico_homologs = {}

    for target in list_orthologs:
        # assess the overlap between the query and its targets
        overlap = assess_overlap(dico_all_coord[query], dico_all_coord[target])
        if overlap is True:
            if query not in correct_dico_homologs:
                correct_dico_homologs[query] = [target]
            else:
                correct_dico_homologs[query].append(target)
    # if the query has no homolog or no homolog true : it becomes associated to an empty list.
    if len(correct_dico_homologs) == 0:
        correct_dico_homologs[query] = []
    return correct_dico_homologs


def create_dico_all_prot_hit(my_file: list) -> dict:
    """
    This function creates a dictionary where each key is a query (denovo protein) and each value is another
    dictionary representing its BLAST hits (targets). Hits to the query itself are removed.

    Parameters:
        my_file (file object in readline mode): A file object containing the BLAST results, where each line represents a hit.
                               The first column is the query, and the second column is the target.

    Returns:
        dict: A dictionary with queries as keys and their corresponding BLAST hits as values. Each value is
              a dictionary where the target proteins are keys and all have a value of 1 (indicating a hit).
    """
    dico_all_hits = {}
    for line in my_file:
        elts_line = line.split()
        query = elts_line[0]
        target = elts_line[1]
        if query in dico_all_hits:
            dico_all_hits[query][target] = 1
        else:
            dico_all_hits[query] = {target: 1}
    # remove hits to itself
    for query in dico_all_hits:
        if query in dico_all_hits[query]:
            del dico_all_hits[query][query]
    return dico_all_hits


def fill_dico_orthogroups(dico_my_hits: dict, dico_orthogroups: dict) -> None:
    """
    Populates the dictionary of orthogroups with new orthogroups derived from BLAST hits.

    Parameters:
        dico_my_hits (dict): A dictionary where keys are query proteins, and values are lists of their BLAST hits.
        dico_orthogroups (dict): A dictionary storing orthogroups, where keys are orthogroup names (e.g., "orthogroup1"),
                                 and values are lists of protein sequences belonging to that group.

    Returns:
        None: The function modifies dico_orthogroups in place by adding new orthogroups.
    """
    nb = len(dico_orthogroups) + 1
    for denovo in dico_my_hits:
        my_orthogroup = dico_my_hits[denovo]
        my_orthogroup.append(denovo)
        newname = "orthogroup" + str(nb)
        dico_orthogroups[newname] = my_orthogroup


def handle_blast_inputs_transc(
    dico_orthogroups: dict,
    dico_name_size_denovo: dict,
    name_intermediate_directory: str,
    name_output_directory: str,
) -> None:
    """
    Processes BLAST results to identify and categorize ORFS orthogroups.

    This function:
    - Extracts unspliced coordinates of all de novo ORFs.
    - Parses BLAST results to determine homologous relationships.
    - Ensures reciprocal BLAST validation for ortholog assignment.
    - Filters orthologous groups based on genomic overlap.
    - Removes unnecessary intermediary files.
    - Populates the orthogroup dictionary with validated orthologs.

    Parameters:
        dico_orthogroups (dict): A dictionary storing orthogroups.
        dico_name_size_denovo (dict): A dictionary mapping de novo transcript names to their sizes.
        name_intermediate_directory (str): Path to the intermediate directory containing BLAST results.
        name_output_directory (str): Path to the outgroup directory with de novo transcript data.

    Returns:
        None: Modifies dico_orthogroups in place by adding validated ortholog groups.
    """
    link_info_file = name_output_directory + "/information_file.txt"
    # extracts unspliced coordinates of all denovo
    dico_all_coord = extract_coord(link_info_file)

    dico_already_attributed_hit = {}
    link_diamond_transc_prot = name_intermediate_directory + "/diamond_transc_prot.out"
    check_file = os.stat(link_diamond_transc_prot).st_size
    # if the blast between query and target protein has output (it should mandatorily because the query and targte contain the same proteins)
    if check_file > 0:
        file_blast = openFile(link_diamond_transc_prot)
        # create a dico with query as keys and their hts as target (hit to themselves are removed).
        dico_all_hits = create_dico_all_prot_hit(file_blast)
        # remove files that became unusefull.
        os.system("rm " + name_intermediate_directory + "/diamond_transc_prot.out")
        os.system(
            "rm "
            + name_intermediate_directory
            + "/Intermediate_prot_BLAST/query_prot.fa"
        )
        os.system(
            "rm "
            + name_intermediate_directory
            + "/Intermediate_prot_BLAST/target_prot.fa"
        )
        for query in dico_all_hits:  ## here retrieve reciprocal blast
            if query not in dico_already_attributed_hit:
                suspected_orthogroup = []
                # generate a dic with all query hits
                dico_homologs_query = dico_all_hits[query]
                for homolog in dico_homologs_query:
                    # here make sure of the reciprocal blast and the fact that the homolog does not already belong to a group where the query is not
                    if (
                        homolog not in dico_already_attributed_hit
                        and query in dico_all_hits[homolog]
                    ):
                        suspected_orthogroup.append(homolog)
                ## reduce orthogroup by overlapping location. If suspected_orthogroup suspected orthogroup is empty, the query is associated to an empty list in the dictionary, as well as if no homologs fits.
                dico_correct_orthogroup = {}
                dico_correct_orthogroup = reduce_by_location(
                    suspected_orthogroup, query, dico_all_coord
                )
                dico_already_attributed_hit[query] = 1
                if len(dico_correct_orthogroup[query]) > 0:
                    for target in dico_correct_orthogroup[query]:
                        dico_already_attributed_hit[target] = 1
                fill_dico_orthogroups(dico_correct_orthogroup, dico_orthogroups)

    dico_no_hit_denovo = {}
    for denovo in dico_name_size_denovo:
        if denovo not in dico_already_attributed_hit:
            dico_no_hit_denovo[denovo] = []
    if len(dico_no_hit_denovo) > 0:
        for i in (
            dico_no_hit_denovo
        ):  # ChangeMarie ##Prevent overwriting the elements in that dictionary
            dico_no_hit_small = {
                i: dico_no_hit_denovo[i]
            }  # ChangeMarie ##A really dumb solution but it works
            fill_dico_orthogroups(
                dico_no_hit_small, dico_orthogroups
            )  # ChangeMarie ##Adapt the dictionary name


def make_dico_name_denovo_strat2(name_output_directory: str) -> dict:
    """
    Extracts the names of all de novo proteins.

    This function parses a FASTA file containing de novo protein sequences and
    stores their identifiers in a dictionary.

    Parameters:
        name_output_directory (str): Path to the directory containing the 'denovo_protein.fa' file.

    Returns:
        dict: A dictionary where keys are de novo protein names (sequence IDs) and values are set to 1.
    """
    my_denovo_file = name_output_directory + "/denovo_protein.fa"
    dico_name_denovo = {}
    for seq_record in SeqIO.parse(my_denovo_file, "fasta"):
        dico_name_denovo[str(seq_record.id)] = 1
    return dico_name_denovo


def build_final_file_strat2_step3(
    dico_orthogroups: dict, name_intermediate_directory: str, name_output_directory: str
) -> None:
    """
    Builds the final output file for step 3 of strategy 2.

    This function removes the temporary directory used for intermediate protein BLAST results
    and writes the final orthogroup assignments to an output file.

    Parameters:
        dico_orthogroups (dict): A dictionary containing orthogroups, where keys are orthogroup names
                                 and values are lists of protein names.
        name_intermediate_directory (str): Path to the intermediate directory.
        name_output_directory (str): Path to the outgroup directory where the final file is stored.

    Returns:
        None: The function writes the results to 'Orthogroup_output_step3.txt' in the outgroup directory.
    """
    link_dir = name_intermediate_directory + "/Intermediate_prot_BLAST"
    os.system("rm -r " + link_dir)
    link_new_output_file = name_output_directory + "/Orthogroup_output_step3.txt"
    new_file = open(link_new_output_file, "w")
    for orthogroup in dico_orthogroups:
        new_file.write(orthogroup)
        for names in dico_orthogroups[orthogroup]:
            if names == dico_orthogroups[orthogroup][0]:
                new_file.write(":" + names)
            else:
                new_file.write("," + names)
        new_file.write("\n")
    new_file.close()
