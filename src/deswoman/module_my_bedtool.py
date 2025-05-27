from deswoman.module_sort_gff_gene_order import sort_gff_dic
from deswoman.module_colors import *


def assess_gff_markers(my_gff: list) -> (str, str):
    """
    Assesses the markers in a GFF file to determine the gene and exon markers.
    Parameters:
    my_gff (List[str]): A list of strings, each representing a line from a GFF file.
    Returns:
    str, str: The gene marker and the exon marker.
    """
    list_all_elts_2 = []
    # goes accross all lines of the genome annotation file
    for line in my_gff:
        if line[0] != "#":
            elts_line = line.split()
            if len(elts_line) >= 2:
                elts_2 = elts_line[2]
                # insert all elements from the genomic file in the list without redondancy
                if elts_2 not in list_all_elts_2:
                    list_all_elts_2.append(elts_2)
    gene_marker = "gene"
    # if exon in the list : marker is exon. Otherwise, it is supposed to be CDS
    if "exon" in list_all_elts_2:
        exon_marker = "exon"
    else:
        exon_marker = "CDS"
    return gene_marker, exon_marker


def build_object_ID(elts: list) -> list:
    """
    Constructs a list of identifiers from genomic feature information.

    This function takes an input list `elts`, which contains genomic feature data,
    and extracts specific elements to build a sublist representing a genomic feature
    (e.g., an exon or intron). The sublist contains the chromosome/scaffold name,
    start position, stop position, and orientation of the feature.

    Parameters:
    -----------
    elts : list
        A list containing the genomic feature details. It is assumed that `elts` has
        the following structure:
        - `elts[0]`: chromosome or scaffold name (string)
        - `elts[3]`: start position (string or integer convertible)
        - `elts[4]`: stop position (string or integer convertible)
        - `elts[6]`: orientation (string, typically '+' or '-')

    Returns:
    --------
    list
        A list with 4 elements:
        - Chromosome or scaffold name (string)
        - Start position (integer)
        - Stop position (integer)
        - Orientation (string, either '+' or '-')
    """
    chrom_scaffold = elts[0]
    start = int(elts[3])
    stop = int(elts[4])
    orientation = elts[6]
    sublist = [chrom_scaffold, start, stop, orientation]
    return sublist


def extract_introns(list_exons: list) -> list:
    """
    Extracts the intron coordinates from a list of exons.

    This function takes a list of exons (each represented by a sublist with chromosomal
    coordinates and orientation) and calculates the intronic regions between consecutive
    exons. The function supports both forward and reverse orientation, adjusting the
    intron coordinates accordingly. It returns a list of introns, where each intron is
    represented by a sublist with the chromosome/scaffold, start position, stop position,
    and orientation.

    Parameters:
    -----------
    list_exons : list of lists
        A list where each element represents an exon. Each exon is a list with the following elements:
        - `elts[0]`: chromosome or scaffold name (string)
        - `elts[1]`: start position (integer)
        - `elts[2]`: stop position (integer)
        - `elts[3]`: orientation (string, either '+' or '-')

    Returns:
    --------
    list of lists
        A list of intron coordinates, where each intron is represented as a list:
        - Chromosome or scaffold name (string)
        - Start position (integer)
        - Stop position (integer)
        - Orientation (string, either '+' or '-')
    """
    list_introns = []
    # if the exons are forward.
    if list_exons[0][3] == "+":
        for my_index in range(len(list_exons) - 1):
            start_intron = list_exons[my_index][2] + 1
            stop_intron = list_exons[my_index + 1][1] - 1
            sublist_intron = [
                list_exons[my_index][0],
                start_intron,
                stop_intron,
                list_exons[my_index][3],
            ]
            list_introns.append(sublist_intron)
    # if the exons are reverse.
    else:
        for my_index in range(len(list_exons) - 1):
            stop_intron = list_exons[my_index][1] - 1
            start_intron = list_exons[my_index + 1][2] + 1
            sublist_intron = [
                list_exons[my_index][0],
                start_intron,
                stop_intron,
                list_exons[my_index][3],
            ]
            list_introns.append(sublist_intron)

    return list_introns


def dico_object_per_chrom(list_elts: list) -> dict:
    """
    Create a dictionary with chromosomes/scaffolds as keys and lists of exon or intron properties as values.

    This function takes a list of elements where each element contains information about
    a genomic feature (e.g., exon or intron) such as the chromosome/scaffold name, start
    position, stop position, and orientation. It then organizes these elements into a dictionary
    where the keys are the chromosome or scaffold names and the values are lists of sublists,
    each containing the start position, stop position, and orientation of the features.

    Parameters:
    -----------
    list_elts : list of lists
        A list of elements where each element represents a genomic feature (exon or intron) and is
        represented as a sublist with the following structure:
        - `elts[0]`: chromosome or scaffold name (string)
        - `elts[1]`: start position (integer)
        - `elts[2]`: stop position (integer)
        - `elts[3]`: orientation (string, either '+' or '-')

    Returns:
    --------
    dict
        A dictionary where the keys are chromosome or scaffold names (strings) and the values are
        lists of sublists representing the genomic features. Each sublist contains:
        - Start position (integer)
        - Stop position (integer)
        - Orientation (string, either '+' or '-')

    Example:
    --------
    input = [["chrom1", 100, 200, "+"], ["chrom1", 300, 400, "+"], ["chrom2", 150, 250, "-"]]
    output = dico_object_per_chrom(input)
    print(output)
    # Output: {'chrom1': [[100, 200, '+'], [300, 400, '+']], 'chrom2': [[150, 250, '-']]}
    """
    dico = {}
    for sublist in list_elts:
        new_sublist = [sublist[1], sublist[2], sublist[3]]
        if sublist[0] not in dico:
            dico[sublist[0]] = [new_sublist]
        else:
            dico[sublist[0]].append(new_sublist)
    return dico


def define_start_gtf(my_gff: list) -> int:
    """
    Find the line number where the actual data starts in a GFF/GTF file, skipping header lines.

    This function reads a GFF or GTF file, line by line, and determines the first line that contains
    data (i.e., not a header line starting with "#"). The function returns the line number of the first
    data line, allowing for proper parsing of the file.

    Parameters:
    -----------
    my_gff : list of str
        A list of strings where each string represents a line from the GFF or GTF file.
        The file may contain header lines (lines starting with '#') before the actual data.

    Returns:
    --------
    int
        The line number (0-based index) of the first data line in the GFF/GTF file.
        This is the line after the header lines, or `0` if no headers are present.
    """
    l_number = 0
    for line in my_gff:
        if line[0] != "#":
            break
        else:
            l_number += 1
    return l_number


def extract_gff_elts(link_to_gff: str) -> (dict, dict):
    """
    Extracts exon and intron information from a GFF file and returns dictionaries containing their positions.

    This function processes a GFF file to extract data about exons and introns. The exons are grouped by chrom/scaffold,
    and the introns are calculated based on the gaps between successive exons. The function returns two dictionaries:
    one for exons and one for introns. Each dictionary contains the chromosome as the key, and a list of lists with
    the start, stop, and orientation for each exon or intron.

    Parameters:
    -----------
    link_to_gff : str
        The path to the GFF file that contains genomic feature annotations, including exons and introns.

    Returns:
    --------
    tuple:
        - dico_chrom_exon_pos (dict): A dictionary where the key is the chromosome/scaffold name and the value
          is a list of exons in the format [start, stop, orientation].
        - dico_chrom_intron_pos (dict): A dictionary where the key is the chromosome/scaffold name and the value
          is a list of introns in the format [start, stop, orientation].

    Example:
    --------
    input = "path/to/file.gff"
    exons, introns = extract_gff_elts(input)
    print(exons)
    # Output: {'chr1': [[100, 200, '+'], [300, 400, '-']]}
    print(introns)
    # Output: {'chr1': [[201, 299, '+']]}
    """
    # generate lists
    list_intron_data = []
    list_exon_data = []
    my_gff = openFile(link_to_gff)
    list_exons_my_gene = []

    # exon markers determines if the exons are called cds or exons in the gff file
    gene_marker, exon_marker = assess_gff_markers(
        my_gff
    )  # i do not use gene marker anymore.
    # goes through all lines of the gtf
    for line in my_gff:
        if line[0] != "#":
            elts = line.split()
            if len(elts) >= 2:
                # ones it detects an exon, it take the list of all next exons in the gene
                # if then the exon list is over, it takes all introns between this succession of introns.
                # The list of introns and exon is reinitialised each time the exon succession is stoped by an other element.

                # if we do not detect an exon
                if elts[2] != exon_marker:  # before was gene marker.
                    # if before that exons were detected
                    if len(list_exons_my_gene) > 1:
                        # we extract introns between the exons from the previous gene/mrna/transcript/whatever it is.
                        list_my_gene_introns = extract_introns(list_exons_my_gene)
                        for intron_list in list_my_gene_introns:
                            list_intron_data.append(intron_list)
                    # we empty the exon list if we meet an other element than an exon
                    list_exons_my_gene = []
                # if there is an exon
                elif elts[2] == exon_marker:
                    # the exon object is built ( a list : [chrom_scaffold, start, stop, orientation])
                    sublist_exon = build_object_ID(elts)
                    # the exon list is implemented in list_exon_data with all exons and list_exons_my_gene for all exons of this specific gene/mRNA/transcripts/whatever it is.
                    list_exon_data.append(sublist_exon)
                    list_exons_my_gene.append(sublist_exon)
    # take care of the last exon set
    if len(list_exons_my_gene) > 1:
        list_my_gene_introns = extract_introns(list_exons_my_gene)
        for intron_list in list_my_gene_introns:
            list_intron_data.append(intron_list)
    dico_chrom_exon_pos = dico_object_per_chrom(list_exon_data)
    dico_chrom_intron_pos = dico_object_per_chrom(list_intron_data)
    return dico_chrom_exon_pos, dico_chrom_intron_pos


def my_bedtool_genes_introns_ordered(link_to_query_genome_gff: str) -> (dict, dict):
    """
    Creates dictionaries of all exons and introns from a GFF file, sorted by their start position.

    This function processes a GFF file containing gene annotations to extract information about exons and introns.
    The exons and introns are extracted into two separate dictionaries, where each key is a chromosome (or scaffold)
    and the value is a list of start, stop, and orientation for each feature. The function then sorts these features
    by their start position in increasing order within each chromosome.

    Parameters:
    -----------
    link_to_query_genome_gff : str
        The path to the GFF file containing genomic feature annotations, including exons and introns.

    Returns:
    --------
    tuple:
        - dico_sorted_chrom_exon_pos (dict): A dictionary where the key is the chromosome/scaffold name and
          the value is a sorted list of exons in the format [start, stop, orientation].
        - dico_sorted_chrom_intron_pos (dict): A dictionary where the key is the chromosome/scaffold name and
          the value is a sorted list of introns in the format [start, stop, orientation].

    Example:
    --------
    input = "path/to/query_genome.gff"
    sorted_exons, sorted_introns = my_bedtool_genes_introns_ordered(input)
    print(sorted_exons)
    # Output: {'chr1': [[100, 200, '+'], [300, 400, '-']]}
    print(sorted_introns)
    # Output: {'chr1': [[201, 299, '+']]}
    """
    # create 2 dictionaries with exons and introns coordinates and chro and orientation ({chrom1 : [[startexon1, stopexon1, orientexon1], [startexon2, stopexon2, orientexon2] ...]})
    dico_chrom_exon_pos, dico_chrom_intron_pos = extract_gff_elts(
        link_to_query_genome_gff
    )
    # sort the list by start position increasing in each chrom list of the dictionaries
    dico_sorted_chrom_exon_pos = sort_gff_dic(
        dico_chrom_exon_pos
    )  # from module_sort_gff_gene_order
    dico_sorted_chrom_intron_pos = sort_gff_dic(
        dico_chrom_intron_pos
    )  # from module_sort_gff_gene_order
    return dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos


def find_intronic_overlap(
    dico_intron: dict,
    dict_transcripts_properties: dict,
    unknown_transcripts: list,
    antisense_transcript: list,
    sign_dot: str,
) -> (list, list):
    """
    This function studies the overlap between transcripts and introns, categorizing them as intronic or antisense based on their position and orientation.

    Transcripts that are located inside an intron in the orientation of transcription are classified as intronic.
    Transcripts that are inside an intron but in the reverse orientation of transcription are classified as antisense.

    Parameters:
    -----------
    dico_intron : dict
        A dictionary where the key is the chromosome name and the value is a list of intron coordinates.
        Each intron is represented by a list: [start, stop, orientation].

    dict_transcripts_properties : dict
        A dictionary containing the properties of transcripts. Each key is a transcript ID, and the value is a list
        with the following information: [chromosome, orientation, start, stop].

    unknown_transcripts : list
        A list of transcript IDs to be checked for overlap with introns.

    antisense_transcript : list
        A list that will be updated with transcript IDs that are classified as antisense (located in reverse orientation
        of transcription relative to the intron).

    sign_dot : str
        A string representing the orientation when a transcript has a dot (`"."`) in the `dict_transcripts_properties` value,
        which is used to define its transcription direction.

    Returns:
    --------
    tuple:
        - list_intronic_transcripts (list): A list of transcript IDs that are intronic (located inside an intron in the
          orientation of transcription).
        - antisense_transcript (list): The updated list of transcript IDs that are classified as antisense (located inside
          an intron but in the reverse orientation of transcription).
    """
    list_intronic_transcripts = []
    for transcript in unknown_transcripts:
        if transcript in dict_transcripts_properties:
            chrom_transcript = dict_transcripts_properties[transcript][0]
            orientation_transcript = dict_transcripts_properties[transcript][1]
            if orientation_transcript == ".":
                orientation_transcript = sign_dot
            start_transcript = int(dict_transcripts_properties[transcript][2])
            stop_transcript = int(dict_transcripts_properties[transcript][3])
            if chrom_transcript in dico_intron:
                list_intron_chrom = dico_intron[chrom_transcript]
                intronic = False
                antisense = False
                for intron_spe_list in list_intron_chrom:
                    if (
                        intron_spe_list[0] <= start_transcript
                        and intron_spe_list[1] >= stop_transcript
                    ):
                        if orientation_transcript == intron_spe_list[2]:
                            intronic = list_intronic_transcripts.append(transcript)
                        else:
                            antisense = antisense_transcript.append(transcript)
                        break
    return list_intronic_transcripts, antisense_transcript


def find_genic_overlap(
    dico_exon,
    dict_transcripts_properties: dict,
    list_accepted_transcript: list,
    sign_dot: str,
) -> (list, list):
    """
    This function detects the overlap between transcripts and exons. Transcripts that overlap with an exon or part of an exon
    in the same orientation as the exon are considered genic, while the same overlap in the reverse orientation (antisense)
    makes the transcript antisense.

    Transcripts that overlap with multiple exons or genes in either sense or antisense orientation are still classified as genic.
    Transcripts that do not overlap with any exons are not classified.

    Parameters:
    -----------
    dico_exon : dict
        A dictionary where the key is the chromosome name, and the value is a list of exons for that chromosome. Each exon
        is represented by a list: [start, stop, orientation].

    dict_transcripts_properties : dict
        A dictionary containing the properties of transcripts. Each key is a transcript ID, and the value is a list of properties:
        [chromosome, orientation, start, stop].

    list_accepted_transcript : list
        A list of transcript IDs that are to be checked for overlap with exons.

    sign_dot : str
        A string representing the orientation when a transcript has a dot (`"."`) in the `dict_transcripts_properties` value,
        used to define its transcription direction.

    Returns:
    --------
    tuple:
        - list_genic (list): A list of transcript IDs that are classified as genic, meaning they overlap with an exon or part
          of an exon in the same orientation.
        - list_antisense (list): A list of transcript IDs that are classified as antisense, meaning they overlap with an exon
          or part of an exon in the reverse orientation.
    """
    list_genic = []
    list_antisense = []
    for transcript in list_accepted_transcript:
        if transcript in dict_transcripts_properties:
            chrom_transcript = dict_transcripts_properties[transcript][0]
            orientation_transcript = dict_transcripts_properties[transcript][1]
            if orientation_transcript == ".":
                orientation_transcript = sign_dot
            start_transcript = int(dict_transcripts_properties[transcript][2])
            stop_transcript = int(dict_transcripts_properties[transcript][3])
            if chrom_transcript in dico_exon:
                list_exon_chrom = dico_exon[chrom_transcript]
                genic = False
                antisense = False
                for exon_spe_list in list_exon_chrom:
                    if (
                        start_transcript <= exon_spe_list[0]
                        and stop_transcript > exon_spe_list[0]
                    ):
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                    if (
                        start_transcript < exon_spe_list[1]
                        and stop_transcript >= exon_spe_list[1]
                    ):
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                    if (
                        start_transcript >= exon_spe_list[0]
                        and stop_transcript <= exon_spe_list[1]
                    ):
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True

                    # New condition to account for the 1 nt overlaps
                    if stop_transcript == exon_spe_list[0]:
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                    if start_transcript == exon_spe_list[1]:  # Change2 is here
                        if orientation_transcript == exon_spe_list[2]:
                            genic = True
                        else:
                            antisense = True
                if genic == True:
                    list_genic.append(transcript)
                else:
                    if antisense == True:
                        list_antisense.append(transcript)
    return list_genic, list_antisense


def my_bedtool_extract_overlap(
    dict_transcripts_properties: dict,
    dict_gene_transcripts: dict,
    dico_sorted_chrom_exon_pos: dict,
    dico_sorted_chrom_intron_pos: dict,
    name_intermediate_directory: str,
) -> list:
    """
    This function explores the overlap between transcripts and genic elements (exons and introns) to classify the transcripts based on their genomic locations.
    It returns lists of transcripts categorized as intergenic, genic, antisense, or intronic based on their overlap with exons and introns.

    The function first loads and processes information about accepted transcripts (based on a TPM threshold). It then identifies the type of each transcript
    by checking whether it overlaps with exons (genic), is on the opposite strand of exons (antisense), or overlaps with introns (intronic). Transcripts that
    do not overlap with exons or introns are considered intergenic.

    Parameters:
    -----------
    dict_transcripts_properties : dict
        A dictionary containing the properties of each transcript. The key is the transcript ID, and the value is a list containing the following:
        [chromosome, orientation, start, stop]. This information is used to determine the location and direction of each transcript.

    dict_gene_transcripts : dict
        A dictionary containing genes and their associated transcripts. The key is the gene ID, and the value is a list of transcript IDs that are associated
        with that gene. Only transcripts above a specified TPM threshold are included in this list.

    dico_sorted_chrom_exon_pos : dict
        A dictionary containing the positions of exons on each chromosome. The key is the chromosome name, and the value is a list of exons represented as
        [start, stop, orientation].

    dico_sorted_chrom_intron_pos : dict
        A dictionary containing the positions of introns on each chromosome. The key is the chromosome name, and the value is a list of introns represented as
        [start, stop, orientation].

    name_intermediate_directory : str
        The path to the directory where intermediate files (such as orientation data) are stored. This is used to determine the orientation of transcripts.

    Returns:
    --------
    list:
        A list of four lists:
        - intergenic_transcripts (list): A list of transcript IDs that are located in intergenic regions (not overlapping with exons or introns).
        - intronic_transcripts (list): A list of transcript IDs that overlap with introns in the gene frame.
        - antisense_transcript (list): A list of transcript IDs that overlap with exons or introns in the reverse orientation (antisense).
        - genic_transcripts (list): A list of transcript IDs that overlap with exons in the same orientation (genic).
    """
    # Initialize a list to store all transcripts that will be used, based on if they are present in the dict_gene_transcripts, which contain only the transcripts validated previously with a specific TPM value and putatively TE.
    list_accepted_transcript = []
    link_file_dot = name_intermediate_directory + "/Dot_indication.txt"
    file_orientation_dot_transc = openFile(link_file_dot)
    sign_dot = file_orientation_dot_transc[0].split("\n")[0]
    # Iterate through genes and their transcripts to build the list of accepted transcripts
    for gene in dict_gene_transcripts:
        for transcript in dict_gene_transcripts[gene]:
            list_accepted_transcript.append(transcript)

    # Initialize lists for different types of transcripts
    intergenic_transcripts = []
    genic_transcripts = []  # These one ovrlap with gens in the same direction and are not inside an intron
    antisense_transcript = []
    intronic_transcripts = []
    print(
        BRIGHT_BLUE
        + "Total number of transcripts above the specified TPM thresshold : "
        + RESET
        + str(len(list_accepted_transcript))
    )
    genic_transcripts, antisense_transcript = find_genic_overlap(
        dico_sorted_chrom_exon_pos,
        dict_transcripts_properties,
        list_accepted_transcript,
        sign_dot,
    )
    unknown_transcripts = [
        x
        for x in list_accepted_transcript
        if x not in genic_transcripts and x not in antisense_transcript
    ]
    intronic_transcripts, antisense_transcript = find_intronic_overlap(
        dico_sorted_chrom_intron_pos,
        dict_transcripts_properties,
        unknown_transcripts,
        antisense_transcript,
        sign_dot,
    )
    intergenic_transcripts = [
        x
        for x in unknown_transcripts
        if x not in antisense_transcript and x not in intronic_transcripts
    ]
    print(
        BRIGHT_BLUE
        + "Number of genic transcripts (not denovo) : "
        + RESET
        + str(len(genic_transcripts))
    )
    print(
        BRIGHT_BLUE
        + "Number of intergenic transcripts : "
        + RESET
        + str(len(intergenic_transcripts))
    )
    print(
        BRIGHT_BLUE
        + "Number of antisense transcripts : "
        + RESET
        + str(len(antisense_transcript))
    )
    print(
        BRIGHT_BLUE
        + "Number of intronic transcripts (in gene frame) : "
        + RESET
        + str(len(intronic_transcripts))
    )
    return [
        intergenic_transcripts,
        intronic_transcripts,
        antisense_transcript,
        genic_transcripts,
    ]


def index_choice_transcript_overlap(list_choice_denovo: list) -> list:
    """
    This function maps user choices for transcript overlap types to corresponding indices.
    It takes a list of user choices and returns a list of indices based on the presence of specific choices.

    Parameters:
    -----------
    list_choice_denovo : list
        A list of strings where each string represents a type of transcript overlap the user is interested in.
        The possible values are "intergenic", "intronic", "antisense", and "genic".

    Returns:
    --------
    list:
        A list of integers representing the indices for each selected transcript overlap type.
        The mapping is as follows:
        - "intergenic" -> index 0
        - "intronic" -> index 1
        - "antisense" -> index 2
        - "genic" -> index 3
    """
    # Initialize an empty list to store indices
    list_index = []

    # Check each type of transcript overlap and append the corresponding index if present in the user's choices
    if "intergenic" in list_choice_denovo:
        list_index.append(0)
    if "intronic" in list_choice_denovo:
        list_index.append(1)
    if "antisense" in list_choice_denovo:
        list_index.append(2)
    if "genic" in list_choice_denovo:
        list_index.append(3)

    # Return the list of indices
    return list_index


def extract_denovo_transcripts(
    sublist_transcripts_overlap: list,
    list_choice_denovo: list,
    dict_transcript_fasta: dict,
    dict_gene_transcripts: dict,
    name_output_directory: str,
) -> (dict, dict, dict, dict):
    """
    Extracts denovo transcripts based on user-defined overlap types and generates output files and dictionaries.

    This function processes the provided lists of transcripts that overlap specific genomic features (such as intergenic, intronic, antisense, or genic).
    It organizes the transcripts based on user choices, writes the corresponding FASTA sequences to a file,
    and creates dictionaries mapping gene names to their associated denovo transcripts and other relevant data.

    Parameters:
    ----------
    sublist_transcripts_overlap : list of lists
        A list of lists, where each sublist contains transcripts that correspond to a specific overlap type (e.g., intergenic, intronic).
        Each list represents transcripts that overlap a specific genomic feature type.

    list_choice_denovo : list of str
        A list of strings where each string represents a chosen overlap type (such as "intergenic", "intronic", "antisense", "genic")
        indicating which types of overlaps the user is interested in for denovo transcript extraction.

    dict_transcript_fasta : dict
        A dictionary where keys are transcript names and values are the corresponding FASTA sequences of the transcripts.

    dict_gene_transcripts : dict
        A dictionary where keys are gene names and values are lists of transcript names that belong to that gene.
        It contains the reference or existing transcripts for each gene.

    name_output_directory : str
        The directory path where the output files, specifically the FASTA file of denovo transcripts, will be saved.

    Returns:
    --------
    dict_denovo_gene_transcripts : dict
        A dictionary where keys are gene names and values are lists of transcript names that are considered denovo for each gene.

    dict_denovo_transcript_fasta : dict
        A dictionary where keys are denovo transcript names and values are their corresponding FASTA sequences.

    dict_transcrip_status : dict
        A dictionary where keys are transcript names and values are the status of the transcript (e.g., "intergenic", "intronic", "antisense", "genic").

    dict_gene_status : dict
        A dictionary where keys are gene names and values are the status of the gene, which is either "genic" (existing transcripts) or "denovo" (newly detected transcripts).

    Notes:
    -----
    - The FASTA file of denovo transcripts is saved in the specified `name_output_directory` as `denovo_transcript.fa`.
    - Genes with a greater number of denovo transcripts compared to their known transcripts are labeled as "denovo".
    - The `list_choice_denovo` parameter determines which types of transcript overlaps to include in the denovo analysis.
    """
    # Initialize dictionaries to store denovo gene transcripts, denovo transcript FASTA sequences, transcript status, and gene status
    dict_denovo_gene_transcripts = {}
    dict_denovo_transcript_fasta = {}
    dict_transcrip_status = {}
    dict_gene_status = {}

    # Open a FASTA file to write denovo transcripts
    link_file_denovo_transcripts = name_output_directory + "/denovo_transcript.fa"
    fasta_denovo = open(link_file_denovo_transcripts, "w")

    # Map user choices to corresponding indices
    dict_index_value = {0: "intergenic", 1: "intronic", 2: "antisense", 3: "genic"}
    list_indexe_overlap = index_choice_transcript_overlap(list_choice_denovo)

    # Iterate through each chosen transcript overlap type
    for index in list_indexe_overlap:
        my_transcript_list = sublist_transcripts_overlap[index]

        # Iterate through transcripts in the chosen overlap type
        for transcript_name in my_transcript_list:
            gene_name = (
                transcript_name.split(".")[0] + "." + transcript_name.split(".")[1]
            )

            # Update dictionaries with denovo gene transcripts and denovo transcript FASTA sequences
            if gene_name in dict_denovo_gene_transcripts:
                dict_denovo_gene_transcripts[gene_name].append(transcript_name)
            else:
                dict_denovo_gene_transcripts[gene_name] = [transcript_name]

            dict_denovo_transcript_fasta[transcript_name] = dict_transcript_fasta[
                transcript_name
            ]
            dict_transcrip_status[transcript_name] = dict_index_value[index]

            # Write denovo transcripts to the FASTA file
            fasta_denovo.write(">" + transcript_name + "\n")
            fasta_denovo.write(dict_denovo_transcript_fasta[transcript_name] + "\n")

    # Close the FASTA file
    fasta_denovo.close()

    # Update gene status based on the number of denovo transcripts
    for gene in dict_denovo_gene_transcripts:
        if len(dict_denovo_gene_transcripts[gene]) < len(dict_gene_transcripts[gene]):
            dict_gene_status[gene] = "genic"
        else:
            dict_gene_status[gene] = "denovo"

    # Return dictionaries containing denovo gene transcripts, denovo transcript FASTA sequences, transcript status, and gene status
    return (
        dict_denovo_gene_transcripts,
        dict_denovo_transcript_fasta,
        dict_transcrip_status,
        dict_gene_status,
    )
