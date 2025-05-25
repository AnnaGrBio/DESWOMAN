import os
from Bio import SeqIO
from deswoman.module_transcripts_data_and_threeshold import get_path_gff
from deswoman.module_colors import openFile


__author__ = "Anna Grandchamp"
__contributor__ = ""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def store_genome_informations(species_name: str, path_genome: str) -> dict:
    """
    Retrieve gene information from the GFF file associated with the genome of the specified species.

    This function reads a GFF (General Feature Format) file corresponding to the genome of a species
    to extract details about genes annotated as 'mRNA'. For each gene, the function retrieves the gene ID,
    chromosome, start, and stop positions. The gene information is stored in a dictionary, with chromosomes as
    keys and lists of gene coordinates and IDs as values.

    Parameters:
    -----------
    species_name : str
        The name of the species. This is used to retrieve the appropriate GFF file for the species.

    path_genome : str
        The path to the directory where the genome files are stored. This is used to construct the full path
        to the species' GFF file.

    Returns:
    --------
    dico_gff : dict
        A dictionary containing gene information organized by chromosome. The keys are chromosome names (str),
        and the values are lists of lists, where each inner list contains the start position (str), stop position
        (str), and gene ID (str) for a gene on that chromosome.
    """
    link_to_gff = get_path_gff(path_genome, species_name)
    my_file = openFile(link_to_gff)
    dico_gff = {}
    for line in my_file:
        list_elts = line.split()
        if len(list_elts) > 2 and list_elts[2] == "mRNA":
            chromosome = list_elts[0]
            sublist = list_elts[8].split(";")[0]
            gene_ID = sublist.split("=")[1]
            start = list_elts[3]
            stop = list_elts[4]
            if chromosome not in dico_gff:
                dico_gff[chromosome] = [[start, stop, gene_ID]]
            else:
                dico_gff[chromosome].append([start, stop, gene_ID])
    return dico_gff


def sort_gff_dic(dico_gff: dict) -> dict:
    """
    Sort the genes within each chromosome in increasing order of their start position.

    This function takes a dictionary containing gene information (organized by chromosome)
    and sorts the genes for each chromosome in ascending order based on their start position.
    If two genes share the same start position and direction, the longer gene is retained.

    Parameters:
    -----------
    dico_gff : dict
        A dictionary where keys are chromosome names (str), and values are lists of lists,
        where each inner list contains the start position (str), stop position (str),
        and gene ID (str) for a gene on that chromosome.

    Returns:
    --------
    new_dico_gff : dict
        A new dictionary with the same structure as `dico_gff`, but with genes sorted in increasing
        order of their start positions for each chromosome. If two genes have the same start position
        and direction, the longest gene is retained.

    Example:
    --------
    dico_gff = {
        "chr1": [["100", "200", "gene1"], ["150", "250", "gene2"], ["50", "100", "gene3"]],
        "chr2": [["300", "400", "gene4"], ["350", "450", "gene5"]]
    }
    sorted_dico_gff = sort_gff_dic(dico_gff)
    # sorted_dico_gff will have genes sorted by their start positions on each chromosome.
    """
    new_dico_gff = {}
    # c = 0
    for chromosome in dico_gff:
        new_dico_gff[chromosome] = []
        # sort only if there is more than 1 elt in the chromosome.
        if len(dico_gff[chromosome]) > 1:
            dico_start_assoc_struct = {}
            list_genes_pos = dico_gff[chromosome]
            list_ordered_start = []
            for elt in list_genes_pos:
                list_ordered_start.append(int(elt[0]))
                # just implemented. If two elements have the same start and same direction, we only keep the longest.
                if (
                    int(elt[0]) in dico_start_assoc_struct
                    and elt[2] == dico_start_assoc_struct[int(elt[0])][2]
                ):
                    if int(elt[1]) > int(dico_start_assoc_struct[int(elt[0])][1]):
                        dico_start_assoc_struct[int(elt[0])] = elt
                else:
                    dico_start_assoc_struct[int(elt[0])] = elt
            list_ordered_start.sort()
            for start in list_ordered_start:
                new_dico_gff[chromosome].append(dico_start_assoc_struct[start])
        else:
            new_dico_gff[chromosome] = dico_gff[chromosome]
    return new_dico_gff
