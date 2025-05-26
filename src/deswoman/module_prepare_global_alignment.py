import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from deswoman.module_colors import openFile


__author__ = "Anna Grandchamp"
__contributor__ = "Marie Lebherz"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def build_dico_seq_denovo_filtered_candidates(
    dico_denovo_best_hit: dict, name_output_directory: str
) -> dict:
    """
    Builds a dictionary containing denovo sequences, including introns for sequences that have them.

    This function creates a dictionary where each key is the ID of a denovo sequence and each value is the corresponding
    sequence string (including introns if available). It looks for sequence IDs that are present in the `dico_denovo_best_hit`
    dictionary and then fetches the corresponding sequences from the file `denovo_unspliced_lowered_introns.fa` located in the
    provided `name_output_directory`.

    Parameters:
    -----------
    dico_denovo_best_hit : dict
        A dictionary where keys are sequence IDs (denovo hits) and values are the best match or hit for that sequence.
        Only sequences that have a match in this dictionary will be included in the output.

    name_output_directory : str
        Path to the output directory where the `denovo_unspliced_lowered_introns.fa` file is located. This file should
        contain denovo sequences, potentially with introns, in FASTA format.

    Returns:
    --------
    dict
        A dictionary (`dico_denovo_seq_intron`) where keys are sequence IDs (from `dico_denovo_best_hit`) and values are
        the corresponding denovo sequence strings (including introns if they are present).

    Example:
    --------
    dico_denovo_best_hit = {"seq1": "best_hit1", "seq2": "best_hit2"}
    name_output_directory = "genome/outgroup"
    denovo_sequences = build_dico_seq(dico_denovo_best_hit, name_output_directory)
    """
    dico_denovo_seq_intron = {}
    link_my_file = name_output_directory + "/denovo_unspliced_lowered_introns.fa"
    for seq_record in SeqIO.parse(link_my_file, "fasta"):
        if str(seq_record.id) in dico_denovo_best_hit:
            dico_denovo_seq_intron[str(seq_record.id)] = str(seq_record.seq)
    return dico_denovo_seq_intron


def build_dico_seq_NcHomologs(
    dico_denovo_best_hit: dict, link_to_genome: str, pop_name: str
) -> dict:
    """
    Extracts the sequences of non-coding homologous hits from a genome, slightly extended upstream and downstream.

    This function creates a dictionary where each key is a denovo sequence ID, and each value is the sequence of its
    corresponding non-coding homolog, extended by 2 bases upstream and downstream. If the homologous hit is located on the
    reverse strand, the extracted sequence is reverse complemented.

    Parameters:
    -----------
    dico_denovo_best_hit : dict
        A dictionary where keys are denovo sequence IDs, and values are strings representing the coordinates of
        the homologous hit in the format "start-stop-chromosome-strand".
        The strand is either "+" (forward) or "r" (reverse).

    link_to_genome : str
        The path to the genome file in FASTA format. This file should contain chromosome sequences, which are
        indexed by chromosome name.

    pop_name : str
        The name of the population (currently not used in the function but can be used for annotation or logging purposes).

    Returns:
    --------
    dict
        A dictionary (`dico_NcHit`) where the keys are denovo sequence IDs, and the values are the corresponding
        extended non-coding homologous sequences. The sequences are either taken as-is or reverse complemented
        depending on the strand of the homologous hit.

    Example:
    --------
    dico_denovo_best_hit = {
        "hit1": "1000-1020-chrom1-+",
        "hit2": "2000-2020-chrom2-r"
    }
    link_to_genome = "genome.fasta"
    pop_name = "pop1"

    non_coding_hits = build_dico_seq_NcHomologs(dico_denovo_best_hit, link_to_genome, pop_name)
    """
    dico_NcHit = {}
    dico_genome = {}
    # Store the genome in FASTA format, with chromosome as keys
    for seq_record in SeqIO.parse(link_to_genome, "fasta"):
        dico_genome[str(seq_record.id)] = str(seq_record.seq)

    # Then we extract the sequence of the homologous hits, based on the genome and the coordinate of the homologous hits
    # If the homologous hit is reverse, it is reverse transcribed
    for denovo in dico_denovo_best_hit:
        hit_coords = dico_denovo_best_hit[denovo]
        chrom_scaf = hit_coords.split("-")[2]
        start = int(hit_coords.split("-")[0]) - 1  # changed by Marie, previous -2
        stop = int(hit_coords.split("-")[1])  # changed by Marie, previous +2
        seq_hit_extended = dico_genome[chrom_scaf][start:stop]
        # next line in case the -2 and +2 adding leads out of the boundaries of the chromosome and gives an empty sequence.
        if len(seq_hit_extended) == 0:
            seq_hit_extended = dico_genome[chrom_scaf][start + 2 : stop - 2]
        if hit_coords.split("-")[3] == "r":
            seq_hit_extended = Seq(seq_hit_extended)
            seq_hit_extended_rev_comp = str(seq_hit_extended.reverse_complement())
            dico_NcHit[denovo] = seq_hit_extended_rev_comp
        else:
            dico_NcHit[denovo] = seq_hit_extended
    return dico_NcHit
