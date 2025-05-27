import os
from Bio import SeqIO
from deswoman.module_orfs import build_dict_orfs_per_transcript
from deswoman.module_colors import openFile


def build_dict_unspliced_seqs(opened_gtf_file: list) -> dict:
    """
    This function parses a GTF file containing transcriptome assembly information.
    It creates a dictionary where each transcript name is associated with a list containing
    chromosome, direction of transcription, and a sublist with exon start and stop positions.

    Parameters:
    - opened_gtf_file: The opened GTF file for transcriptome assembly.

    Returns:
    A dictionary where keys are transcript names, and values are lists containing chromosome,
    direction of transcription, and a sublist with exon start and stop positions.
    """

    # Initialize an empty dictionary to store transcript information
    dict_gene_exon_pos = {}
    transcript_name = ""
    # Iterate through each line in the opened GTF file
    for line in opened_gtf_file:
        if line[0] != "#":
            # Split the line into elements
            list_elts_line = line.split("	")

            # Check if the line represents a transcript entry
            if list_elts_line[2] == "transcript":
                # If a previous transcript was being processed, store its information in the dictionary
                if transcript_name != "":
                    transcript_info = [chrom, direction, list_exons]
                    dict_gene_exon_pos[transcript_name] = transcript_info

                # Reset variables for the new transcript
                list_exons = []
                chrom = list_elts_line[0]
                direction = list_elts_line[6]

                # Extract transcript name from the line
                subLigne = list_elts_line[8].split(";")
                subsubLigne = subLigne[1].split(" ")
                transcript_name = subsubLigne[2][1 : len(subsubLigne[2]) - 1]

            elif list_elts_line[2] == "exon":
                # Extract start and end positions for exons
                start_exon = int(list_elts_line[3])
                end_exon = int(list_elts_line[4])
                subList_exon = [start_exon, end_exon]
                list_exons.append(subList_exon)

    # Store the information for the last transcript in the dictionary
    transcript_info = [chrom, direction, list_exons]
    dict_gene_exon_pos[transcript_name] = transcript_info
    return dict_gene_exon_pos


def build_dict_chrom(genome_fasta: str) -> dict:
    """
    Parses a genome FASTA file and builds a dictionary of chromosomes.

    This function reads a FASTA file containing genome sequences and creates a dictionary where each key is the chromosome
    name (ID) and each value is the corresponding nucleotide sequence.

    Args:
        genome_fasta (file object): The opened FASTA file containing the genome sequences.

    Returns:
        dict: A dictionary where each key is the chromosome name (ID) and each value is the chromosome's nucleotide sequence.
    """
    dict_chrom = {}
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        dict_chrom[str(seq_record.id)] = str(seq_record.seq)
    return dict_chrom


def reverse_sequence(transcript: str) -> str:
    """
    Reverse transcribes a nucleotide sequence by replacing each base with its complementary base.

    This function takes a nucleotide sequence (in FASTA format) and returns its reverse complement
    by swapping each nucleotide with its complement (A <-> T, G <-> C, etc.), including handling
    both uppercase and lowercase letters. The reverse sequence is returned as a string.

    Args:
        transcript (str): A nucleotide sequence represented as a string, typically in FASTA format.

    Returns:
        str: The reverse complement of the input nucleotide sequence.
    """
    reverse = ""
    Dico = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "n": "n",
        "y": "y",
        "Y": "Y",
        "U": "U",
        "u": "u",
        "R": "R",
        "r": "r",
        "K": "K",
        "k": "k",
        "M": "M",
        "m": "m",
        "S": "S",
        "s": "s",
        "W": "W",
        "w": "w",
        "B": "B",
        "b": "b",
        "D": "D",
        "d": "d",
        "H": "H",
        "h": "h",
        "V": "V",
        "v": "v",
    }
    L = list(transcript.strip())
    for i in reversed(L):
        new_nucl = Dico[i]
        reverse += new_nucl
    return reverse


def build_new_spliced_transcript(
    list_transcript_elts: list,
    dict_chromosomes_fasta: dict,
    direction_transcription: str,
) -> (str, dict):
    """
    Builds a new spliced transcript sequence by extracting exons from the genome sequence and extending
    the transcript slightly in both directions, considering the transcription direction.

    Parameters:
    - list_transcript_elts (list): A list containing properties of the transcript:
        0: Name of the chromosome where the transcript is located.
        1: Direction of transcription (either '+' or '-').
        2: List of exons with their start and stop positions in the genome.
    - dict_chromosomes_fasta (dict): A dictionary containing the genome sequences in FASTA format,
      with chromosomes as keys and their sequences as values.
    - direction_transcription (str): The direction of transcription ('+' for forward, '-' for reverse).

    Returns:
    - tuple:
        - str: The newly spliced transcript sequence, extended slightly in both the 5' and 3' directions.
        - dict: A dictionary mapping the position in the spliced transcript to the corresponding position in the genome.
    """
    chrom_name = list_transcript_elts[0]
    list_exons_pos = list_transcript_elts[2]
    chrom_sequence = dict_chromosomes_fasta[chrom_name]

    c = 3  # Counter for position in the spliced transcript
    dict_pos_spliced_transcript = {
        1: list_exons_pos[0][0] - 3,
        2: list_exons_pos[0][0] - 2,
    }

    if direction_transcription == "+":
        new_spliced_transcript = ""

        # Extend in the 5' direction
        new_spliced_transcript += chrom_sequence[
            list_exons_pos[0][0] - 3 : list_exons_pos[0][0] - 1
        ]

        # Add exonic sequences
        for sublist_exon in list_exons_pos:
            my_exon_seq = chrom_sequence[sublist_exon[0] - 1 : sublist_exon[1]]
            for j in range(sublist_exon[0] - 1, sublist_exon[1]):
                dict_pos_spliced_transcript[c] = j
                c += 1
            new_spliced_transcript += my_exon_seq

        # Extend in the 3' direction
        if len(chrom_sequence) > list_exons_pos[len(list_exons_pos) - 1][1]:
            new_spliced_transcript += chrom_sequence[
                list_exons_pos[len(list_exons_pos) - 1][1]
            ]
        dict_pos_spliced_transcript[c] = list_exons_pos[len(list_exons_pos) - 1][1]

    elif direction_transcription == "-":  # Decouverte donc, les "." sont des reverse
        new_spliced_transcript = ""

        # Extend in the 5' direction
        new_spliced_transcript += chrom_sequence[
            list_exons_pos[0][0] - 3 : list_exons_pos[0][0] - 1
        ]

        # Add exonic sequences
        for sublist_exon in list_exons_pos:
            my_exon_seq = chrom_sequence[sublist_exon[0] - 1 : sublist_exon[1]]
            for j in range(sublist_exon[0] - 1, sublist_exon[1]):
                dict_pos_spliced_transcript[c] = j
                c += 1
            new_spliced_transcript += my_exon_seq

        # Extend in the 3' direction
        if len(chrom_sequence) > list_exons_pos[len(list_exons_pos) - 1][1]:
            new_spliced_transcript += chrom_sequence[
                list_exons_pos[len(list_exons_pos) - 1][1]
            ]
        dict_pos_spliced_transcript[c] = list_exons_pos[len(list_exons_pos) - 1][1]

        # Reverse the transcript if direction is "-" or "."
        new_spliced_transcript = reverse_sequence(new_spliced_transcript)

    return new_spliced_transcript, dict_pos_spliced_transcript


def max_dict(my_dict: dict) -> int:
    """
    This function finds the maximum key in a dictionary based on its key values.

    Parameters:
    - my_dict (dict): A dictionary where the keys are comparable values (e.g., integers, floats).

    Returns:
    - The maximum key in the dictionary. If the dictionary is empty, the function returns 0.
    """
    max_detected = 0
    for i in my_dict:
        if i > max_detected:
            max_detected = i
    return max_detected


def rebuild_spliced_orfs(
    genome_file: str,
    dict_gene_structure: dict,
    dict_orfs: dict,
    name_intermediate_directory: str,
) -> dict:
    """
    Reconstructs spliced transcripts based on GTF information, extracts spliced ORFs, and validates the reconstructed ORFs.
    If valid, it retrieves the stop codon of the spliced ORF in the genome and the coordinates of the unspliced ORF.
    The function returns a dictionary with validated ORFs and their coordinates in the genome.

    Parameters:
    - genome_file (str): Path to the genome FASTA file containing the genome sequences.
    - dict_gene_structure (dict): A dictionary where keys are transcript names, and values are lists containing
                                  transcript properties such as chromosome, orientation, and exons.
    - dict_orfs (dict): A dictionary where keys are ORF names, and values are the sequences of the corresponding ORFs.
    - name_intermediate_directory (str): Directory path where intermediate files (such as a dot file indicating transcript orientation) are stored.

    Returns:
    - dict_orf_unspliced_start_end_stop (dict): A dictionary where the keys are ORF names, and the values are lists
                                                 containing the start coordinate, stop coordinate (including the stop codon),
                                                 and stop codon sequence for each validated ORF.

    Notes:
    - The function checks that the reconstructed spliced ORF matches the original ORF.
    - It also verifies that the ORF contains a valid stop codon from the list ["TAG", "TAA", "TGA"].
    - The function handles both positive and negative transcription orientations by reconstructing the spliced transcripts accordingly.
    - The dictionary returned contains ORFs that passed the validation criteria for stop codon presence and sequence consistency.
    """
    dict_orf_unspliced_start_end_stop = {}
    dict_chromosomes_fasta = build_dict_chrom(genome_file)
    dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_orfs)
    stop_codon_list = ["TAG", "TAA", "TGA"]
    compteur_problematic_seq = 0
    dic_nucl_translation = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "n": "n",
        "y": "y",
        "Y": "Y",
        "U": "U",
        "u": "u",
        "R": "R",
        "r": "r",
        "K": "K",
        "k": "k",
        "M": "M",
        "m": "m",
        "S": "S",
        "s": "s",
        "W": "W",
        "w": "w",
        "B": "B",
        "b": "b",
        "D": "D",
        "d": "d",
        "H": "H",
        "h": "h",
        "V": "V",
        "v": "v",
    }
    total = 0
    dot_file_link = name_intermediate_directory + "/Dot_indication.txt"
    file_orientation_dot_transc = openFile(dot_file_link)
    # this next line to acces the direction of transcripts with a dot
    sign_dot = file_orientation_dot_transc[0].split("\n")[0]
    for transcript_name in dict_gene_structure:
        if (
            transcript_name in dict_transcrit_all_orfs
        ):  # remove transcripts that do not have an ORF
            list_transcript_elts = dict_gene_structure[transcript_name]
            # sequence_to_be_taken = True

            for orf_name in dict_transcrit_all_orfs[transcript_name]:
                sequence_to_be_taken = True
                # The next block extracts the ORFs and rebuild the spliced sequences based on the informations of the unspliced reconstructed transcritps
                if orf_name in dict_orfs:
                    total += 1
                    list_elts_orf_name = orf_name.split("_")
                    start_orf = int(list_elts_orf_name[2])
                    end_orf = int(list_elts_orf_name[3])
                    spliced_orf = dict_orfs[orf_name]
                    orientation_transcript = list_transcript_elts[1]
                    if orientation_transcript == ".":
                        orientation_transcript = sign_dot

                    if orientation_transcript == "+":
                        # reconstruct spliced transcript from the genome based on the exon intron structure of the transcript contrinaed in the gtf
                        (
                            reconstructed_spliced_transcript,
                            dict_pos_spliced_transcript,
                        ) = build_new_spliced_transcript(
                            list_transcript_elts, dict_chromosomes_fasta, "+"
                        )
                        reconstructed_spliced_orf = reconstructed_spliced_transcript[
                            start_orf + 1 : end_orf + 2
                        ]
                        chromosome_fasta = dict_chromosomes_fasta[
                            list_transcript_elts[0]
                        ]
                        start_orf_in_genome = dict_pos_spliced_transcript[start_orf + 2]
                        end_orf_in_genome = dict_pos_spliced_transcript[end_orf + 2]
                        # The reason to the fact that the ORFs are added the stop codon is that the original software getORF did not extract the stop codons from the sequences
                        # THe getORF functions that we built does then the same

                        end_orf_in_genome_with_stop = (
                            dict_pos_spliced_transcript[end_orf + 2] + 3
                        )
                        orf_stop_codon = chromosome_fasta[
                            end_orf_in_genome + 1 : end_orf_in_genome + 4
                        ].upper()

                        if orf_stop_codon not in stop_codon_list:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False

                        if spliced_orf != reconstructed_spliced_orf:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False

                    elif (
                        orientation_transcript == "-"
                    ):  # Decouverte donc, les "." sont des forward ou des reverse selon les transcripts
                        # reconstruct spliced transcript from the genome based on the exon intron structure of the transcript contrinaed in the gtf
                        (
                            reconstructed_spliced_transcript,
                            dict_pos_spliced_transcript,
                        ) = build_new_spliced_transcript(
                            list_transcript_elts, dict_chromosomes_fasta, "-"
                        )
                        reconstructed_spliced_orf = reconstructed_spliced_transcript[
                            start_orf : end_orf + 1
                        ]
                        chromosome_fasta = dict_chromosomes_fasta[
                            list_transcript_elts[0]
                        ]
                        start_orf_in_genome = dict_pos_spliced_transcript[
                            max_dict(dict_pos_spliced_transcript) - start_orf
                        ]
                        end_orf_in_genome = dict_pos_spliced_transcript[
                            max_dict(dict_pos_spliced_transcript) - end_orf
                        ]
                        end_orf_in_genome_with_stop = end_orf_in_genome - 3
                        orf_stop_codon = reverse_sequence(
                            chromosome_fasta[
                                end_orf_in_genome_with_stop:end_orf_in_genome
                            ]
                        ).upper()

                        if spliced_orf != reconstructed_spliced_orf:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False
                        if orf_stop_codon not in stop_codon_list:
                            compteur_problematic_seq += 1
                            sequence_to_be_taken = False
                    if sequence_to_be_taken is True:
                        dict_orf_unspliced_start_end_stop[orf_name] = [
                            start_orf_in_genome,
                            end_orf_in_genome_with_stop,
                            orf_stop_codon,
                        ]
    return dict_orf_unspliced_start_end_stop


def add_stop_codon_to_spliced_orfs(
    dict_pos_unspliced_ORFs_plus_stop: dict, dict_all_ORFs: dict
) -> dict:
    """
    Adds the corresponding stop codon to the spliced ORFs that passed the stop codon and identity test.
    This function creates a dictionary of validated ORFs, appending the stop codon to the ORF sequences.

    Parameters:
    - dict_pos_unspliced_ORFs_plus_stop (dict): A dictionary where keys are ORF names and values are lists
                                                 containing the ORF start position, stop position, and stop codon sequence.
    - dict_all_ORFs (dict): A dictionary where keys are ORF names and values are the corresponding ORF sequences in FASTA format.

    Returns:
    - dict_all_orfs_with_stop (dict): A dictionary where keys are ORF names and values are ORF sequences with the
                                      appended stop codon.

    Notes:
    - This function ensures that each ORF sequence is properly terminated with the correct stop codon as validated earlier.
    """
    dict_all_orfs_with_stop = {}
    for orf_name in dict_pos_unspliced_ORFs_plus_stop:
        orf_seq = dict_all_ORFs[orf_name]
        orf_seq_with_stop = orf_seq + dict_pos_unspliced_ORFs_plus_stop[orf_name][2]
        dict_all_orfs_with_stop[orf_name] = orf_seq_with_stop

    return dict_all_orfs_with_stop


def get_dot_orientation(
    genome_file: str,
    dict_transcripts_properties: dict,
    dict_transcript_fasta_all: dict,
    gtf_file: str,
    name_intermediate_directory: str,
) -> None:
    """
    Determines the orientation of transcripts with a dot (".") based on reconstructed spliced sequences.

    This function reconstructs spliced transcripts from the provided genome and GTF files and compares them to the
    transcript sequences from a transcriptome assembly. If the reconstructed transcript matches the transcript from
    the assembly, the orientation of the dot (".") is set to "+". Otherwise, it is set to "-".

    The function analyzes all transcripts with an unassigned orientation (".") and determines the majority orientation
    among those. The final orientation is written to a file in the intermediate directory.

    Parameters:
    - genome_file (str): Path to the genome file in FASTA format.
    - dict_transcripts_properties (dict): A dictionary where keys are transcript names and values are lists containing
                                          transcript properties, including orientation.
    - dict_transcript_fasta_all (dict): A dictionary where keys are transcript names and values are the corresponding
                                        transcript sequences from transcriptome assembly.
    - gtf_file (str): Path to the GTF file containing transcript information.
    - name_intermediate_directory (str): Path to the directory where intermediate files are stored.

    Returns:
    - None: The function writes the determined orientation (either "+" or "-") to a file in the intermediate directory.

    Notes:
    - If there are no transcripts with an unassigned orientation, the dot (".") is arbitrarily set to "+".
    - The majority orientation is determined based on the number of correctly reconstructed spliced transcripts.
    """
    ID_sign = "."
    nb_plus = 0
    nb_tot = 0
    # create a dictionnary with transcripts as key, associated with a list containing their chrom, orientation and exons pos.
    dict_gene_exon_pos = build_dict_unspliced_seqs(gtf_file)
    # associate chromosomes name (key) to their fasta sequence.
    dict_chromosomes_fasta = build_dict_chrom(genome_file)
    for transcript_name in dict_gene_exon_pos:
        # if the transcript is unroriented
        if (
            transcript_name in dict_transcripts_properties
            and dict_transcripts_properties[transcript_name][1] == "."
        ):
            # extract the transcript list containing chrom, orientation and exons pos
            list_transcript_elts = dict_gene_exon_pos[transcript_name]
            # here it reconstructs the unoriented spliced transcript in the + position arbitrarily
            reconstructed_spliced_transcript, dict_pos_spliced_transcript = (
                build_new_spliced_transcript(
                    list_transcript_elts, dict_chromosomes_fasta, "+"
                )
            )
            reconstructed_spliced_transcript = reconstructed_spliced_transcript[
                2 : len(reconstructed_spliced_transcript) - 1
            ]
            # if the recontructed transcript in forward position corresponds to the transcripts : nb_plus positif takes 1
            if (
                reconstructed_spliced_transcript
                == dict_transcript_fasta_all[transcript_name]
            ):
                nb_plus += 1
            nb_tot += 1
    # if there is no dot, "." orientation is +. Otherwise, we calculate the number of "." transcripts for which the + splicing worked.
    if nb_tot > 0:
        percentage_correct = 100 * nb_plus / nb_tot
    else:  # if there is no transcript in "." orientation, "." will be considered as "+", arbitrarily, the file will not get read anyway
        percentage_correct = 100
    if percentage_correct > 50.0:
        ID_sign = "+"
    else:
        ID_sign = "-"
    # we attribute the sign "+" or "-" according to the majority.
    dot_file_link = name_intermediate_directory + "/Dot_indication.txt"
    os.system("touch " + dot_file_link)
    dot_indication = open(dot_file_link, "w")
    dot_indication.write(ID_sign)
    dot_indication.close()
