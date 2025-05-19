import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Align # ***
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from module_colors import openFile


__author__ = "Anna Grandchamp"
__contributor__="Marie Lebherz"
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def check_ali(seq : str) -> bool:
    """
    This function checks if the start codon (ATG) and stop codon (TAG, TAA, TGA) are intact in the alignment.
    If both the start and stop codons are not interrupted by gaps ('-'), the alignment is considered valid.
    Returns True if start and stop codons are intact, otherwise False.

    Args:
        seq (str): The sequence to check, typically a gene alignment.

    Returns:
        bool: True if start and stop codons are intact, False otherwise.
    """
    # the function is broken when the number reach the index of the first nucl. numbers variable keeps this index
    for numbers in range(0,len(seq)):
        if seq[numbers] != "-":
            break
    start = seq[numbers:numbers+3].upper()
    # the function is broken when the number reach the index of the last nucl. numbers variable keeps this index
    for numbers in range(len(seq)-1,0,-1):
        if seq[numbers] != "-":
            break
    stop = seq[numbers-2:numbers+1].upper()
    if start == "ATG" and stop == "TAG" or start == "ATG" and stop == "TAA" or start == "ATG" and stop == "TGA":
        return True
    else:
        return False


def look_ATG_frameshift(denovo : str, nchit : str) -> (str, int):
    """
    This function checks if an ATG start codon is present within the first 20 nucleotides of a given aligned homolog sequence.

    The function scans the first 20 nucleotides (or the length of the sequence, whichever is shorter) of the alignment to see if 
    any 3 consecutive nucleotides form the start codon 'ATG'. If an 'ATG' is found, the function updates the status to 'S' 
    (indicating the presence of an ATG) and returns the position of the first occurrence of the 'ATG' codon. 

    If no 'ATG' is found within the first 20 nucleotides, the function returns the default status 'A' (absence of ATG) and 
    the position is set to -20, indicating that no valid start codon was detected.

    Parameters:
    -----------
    denovo : str
        An identifier or name for the denovo gene or sequence being analyzed (not used in the logic but may be for tracking).
    
    nchit : str
        The aligned nucleotide sequence (as a string), where '-' represents gaps in the alignment and other characters represent 
        the actual nucleotides.

    Returns:
    --------
    tuple
        A tuple with two elements:
        - presence_atg (str): 'S' if an ATG start codon is found within the first 20 nucleotides, 'A' if not.
        - pos_ATG_in_ali (int): The position of the first 'ATG' codon in the sequence (0-based index), or -20 if no ATG is found.
    """
    presence_atg = "A"
    compteur_nucl = -1
    pos_ATG_in_ali = -20
    for compteur_pos in range(0,len(nchit)):
        if nchit[compteur_pos] != "-":
            compteur_nucl += 1
            nucl = ""
            nucl += nchit[compteur_pos].upper()
            # Look for Start
            if nucl == "A":
                for compteur_pos_next in range(compteur_pos + 1,len(nchit)):
                    if nchit[compteur_pos_next] != "-":
                        nucl += nchit[compteur_pos_next].upper()
                    if len(nucl) == 3:
                        if nucl == "ATG":
                            pos_ATG_in_ali = compteur_pos
                            presence_atg = "S"
                        break
            if presence_atg == "S" or compteur_nucl > 20:
                break
    return presence_atg, pos_ATG_in_ali


def searchATG(alignment : list) -> (str, int):
    """
    Search for the presence of the ATG start codon in an alignment of sequences.

    This function examines an alignment of two sequences: a "denovo" gene sequence and a corresponding homologous sequence ("nchit").
    It looks for the presence of the ATG start codon within the homologous sequence (nchit), accounting for gaps in the alignment.
    If the first non-gap nucleotides in the homologous sequence form an ATG, the function returns "P" (for presence).
    If the start codon is shifted due to a frameshift, the function will call the `look_ATG_frameshift` function to investigate further.
    If no ATG is found, the function returns "A" (for absence).

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two elements:
        - denovo (str): A string representing the denovo gene sequence (used to compare the start codon position).
        - nchit (str): A string representing the homologous sequence (may contain gaps represented by '-').

    Returns:
    --------
    tuple
        A tuple containing:
        - presence_ATG (str): "P" if ATG is found at the beginning of the homologous sequence, "A" if not, or "S" if frameshifted ATG is found.
        - start_atg_in_ali (int): The position of the start ATG codon in the homologous sequence, or -20 if no valid ATG is found.
    """
    presence_ATG = "A"
    start_atg_in_ali = -20

    denovo = alignment[0]
    nchit = alignment[1]
    for numbers_1 in range(0,len(denovo)):
        if denovo[numbers_1] != "-":
            break
    start_nchit = ""
    for numbers_2 in range(numbers_1,len(nchit)):
        if nchit[numbers_2] != "-":
            start_nchit += nchit[numbers_2].upper()
            if len(start_nchit) == 3:
                break
    if start_nchit != "ATG":
        presence_ATG, start_atg_in_ali = look_ATG_frameshift(denovo, nchit)
    else:
        presence_ATG = "P"
    
    return presence_ATG, start_atg_in_ali


def searchStop(alignment : list) -> str:
    """
    Search for the presence of a stop codon in the homologous sequence (nchit).

    This function examines an alignment of two sequences: a "denovo" gene sequence and its homologous sequence ("nchit").
    It searches for the presence of a stop codon (TAG, TAA, or TGA) in the homologous sequence. The stop codon is searched 
    starting from the end of the denovo sequence, accounting for gaps in the alignment. If a stop codon is found in 
    the homologous sequence, the function returns "P" (for presence). If no stop codon is found, it returns "A" (for absence).

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two elements:
        - denovo (str): A string representing the denovo gene sequence (used to compare the stop codon position).
        - nchit (str): A string representing the homologous sequence (may contain gaps represented by '-').

    Returns:
    --------
    presence_stop (str):
        - "P" if a stop codon (TAG, TAA, or TGA) is found in the homologous sequence.
        - "A" if no stop codon is found in the homologous sequence.
    """
    presence_stop = "A"
    list_stops = ["TAG", "TAA", "TGA"]
    denovo = alignment[0]
    nchit = alignment[1]
    for numbers_1 in range(len(denovo)-1,0,-1):
        if denovo[numbers_1] != "-":
            break
    stop_nchit = ""
    for numbers_2 in range(numbers_1,0,-1):
        if nchit[numbers_2] != "-":
            stop_nchit = nchit[numbers_2].upper() + stop_nchit
            if len(stop_nchit) == 3:
                break
    if stop_nchit not in list_stops:
        presence_stop = "A"
    else:
        presence_stop = "P"
    return presence_stop


def calculate_with_schmidt_similarity_score(ali_denovo : str, ali_homolog : str) -> float:
    """
    Calculate the Aaron Schmidt similarity score to assess the percentage of frameshift alignment 
    between two sequences.

    This function calculates the similarity score based on the alignment of two sequences (denovo 
    and homolog) by comparing their nucleotide frames. The score measures how well the frames (1, 2, 
    or 3) of both sequences align, while accounting for gaps ('-') and stop codons (TAA, TAG, TGA).
    
    The algorithm assigns a frame (1, 2, or 3) to each nucleotide in both sequences based on their 
    position in the alignment, then compares these frames between the two sequences to calculate the 
    percentage of positions where the frames match.

    Parameters:
    -----------
    ali_denovo : str
        The aligned denovo sequence, where nucleotides are represented as 'A', 'T', 'G', 'C', or '-'
        for gaps.
    
    ali_homolog : str
        The aligned homologous sequence, which follows the same nucleotide format as ali_denovo.

    Returns:
    --------
    score_schmidt : float
        The Aaron Schmidt similarity score, representing the percentage of matching frames between the 
        denovo and homologous sequences (excluding gaps and stop codons).
    """
    dico_denovo = {}
    dico_homolog = {}
    frame_denovo = 1
    frame_homolog = 1
    # here we attribute a frame 1,2,or 3 for each nucl in both aligned sequences.
    for i in range(0,len(ali_denovo)):
        pos = i + 1
        if ali_denovo[i] != "-":
            dico_denovo[pos] = frame_denovo
            frame_denovo += 1
        else:
            dico_denovo[pos] = "n"
        if ali_homolog[i] != "-":
            dico_homolog[pos] = frame_homolog
            frame_homolog += 1
        else:
            dico_homolog[pos] = "n"
        if frame_denovo > 3:
            frame_denovo = 1
        if frame_homolog > 3:
            frame_homolog = 1

    # calculate the score based on the two dictionaries.
    nb_correct_frame = 0
    nb_position = 0
    triplet = "" 
    stop_scoring = False 
    for i in range(1,len(ali_denovo) + 1):
        # implemet a triplet to double check stops.
        if ali_homolog[i-1] != "-": 
            triplet += ali_homolog[i-1] 
        if triplet == "TAA" or triplet == "TAG" or triplet == "TGA":
            stop_scoring = True
            if dico_denovo[i] == dico_homolog[i]: # add + 1 for the last correct nucl
                nb_correct_frame += 1
        # if for a given position a frame is similar for the two seqs then the score is implemented.
        if dico_denovo[i] == dico_homolog[i] and stop_scoring != True: 
            nb_correct_frame += 1
        if dico_denovo[i] != "n":
            nb_position += 1
        if len(triplet) == 3:
            triplet = "" 
    score_schmidt = nb_correct_frame * 100 / nb_position
    return score_schmidt
    

def calculate_with_schmidt_similarity_score_with_shifted_start(ali_denovo : str, ali_homolog : str, start_atg_in_target : int) -> float:
    """
    Calculate the Aaron Schmidt similarity score for frameshift percentage in case an ATG 
    start codon is shifted within the first 20 nucleotides of the target sequence.

    This function is an adaptation of the Aaron Schmidt similarity score, which compares 
    the frames (1, 2, or 3) of two sequences (denovo and homolog). It handles the case where 
    the target sequence has a shifted ATG within the first 20 nucleotides, by resetting the 
    frame of the homolog sequence starting at the shifted ATG position. The score measures 
    how well the frames of the two sequences align, while excluding gaps ('-') and accounting 
    for stop codons (TAA, TAG, TGA).

    Parameters:
    -----------
    ali_denovo : str
        The aligned denovo sequence, where nucleotides are represented as 'A', 'T', 'G', 'C', 
        or '-' for gaps.
    
    ali_homolog : str
        The aligned homologous sequence, which follows the same nucleotide format as ali_denovo.
    
    start_atg_in_target : int
        The position of the shifted ATG in the target sequence (homolog). This is used to reset 
        the reading frame of the homolog sequence at this position.

    Returns:
    --------
    score_schmidt : float
        The Aaron Schmidt similarity score, representing the percentage of matching frames 
        between the denovo and homolog sequences, excluding gaps and stop codons, and adjusted 
        for a shifted ATG start codon in the target sequence.
    """
    dico_denovo = {}
    dico_homolog = {}
    frame_denovo = 1
    frame_homolog = 1
    # here we attribute a frame 1,2,or 3 for each nucl in both aligned sequences.
    for i in range(0,len(ali_denovo)):
        pos = i + 1
        if ali_denovo[i] != "-":
            dico_denovo[pos] = frame_denovo
            frame_denovo += 1
        else:
            dico_denovo[pos] = "n"
        if i == start_atg_in_target:
            frame_homolog = 1
        if ali_homolog[i] != "-":
            dico_homolog[pos] = frame_homolog
            frame_homolog += 1
        else:
            dico_homolog[pos] = "n"
        if frame_denovo > 3:
            frame_denovo = 1
        if frame_homolog > 3:
            frame_homolog = 1
    # calculate the score based on the two dictionaries.
    nb_correct_frame = 0
    nb_position = 0
    triplet = "" 
    stop_scoring = False 
    for i in range(1,len(ali_denovo) + 1):
        # implemet a triplet to double check stops. In this specific function start building start only where ATG starts.
        if ali_homolog[i-1] != "-" and i > start_atg_in_target : 
            triplet += ali_homolog[i-1] 
        if triplet == "TAA" or triplet == "TAG" or triplet == "TGA":
            stop_scoring = True
            if dico_denovo[i] == dico_homolog[i]: # add + 1 for the last correct nucl
                nb_correct_frame += 1
        # if for a given position a frame is similar for the two seqs then the score is implemented.
        if dico_denovo[i] == dico_homolog[i] and stop_scoring != True: 
            nb_correct_frame += 1
        if dico_denovo[i] != "n":
            nb_position += 1
        if len(triplet) == 3:
            triplet = "" 
    score_schmidt = nb_correct_frame * 100 / nb_position
    return score_schmidt


def searchIndels(alignment : list, start_atg_in_target : int) -> (int, str):
    """
    Search for indels (insertions and deletions) in the alignment and calculate the frameshift score.

    This function identifies the number of indels (gaps) in both the denovo and homolog sequences 
    after removing the unaligned regions (leading and trailing gaps). It then calculates the frameshift 
    percentage score using either the standard Aaron Schmidt similarity score or the adjusted score 
    if a shifted ATG is detected in the target (homolog) sequence. The final output includes the total 
    number of indels and the calculated frameshift score.

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two aligned sequences: the denovo sequence (str) and the homolog sequence (str).
    
    start_atg_in_target : int
        The position of the shifted ATG in the homolog sequence (target). If no shifted ATG is detected, 
        this is set to -20. This value is used to determine how to calculate the frameshift score.

    Returns:
    --------
    total_indels : int
        The total number of indels (gaps) found in both the denovo and homolog sequences after alignment trimming.
    
    score_perc_frameshift : str
        The frameshift percentage score, calculated using the Aaron Schmidt similarity score, adjusted for 
        any detected shifted ATG in the homolog sequence. The score is returned as a string.
    """
    denovo = alignment[0]
    nchit = alignment[1]
    nb_gap_denovo = 0
    for numbers_start in range(0,len(denovo)):
        if denovo[numbers_start] != "-":
            break
    for numbers_stop in range(len(denovo)-1,0,-1):
        if denovo[numbers_stop] != "-":
            break
    seq_denovo_cut = denovo[numbers_start:numbers_stop + 1]
    seq_nchit_cut = nchit[numbers_start:numbers_stop + 1]

    if start_atg_in_target == -20:
        score_perc_frameshift = str(calculate_with_schmidt_similarity_score(seq_denovo_cut, seq_nchit_cut))
    else:
        score_perc_frameshift = str(calculate_with_schmidt_similarity_score_with_shifted_start(denovo, nchit, start_atg_in_target))
    nb_gap_denovo = seq_denovo_cut.count('-')
    nb_gap_nchit = seq_nchit_cut.count('-')
    total_indels = nb_gap_denovo + nb_gap_nchit
    return total_indels, score_perc_frameshift


def searchSubs(alignment : list) -> int:
    """
    Counts the number of substitutions between the denovo sequence and its homolog (nchit) in the given alignment.

    This function compares the aligned regions of the denovo and homolog sequences, and counts the number 
    of positions where the nucleotides are different (substitutions). It ignores gaps in both sequences and 
    only considers positions where nucleotides are present in both sequences.

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two aligned sequences: the denovo sequence (str) and the homolog sequence (str).
    
    Returns:
    --------
    nb_subs : int
        The total number of nucleotide substitutions between the denovo and homolog sequences within the aligned regions.
    """
    nb_subs = 0
    denovo = alignment[0]
    nchit = alignment[1]
    nb_gap_denovo = 0
    for numbers_start in range(0,len(denovo)):
        if denovo[numbers_start] != "-":
            break
    for numbers_stop in range(len(denovo)-1,0,-1):
        if denovo[numbers_stop] != "-":
            break
    for i in range(numbers_start, numbers_stop+1):
        nucl_denovo = denovo[i].upper()
        nucl_nchit = nchit[i].upper()
        if nucl_denovo != "-" and nucl_nchit != "-" and nucl_denovo != nucl_nchit:
            nb_subs += 1
    return nb_subs


def searchPreStop(alignment : list, perc_seq_accepted : float, start_atg_in_target : int) -> (str, int):
    """
    Searches for a premature stop codon in the homolog (nchit) sequence within an aligned region, based on a given threshold.

    This function checks the aligned portion of the homolog sequence (`nchit`) for premature stop codons within a specified percentage
    of the sequence length. The function can also account for the position of a shifted ATG start codon in the homolog sequence, if any.

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two aligned sequences: the denovo sequence (str) and the homolog sequence (str).
    
    perc_seq_accepted : float
        The percentage of the sequence to be considered for the search of premature stop codons. This parameter determines how much of 
        the homolog sequence will be scanned for stop codons.

    start_atg_in_target : int
        The position of the ATG start codon in the homolog sequence, if shifted. If no shift is present, this value should be -20.
    
    Returns:
    --------
    prem_stop_val : str
        Returns "P" if a premature stop codon is found within the specified region, otherwise returns "A" (indicating no premature stop).
    
    pos_premature_stop : int
        The position of the premature stop codon in the homolog sequence if found. If no premature stop is found, returns 0.
    """
    list_stop = ["TAG", "TAA", "TGA"]
    denovo = alignment[0]
    nchit = alignment[1]
    limit_search = int(len(nchit) * perc_seq_accepted / 100)
    premature_stop = ""
    pos_premature_stop = 0
    prem_stop_val = "A"
    for numbers_start in range(0,len(denovo)):
        if denovo[numbers_start] != "-":
            break
    triplet_nchit = ""
    if start_atg_in_target == -20:
        for i in range(numbers_start, limit_search):
            if nchit[i] != "-":
                triplet_nchit += nchit[i]
            if len(triplet_nchit) == 3:
                triplet_nchit.upper()
                if triplet_nchit in list_stop :
                    premature_stop = triplet_nchit
                    pos_premature_stop = i
                    break
                triplet_nchit = ""
        if premature_stop == "":
            prem_stop_val = "A"
        else:
            prem_stop_val = "P"
    else:
        if start_atg_in_target < (limit_search + 3):
            for i in range(start_atg_in_target, limit_search):
                if nchit[i] != "-":
                    triplet_nchit += nchit[i]
                if len(triplet_nchit) == 3:
                    triplet_nchit.upper()
                    if triplet_nchit in list_stop :
                        premature_stop = triplet_nchit
                        pos_premature_stop = i
                        break
                    triplet_nchit = ""
            if premature_stop == "":
                prem_stop_val = "A"
            else:
                prem_stop_val = "P"
    return prem_stop_val, pos_premature_stop

        
def measure_size_intron(denovo_seq : str) -> int:
    """
    Measures the cumulated length of introns in the denovo sequence.

    This function calculates the total length of introns in the given denovo sequence by counting the number of 
    lowercase characters, which represent the intronic regions.

    Parameters:
    -----------
    denovo_seq : str
        A string representing the denovo sequence where lowercase letters indicate introns.

    Returns:
    --------
    size_intron : int
        The cumulative length of the introns, based on the number of lowercase characters in the sequence.
    """
    size_intron = 0
    for i in denovo_seq:
        if i.islower() == True:
            size_intron += 1
    return size_intron


def splice_alignment(alignment : list) -> (list, list):
    """
    Splices alignments that contain introns, extracting the intronic regions and removing them from the alignment.

    This function processes an alignment pair (denovo and nchit) by identifying introns (represented by lowercase letters in the denovo sequence).
    The function removes the intronic regions from both sequences in the alignment and returns the spliced sequences along with a list of the removed introns.

    Parameters:
    -----------
    alignment : list
        A list containing two strings:
        - `denovo`: A denovo sequence where lowercase letters represent introns.
        - `nchit`: A homologous sequence aligned to the denovo sequence.

    Returns:
    --------
    alignment : list
        A list of two strings:
        - The spliced `denovo` sequence with introns removed.
        - The corresponding `nchit` sequence with the introns removed.
    list_introns : list
        A list of strings representing the removed intronic regions from the denovo sequence.
    """
    last_previous = "A"
    denovo = alignment[0]
    nchit = alignment[1]
    new_seq_denovo = ""
    new_seq_hit = ""
    list_introns = []
    seq_intron = ""
    for number in range(0,len(denovo)):
        nucl_denovo = denovo[number]
        if nucl_denovo == "-":
            if last_previous.isupper():
                new_seq_denovo += nucl_denovo
                new_seq_hit += nchit[number]
            else:
                seq_intron += nucl_denovo
        else:
            if nucl_denovo.islower() == False:
                if seq_intron != "":
                    list_introns.append(seq_intron)
                    seq_intron = ""
                last_previous = nucl_denovo
                new_seq_denovo += nucl_denovo
                new_seq_hit += nchit[number]
            else:
                last_previous = nucl_denovo
                seq_intron += nucl_denovo
    if seq_intron != "":
        list_introns.append(seq_intron)
    alignment = [new_seq_denovo,new_seq_hit]
    return alignment,list_introns


def search_transcript_overlap(dico_target_transcript_coordinate : dict, start_nc : int, stop_nc : int, direction_nc : str) -> str:
    """
    Determines the overlap between a transcript's position in the chromosome and the position of the nchit sequence.

    This function checks if a given region in the nchit sequence (defined by `start_nc` and `stop_nc`) overlaps with a transcript's 
    position in the chromosome (given by `dico_target_transcript_coordinate`). It also takes into account the direction of the transcript 
    (`direction_nc`), which can either be 'forward' or 'reverse'.

    Parameters:
    -----------
    dico_target_transcript_coordinate : dict
        A dictionary containing the coordinates of a target transcript on the chromosome. The dictionary has the following structure:
        - [0] : start position of the transcript (integer)
        - [1] : end position of the transcript (integer)
        - [2] : direction of the transcript ('forward' or 'reverse')

    start_nc : int
        The start position of the region of interest in the nchit sequence (chromosome position).

    stop_nc : int
        The end position of the region of interest in the nchit sequence (chromosome position).

    direction_nc : str
        The direction of the nchit sequence, either 'forward' or 'reverse'.

    Returns:
    --------
    value : str
        A string indicating the relationship between the transcript and the region in the nchit sequence:
        - 'forward' if the region overlaps the transcript and the directions match.
        - 'reverse' if the region overlaps the transcript but the directions are opposite.
        - 'A' if there is no overlap between the region and the transcript.
    """
    value = "A"
    if int(dico_target_transcript_coordinate[0]) <= start_nc and int(dico_target_transcript_coordinate[1]) >= start_nc:
        if dico_target_transcript_coordinate[2] == direction_nc:
            value = "forward"
        else:
            value = "reverse"
    if int(dico_target_transcript_coordinate[1]) >= stop_nc and int(dico_target_transcript_coordinate[0]) <= stop_nc:
        if dico_target_transcript_coordinate[2] == direction_nc:
            value = "forward"
        else:
            value = "reverse"
    if int(dico_target_transcript_coordinate[0]) > start_nc and int(dico_target_transcript_coordinate[1]) < stop_nc:
        if dico_target_transcript_coordinate[2] == direction_nc:
            value = "forward"
        else:
            value = "reverse"
    return value


def validateTranscription(name_denovo : str, name_nchom : str, dico_target_transcript_coordinate : dict, dico_target_transcript_blast_hits : dict) -> str:
    """
    This function evaluates if a given nchit (novel transcript) is transcribed based on its overlap with known transcript coordinates
    and the presence of a corresponding blast hit. The function considers whether the transcription is complete or incomplete or other
    based on both the overlap and the BLAST results.

    Parameters:
    -----------
    name_denovo : str
        The identifier of the novel transcript (denovo) to be evaluated.

    name_nchom : str
        The chromosomal position and direction of the homologous sequence, formatted as "start-stop-chrom-direction".

    dico_target_transcript_coordinate : dict
        A dictionary containing the coordinates and directions of known target transcripts. 
        The dictionary should be structured with chromosome names as keys, where the value is another dictionary with 
        transcript IDs as keys and their respective coordinates and direction as values. The coordinates are represented as a tuple:
        (start_position, stop_position, direction).

    dico_target_transcript_blast_hits : dict
        A dictionary containing BLAST hits for the target transcripts. The keys are the denovo transcript names, 
        and the values are dictionaries with transcript IDs as keys and their respective BLAST status (e.g., "complete", "incomplete", "spliced") as values.

    Returns:
    --------
    valueTranscription : str
        A string indicating the transcription status:
        - 'forward' if the region overlaps with the transcript and the directions match.
        - 'reverse' if the region overlaps with the transcript but the directions are opposite.
        - 'spliced' if the transcription is spliced based on the BLAST hits.
        - 'A' if no overlap is detected.
        - 'NA' if the transcription could not be assessed.
    """
    valueTranscription = "NA"
    if len(dico_target_transcript_coordinate) > 0:
        valueTranscription = "A"
        # retrieve the informations about the nchit
        list_elts_nchom = name_nchom.split("-")
        start = int(list_elts_nchom[0])
        stop = int(list_elts_nchom[1])
        chrom = list_elts_nchom[2]
        direction = list_elts_nchom[3]
        # retrieve coordinates of the target transcripts
        if chrom in dico_target_transcript_coordinate: # just added
            sub_dico_transcripts = dico_target_transcript_coordinate[chrom]
            for transcript in sub_dico_transcripts:
                # assess if the transcripts overlap with the coordinates of the nchit
                transc_status = search_transcript_overlap(sub_dico_transcripts[transcript],start,stop,direction)
                if transc_status == "reverse":
                    valueTranscription = transc_status
                    break
                elif transc_status == "forward":
                    # here to see if the transcription is complete or not, taking into account the splicing
                    if name_denovo in dico_target_transcript_blast_hits:
                        dico_denovo_hit = dico_target_transcript_blast_hits[name_denovo]
                        if transcript in dico_denovo_hit:
                            valueTranscription = dico_denovo_hit[transcript]
                            break
                        else:
                            valueTranscription = "spliced"
                    # and of not?? should i then give a transcription status null?
    return valueTranscription


def get_unspliced_seq(alignment : list) -> str:

    """
    This function generates an unspliced sequence by modifying the homologous hit (nchit) alignment. 
    Specifically, it lowers the nucleotides of the homologous sequence (nchit) that align to introns in the denovo sequence. 
    Introns are represented by gaps (i.e., '-') in the denovo sequence, and the corresponding nucleotides in the homologous sequence 
    (nchit) are converted to lowercase.

    Parameters:
    -----------
    alignment : tuple
        A tuple containing two elements:
        - denovo (str): The novel sequence with potential introns represented by '-'.
        - nchit (str): The homologous sequence aligned with the denovo sequence.

    Returns:
    --------
    unspliced_nchit : str
        The unspliced homologous sequence where the nucleotides aligned to introns in the denovo sequence 
        (represented by '-') are converted to lowercase, and other nucleotides are maintained in their original case.
    """
    denovo = alignment[0]
    nchit = alignment[1]
    unspliced_nchit = ""
    for number in range(0,len(denovo)):
        nucl_denovo = denovo[number]
        nucl_hit = nchit[number]
        if nucl_denovo == "-":
            # the next line makes sure we lower the nucleotides that are alined to an intron with a "-" in the denovo
            if len(unspliced_nchit) > 0 and unspliced_nchit[len(unspliced_nchit) - 1].islower() == True:
                new_nucl_hit = nucl_hit.lower()
                unspliced_nchit += new_nucl_hit
            else:
                unspliced_nchit += nucl_hit
        else:
            if nucl_denovo.islower() == False:
                new_nucl_hit = nucl_hit.upper()
            elif nucl_denovo.islower() == True:
                new_nucl_hit = nucl_hit.lower()
            else:
                new_nucl_hit = nucl_hit
            unspliced_nchit += new_nucl_hit
    return unspliced_nchit


def place_lower_letter(alignment : list, denovo_seq_with_intron : str) -> list:
    """
    This function restores the lowercase nucleotides of introns into the denovo sequence, which were removed during the alignment process.
    The lower case letters in the denovo sequence represent the introns, and this function replaces the corresponding positions 
    in the denovo sequence with the lowercase nucleotides from a provided sequence that includes the introns.

    Parameters:
    -----------
    alignment : list
        A list containing two elements:
        - denovo (str): The denovo sequence with gaps ('-') representing introns, without the lowercase nucleotides.
        - nchit (str): The homologous sequence aligned to the denovo sequence (this is not altered).
    
    denovo_seq_with_intron : str
        A string representing the denovo sequence with introns included, where the nucleotides corresponding to introns are in lowercase.
    
    Returns:
    --------
    alignment : list
        A list with the denovo sequence updated to include the lowercase nucleotides in the correct positions. The format of the alignment is unchanged, 
        with the first element being the updated denovo sequence and the second element being the original homologous sequence.
    """
    denovo = alignment[0]
    nchit = alignment[1]
    new_denovo = ""
    position_in_seq = 0
    for nucl in denovo:
        if nucl == "-":
            new_denovo += nucl
        else:
            new_nucl = denovo_seq_with_intron[position_in_seq]
            position_in_seq += 1
            new_denovo += new_nucl
    alignment = [new_denovo,alignment[1]]
    return alignment


def is_ali_empty(alignment : list) -> bool:
    """
    This function checks if the given alignment is empty by attempting to access the first element of the alignment.
    It returns True if the alignment is empty, and False if it contains at least one sequence.

    Parameters:
    -----------
    alignment : list
        A list containing sequences (e.g., denovo and homologous sequences).
        This function checks if the first sequence (alignment[0]) exists.

    Returns:
    --------
    bool
        Returns True if the alignment is empty (i.e., no sequence at index 0), otherwise returns False.
    """
    try:
        _ = alignment[0]
        return False
    except:
        return True


def main_alignment_function(pop_species_name : str, dico_name_size_denovo : dict, dico_denovo_best_hit : dict, dic_denovo : dict, dico_NcHit : dict, dico_target_transcript_coordinate : dict, dico_target_transcript_blast_hits : dict, value_prem_stop : float, dico_format_all_data : dict) -> dict:
    """
    This function performs multiple alignment tasks for denovo sequences, handling potential introns, assessing alignment properties, 
    and compiling the results into a structured dictionary.

    It processes each denovo sequence, aligning it with its homolog (if present), identifying properties such as ATG, stop codon, 
    indels, frameshift, substitutions, premature stop codons, and transcription status. Additionally, it handles cases involving 
    introns by splicing and unsplicing sequences, as well as accounting for various alignment quality criteria.

    Parameters:
    -----------
    pop_species_name : str
        The name of the species being processed.
        
    dico_name_size_denovo : dict
        A dictionary containing the denovo sequence names and their corresponding sequence sizes.
        
    dico_denovo_best_hit : dict
        A dictionary containing the best hit for each denovo sequence.

    dic_denovo : dict
        A dictionary containing the denovo sequences, where the key is the denovo sequence name and the value is the sequence itself.

    dico_NcHit : dict
        A dictionary containing homologous sequences, where the key is the denovo sequence name and the value is the homologous sequence.

    dico_target_transcript_coordinate : dict
        A dictionary containing the target transcript coordinates for each gene or transcript.

    dico_target_transcript_blast_hits : dict
        A dictionary containing blast hits for the target transcript to check for completeness of transcription.

    value_prem_stop : float
        A threshold value to determine premature stop codons in the alignment.

    dico_format_all_data : dict
        A dictionary that stores the results of each alignment and its corresponding properties for each denovo sequence, updated during processing.

    Returns:
    --------
    dict
        The updated `dico_format_all_data` containing the results of the alignment and the associated properties (e.g., ATG position, stop codon, substitutions, etc.) for each denovo sequence.

    Notes:
    ------
    - This function handles sequences with and without introns and performs alignment using a pairwise aligner.
    - If an intron is detected in the sequence, special handling for splicing and unsplicing is performed.
    - The function performs a number of checks to assess alignment quality, such as checking for the presence of ATG, stop codons, and premature stop codons.
    - The function also computes substitution rates, indels, and frameshift information, as well as validating the transcription status based on the alignment.
    """
    score = 0
    aligner = Align.PairwiseAligner() #***
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    for denovo_name in dico_name_size_denovo:
        # if de novo has an homolog
        if denovo_name in dico_NcHit and denovo_name in dic_denovo:
            validated_hit = "P"
            genomic_position_homolog = dico_denovo_best_hit[denovo_name]
            # Here we extract the denovo sequence, with the intron loered in case of intron
            denovo_seq = dic_denovo[denovo_name]
            nchit_seq = dico_NcHit[denovo_name].upper()
            correct_ali = False
            intron = False
            alignment = False
            # Here we assess if there is intron in the seq
            lowercase_letters = [c for c in denovo_seq if c.islower()]
            if len(lowercase_letters) == 0:
                #alignment = pairwise2.align.globalms(denovo_seq, nchit_seq, 2, -1, -1.5, -.1)
                alignment = aligner.align(denovo_seq, nchit_seq) #***
            else:
                size_intron = measure_size_intron(denovo_seq)
                # Only run with sequences with introns smaller than 3000 nucl. Arbitrary.
                if size_intron < 3000:
                    upper_char_denovo_seq = denovo_seq.upper()
                    #alignment = pairwise2.align.globalms(upper_char_denovo_seq, nchit_seq, 2, -1, -1.5, -.1)

                    alignment = aligner.align(upper_char_denovo_seq, nchit_seq) #***
                    intron = True
            
            if alignment != False and is_ali_empty(alignment):
                alignment = False
            
            # if a correct alignment is now stored in the variable "alignment".
            if alignment != False:
                unspliced_target = "NA"
                if len(denovo_seq) < 500:
                    # we sort to get the best alignment only if the sequences are small otherwise it is too long.
                    alignment = sorted(alignment)
                ali_count = 0
                for ali in alignment: #*** added sorted
                    ali_count += 1
                    # an ali is correct if the ATG and stop is not broken by gap in the aligned de novo. Arbitrary.
                    correct_ali = check_ali(ali[0])
                    if correct_ali == True:
                        alignment = ali
                        break
                    if ali_count > 5000:
                        break
                if correct_ali == False:
                    alignment = alignment[0]
                    score += 1
                if intron == True:
                    # In the alignment, lower letters of the denovo where it is intronic.
                    alignment = place_lower_letter(alignment, denovo_seq)
                    # this function in order to get the homologous hit with lower cases when the seq aligns to the denovo intron
                    unspliced_target = get_unspliced_seq(alignment)
                    # here we remove introns from the alignment
                    alignment,list_introns = splice_alignment(alignment)
                spliced_target = alignment[1]
                denovo_in_ali = alignment[0]
                ATG, start_atg_in_target = searchATG(alignment)
                stop = searchStop(alignment)
                indels, frameshift = searchIndels(alignment, start_atg_in_target)
                nb_subs = searchSubs(alignment)
                premature_stop, pos_premature_stop = searchPreStop(alignment,value_prem_stop, start_atg_in_target)
                transcription = validateTranscription(denovo_name, dico_denovo_best_hit[denovo_name], dico_target_transcript_coordinate, dico_target_transcript_blast_hits)
            else:
                ATG = "NA"
                stop = "NA"
                indels, frameshift = "NA","NA"
                nb_subs = "NA"
                premature_stop, pos_premature_stop = "NA","NA"
                transcription = "NA"
                spliced_target = "NA"
                unspliced_target = "NA"
                denovo_in_ali = "NA"
        else:
            validated_hit = "A"
            genomic_position_homolog = "NA"
            ATG = "NA"
            stop = "NA"
            indels, frameshift = "NA","NA"
            nb_subs = "NA"
            premature_stop, pos_premature_stop = "NA","NA"
            transcription = "NA"
            spliced_target = "NA"
            unspliced_target = "NA"
            denovo_in_ali = "NA"
        list_properties = [validated_hit, genomic_position_homolog, ATG, stop, indels, frameshift, nb_subs, premature_stop, pos_premature_stop, transcription, unspliced_target, spliced_target, denovo_in_ali]
        if denovo_name not in dico_format_all_data:
            dico_format_all_data[denovo_name] = {pop_species_name:list_properties}
        else:
            dico_format_all_data[denovo_name][pop_species_name] = list_properties
    return dico_format_all_data         


def get_transcription_indication(dico_target_species_lines : dict, target_name : str, dico_dot_special_case : dict, unoriented_transc_exclusion : str) -> dict:
    """
    Extracts the transcript positions, names, and coordinates from the transcriptome of a specified target species.

    This function parses a GTF file containing transcript annotations (if provided) for the target species, and retrieves 
    the coordinates (start and stop) and directional information for each transcript. The resulting information is stored 
    in a dictionary where the keys are chromosome names and the values are dictionaries of transcript names with their 
    associated coordinates and direction.

    Parameters:
    -----------
    dico_target_species_lines : dict
        A dictionary containing the species data, where each species name maps to a set of annotations.
        
    target_name : str
        The name of the target species for which the transcriptome data should be extracted.
        
    dico_dot_special_case : dict
        A dictionary mapping the GTF direction field (either "+" or "-") to the desired string representation ("f" for forward, 
        "r" for reverse). It is used to handle special cases in the GTF file's orientation field.
        
    unoriented_transc_exclusion : str
        A flag indicating whether to exclude transcripts with an undefined orientation (".") from the extraction.
        If set to "True", such transcripts will be ignored; otherwise, they will be included.

    Returns:
    --------
    dict
        A dictionary (`dico_transcript_position`) where the keys are chromosome names and the values are dictionaries containing
        transcript names as keys, with each transcript's coordinates (start, stop) and orientation (direction) as values.

    Notes:
    ------
    - This function assumes that the transcriptome is provided as a GTF file, which is processed line by line.
    - If a transcript has no orientation (represented by "." in the GTF file), the `unoriented_transc_exclusion` flag controls whether it will be included in the results.
    - The GTF file format is expected to have the following fields: chromosome, feature type (must be "transcript"), start, stop, and direction.
    """
    dico_transcript_position = {}
    # First look at whether there is a provided transcriptome
    if "transcriptome_gtf" in dico_target_species_lines[target_name]:
        my_gtf_file = openFile(dico_target_species_lines[target_name]["transcriptome_gtf"])
        # go through the lines of the transcriptome annotation file
        for line in my_gtf_file:
            if line[0] != "#":
                elts_line = line.split()
                if elts_line[2] == "transcript":
                    chrom = elts_line[0]
                    start = int(elts_line[3])
                    stop = int(elts_line[4])
                    if unoriented_transc_exclusion == "True" and elts_line[6] == ".":
                        continue
                    else:
                        # get the sign direction and correct the "." 
                        signe_direction = dico_dot_special_case[elts_line[6]]
                        if signe_direction == "+":
                            direction = "f"
                        else:
                            direction = "r"
                        transcript_name = elts_line[11].split("\"")[1]
                        if chrom not in dico_transcript_position:
                            dico_transcript_position[chrom] = {transcript_name:[start,stop,direction]}
                        else:
                            dico_transcript_position[chrom][transcript_name] = [start,stop,direction]
    return dico_transcript_position


def create_final_output_file(dico : dict, name_output_directory : str, name_intermediate_directory : str) -> None:
    """
    Creates a final output file summarizing various properties for each gene in the given dictionary.

    This function generates a CSV file where each row corresponds to a gene from the `dico` dictionary.
    For each gene, the relevant properties for each species are written in the format specified. The columns 
    include information about sequence alignment, transcription status, indels, substitutions, premature stops, 
    and more. Additionally, the intermediate files used for processing are deleted.

    Parameters:
    -----------
    dico : dict
        A nested dictionary where genes are keys, and their corresponding values are dictionaries with species names 
        as keys and a list of relevant properties (validated hit, genomic position, start/stop, etc.) as values.
        
    name_output_directory : str
        The directory path where the final output file will be saved.
        
    name_intermediate_directory : str
        The directory path where intermediate files (e.g., .dmnd files) are stored, which will be removed after the process.

    Returns:
    --------
    None
        This function does not return any value, but it writes the processed data to a CSV file at the specified location.

    Side Effects:
    -------------
    - Writes the output data to a file named `Table_output_step3.txt` in the `name_output_directory`.
    - Removes intermediate `.dmnd` files from the `name_intermediate_directory`.

    Notes:
    ------
    - The output file format is CSV, with columns separated by commas.
    - The first row of the CSV file contains headers describing the information in each column.
    - Each subsequent row contains data for a specific gene and species, with the gene name, species name, 
      and its associated properties listed.
    """
    link_final_file = name_output_directory + "/Table_output_step3.txt"
    my_file = open(link_final_file, "w")
    my_file.write("denovo,target_genome,validated_hit,genomic_position_homolog,start,stop,Indels,perc_seq_not_affected_by_frameshift,substitutions,premature_stop,pos_premature_stop,transcription_status,intron_lowered_homolog,homolog_in_ali,denovo_in_ali" + "\n")
    #[validated_hit, ATG, stop, indels, frameshift, nb_subs, premature_stop, pos_premature_stop, transcription, splicingSite, motifs]
    for gene in dico:
        for species in dico[gene]:
            list_elts = dico[gene][species]
            my_file.write(gene + "," + species)
            for elt in list_elts:
                my_file.write("," + str(elt))
            my_file.write("\n")
    my_file.close()
    os.system("rm " + name_intermediate_directory + "/*.dmnd")
