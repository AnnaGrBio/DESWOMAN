import os
from Bio import SeqIO


def create_new_folder_intermediate(
    name_output_directory: str, name_intermediate_directory: str
) -> None:
    """
    Creates new directories for the outgroup and intermediate folders if they do not already exist.

    This function checks if the specified outgroup and intermediate directories exist.
    If any of them do not exist, the function creates them using the `mkdir` command.

    Parameters:
    -----------
    name_output_directory : str
        The path of the directory to be created for the results.

    name_intermediate_directory : str
        The path of the directory to be created for intermediate files.
    """
    if os.path.isdir(name_output_directory) is False:
        os.system("mkdir " + name_output_directory)
    if os.path.isdir(name_intermediate_directory) is False:
        os.system("mkdir " + name_intermediate_directory)


def create_unspliced_orfs(
    name_output_directory: str, dict_unspliced_orfs_fasta_lowered_introns: dict
) -> None:
    """
    Creates a FASTA file containing unspliced ORFs with lowered introns.

    This function generates a FASTA file where each entry corresponds to an unspliced ORF
    with the introns removed. The ORFs are written to a file in the specified directory.

    Parameters:
    -----------
    name_output_directory : str
        The path to the directory where the FASTA file will be created.

    dict_unspliced_orfs_fasta_lowered_introns : dict
        A dictionary containing ORF names as keys and their corresponding sequences (with lowered introns) as values.
    """
    name_file = name_output_directory + "/denovo_unspliced_lowered_introns.fa"
    my_file = open(name_file, "w")
    for name in dict_unspliced_orfs_fasta_lowered_introns:
        my_file.write(">" + name + "\n")
        my_file.write(dict_unspliced_orfs_fasta_lowered_introns[name] + "\n")
    my_file.close()


def create_denovo_nucl_file(
    name_output_directory: str, dict_all_ORFs_filtered_with_stop: dict
) -> None:
    """
    Creates a FASTA file containing filtered ORFs, including stop codons.

    This function generates a FASTA file where each entry corresponds to an ORF sequence that has been
    filtered and includes stop codons. The sequences are written to a file in the specified directory.

    Parameters:
    -----------
    name_output_directory : str
        The path to the directory where the FASTA file will be created.

    dict_all_ORFs_filtered_with_stop : dict
        A dictionary containing ORF names as keys and their corresponding sequences (with stop codons included) as values.
    """
    name_file = name_output_directory + "/denovo_nucl.fa"
    my_file = open(name_file, "w")
    for name in dict_all_ORFs_filtered_with_stop:
        my_file.write(">" + name + "\n")
        my_file.write(dict_all_ORFs_filtered_with_stop[name].upper() + "\n")
    my_file.close()


def create_denovo_transcripts_file(
    name_output_directory: str, dict_transcript_fasta_denovo: dict
) -> None:
    """
    Creates a FASTA file containing de novo transcript sequences.

    This function generates a FASTA file where each entry corresponds to a de novo transcript sequence.
    The sequences are written to a file in the specified directory.

    Parameters:
    -----------
    name_output_directory : str
        The path to the directory where the FASTA file will be created.

    dict_transcript_fasta_denovo : dict
        A dictionary containing transcript names as keys and their corresponding sequences as values.
    """
    name_file = name_output_directory + "/denovo_transcripts.fa"
    my_file = open(name_file, "w")
    for name in dict_transcript_fasta_denovo:
        my_file.write(">" + name + "\n")
        my_file.write(dict_transcript_fasta_denovo[name].upper() + "\n")
    my_file.close()


def create_denovo_prot_file(name_output_directory: str) -> None:
    """
    Creates a FASTA file containing de novo protein sequences translated from nucleotide sequences.

    This function takes a nucleotide FASTA file (denovo_nucl.fa) and translates each nucleotide sequence
    into its corresponding protein sequence. The translated protein sequences are then written to a new
    FASTA file (denovo_protein.fa) in the specified directory.

    Parameters:
    -----------
    name_output_directory : str
        The path to the directory where both the nucleotide FASTA file (denovo_nucl.fa) and
        the resulting protein FASTA file (denovo_protein.fa) are located and will be created.
    """
    name_nucl_file = name_output_directory + "/denovo_nucl.fa"
    name_prot_file = name_output_directory + "/denovo_protein.fa"
    my_file = open(name_prot_file, "w")
    for record in SeqIO.parse(name_nucl_file, "fasta"):
        protein_id = record.id
        protein = str(record.seq.translate(to_stop=True))
        my_file.write(">" + str(protein_id) + "\n")
        my_file.write(protein + "\n")
    my_file.close()


def create_info_file_orfs(
    name_output_directory: str,
    dict_transcripts_properties: dict,
    dict_gene_status: dict,
    dict_gene_transcripts_denovo: dict,
    dict_transcrit_filtered_orfs: dict,
    dict_pos_unspliced_ORFs_plus_stop: dict,
    sublist_transcripts_overlap: list,
) -> None:
    """
    Creates a CSV file that contains detailed information about genes, transcripts, and associated ORFs, including genomic locations and overlap types.

    This function generates a file named `information_file.txt` in the specified output directory, which contains tabular information for each gene's denovo transcripts and the corresponding ORFs.
    The file includes details such as the gene name, chromosome, orientation, gene status, number of denovo transcripts, transcript details (name, start and stop in genome), ORF details (name, genomic start and stop), and whether the transcript overlaps with other genomic regions (genic, intergenic, intronic, or antisense).

    Parameters:
    -----------
    name_output_directory : str
        Path to the directory where the information file will be created (e.g., `genome/outgroup`).

    dict_transcripts_properties : dict
        A dictionary containing transcript properties, where keys are transcript names and values are lists with the following information:
        - Chromosome
        - Orientation
        - Start position in genome
        - Stop position in genome

    dict_gene_status : dict
        A dictionary with gene names as keys and gene statuses (e.g., "genic" or "non-genic") as values.

    dict_gene_transcripts_denovo : dict
        A dictionary where keys are gene names and values are lists of transcripts associated with each gene.

    dict_transcrit_filtered_orfs : dict
        A dictionary where keys are transcript names and values are lists of filtered ORF names associated with each transcript.

    dict_pos_unspliced_ORFs_plus_stop : dict
        A dictionary where keys are ORF names and values are tuples/lists containing the start and stop positions of the ORFs in the genome.

    sublist_transcripts_overlap : list
        A list containing sublists that classify transcript overlap types. The sublists should categorize:
        - [0]: intergenic transcripts
        - [1]: intronic transcripts
        - [2]: antisense transcripts
        - Any transcript not falling into the above categories is assumed to be "genic_transcripts".
    """
    name_info_file = name_output_directory + "/information_file.txt"
    my_file = open(name_info_file, "w")
    my_file.write(
        "gene_name,chromosome,orientation,putative_genic_variant,nb_denovo_transcripts,transcript_name,transcript_start_in_genome_B1,transcript_end_in_genome_B1,nb_filtered_orf,neORF_candidate,neORF_start_in_genome_B0,neORF_end_in_genome_B0,genomic_overlap"
        + "\n"
    )
    for gene_name in dict_gene_transcripts_denovo:
        list_transcripts = dict_gene_transcripts_denovo[gene_name]
        for transcript_name in list_transcripts:
            transcript_overlap = ""
            if transcript_name in sublist_transcripts_overlap[0]:
                transcript_overlap = "intergenic_transcripts"
            elif transcript_name in sublist_transcripts_overlap[1]:
                transcript_overlap = "intronic_transcripts"
            elif transcript_name in sublist_transcripts_overlap[2]:
                transcript_overlap = "antisense_transcript"
            else:
                transcript_overlap = "genic_transcripts"
            if transcript_name in dict_transcrit_filtered_orfs:
                list_orf = dict_transcrit_filtered_orfs[transcript_name]
                for orf_name in list_orf:
                    if orf_name in dict_pos_unspliced_ORFs_plus_stop:
                        my_gene_name = gene_name
                        my_gene_status = dict_gene_status[my_gene_name]
                        my_transcript_name = transcript_name
                        my_chrom = dict_transcripts_properties[my_transcript_name][0]
                        my_direction = dict_transcripts_properties[my_transcript_name][
                            1
                        ]
                        my_nb_denovo_transcripts = str(len(list_transcripts))
                        my_transcript_start = str(
                            dict_transcripts_properties[my_transcript_name][2]
                        )
                        my_transcript_stop = str(
                            dict_transcripts_properties[my_transcript_name][3]
                        )
                        my_nb_filtered_orfs = str(len(list_orf))
                        my_orf_name = orf_name
                        my_orf_start_in_genome = str(
                            dict_pos_unspliced_ORFs_plus_stop[my_orf_name][0]
                        )
                        my_orf_stop_in_genome = str(
                            dict_pos_unspliced_ORFs_plus_stop[my_orf_name][1]
                        )
                        if int(my_orf_start_in_genome) < int(my_orf_stop_in_genome):
                            my_file.write(
                                my_gene_name
                                + ","
                                + my_chrom
                                + ","
                                + my_direction
                                + ","
                                + my_gene_status
                                + ","
                                + my_nb_denovo_transcripts
                                + ","
                                + my_transcript_name
                                + ","
                                + my_transcript_start
                                + ","
                                + my_transcript_stop
                                + ","
                                + my_nb_filtered_orfs
                                + ","
                                + my_orf_name
                                + ","
                                + my_orf_start_in_genome
                                + ","
                                + my_orf_stop_in_genome
                                + ","
                                + transcript_overlap
                                + "\n"
                            )
                        else:
                            my_file.write(
                                my_gene_name
                                + ","
                                + my_chrom
                                + ","
                                + my_direction
                                + ","
                                + my_gene_status
                                + ","
                                + my_nb_denovo_transcripts
                                + ","
                                + my_transcript_name
                                + ","
                                + my_transcript_start
                                + ","
                                + my_transcript_stop
                                + ","
                                + my_nb_filtered_orfs
                                + ","
                                + my_orf_name
                                + ","
                                + my_orf_stop_in_genome
                                + ","
                                + my_orf_start_in_genome
                                + ","
                                + transcript_overlap
                                + "\n"
                            )
    my_file.close()
