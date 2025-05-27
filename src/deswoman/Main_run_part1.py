from deswoman.module_transcripts_data_and_threeshold import (
    get_path_gtf,
    get_path_to_fasta,
    get_path_gff,
    store_transcriptome_seq,
    associate_gene_and_transcripts,
    get_transcripts_properties,
    get_list_name_query_rep,
)
from deswoman.module_orfs import (
    my_get_orfs,
    build_dict_orfs_per_transcript,
    generate_list_genes_objects,
    sort_orfs_by_properties,
)
from deswoman.module_unspliced_orfs import (
    get_dot_orientation,
    build_dict_unspliced_seqs,
    rebuild_spliced_orfs,
    add_stop_codon_to_spliced_orfs,
)
from deswoman.module_lower_intron import build_unspliced_fasta_seq_lowered_intron
from deswoman.module_output_part1 import (
    create_new_folder_intermediate,
    create_unspliced_orfs,
    create_info_file_orfs,
    create_denovo_nucl_file,
    create_denovo_prot_file,
    create_denovo_transcripts_file,
)
from deswoman.module_my_bedtool import (
    my_bedtool_genes_introns_ordered,
    my_bedtool_extract_overlap,
    extract_denovo_transcripts,
)
from deswoman.module_reformat_dico_with_IDs import (
    implement_dict_transcripts_properties,
    implement_dict_gene_status,
    implement_dict_gene_transcripts_denovo,
    implement_dict_transcrit_filtered_orfs,
    implement_dict_pos_unspliced_ORFs_plus_stop,
    implement_sublist_transcripts_overlap,
    implement_dict_all_ORFs_filtered_with_stop,
    implement_dict_transcript_fasta_denovo,
)
from deswoman.module_colors import *


def run_part_1_strat1(dico_variables: dict) -> bool:
    """
    Executes all steps of Part 1 for Strategy 1 in the DESwoMAN pipeline.

    This function processes transcriptomic and genomic data to extract, filter,
    and analyze de novo transcripts and ORFs (Open Reading Frames). It performs
    the following main steps:

    1. **Data Retrieval**: Loads transcript and genome datasets in various formats.
    2. **Transcript Processing**: Extracts and filters transcripts based on TPM
       threshold and user-defined criteria.
    3. **Genomic Overlap Analysis**: Identifies transcripts overlapping with genomic
       features and filters out unwanted candidates.
    4. **ORF Extraction & Filtering**: Retrieves ORFs from de novo transcripts,
       sorts them based on specific properties, and validates their genomic presence.
    5. **Unspliced ORF Reconstruction**: Aligns ORFs with genome data to ensure
       correctness and adds stop codons if necessary.
    6. **Output Generation**: Creates various output files including FASTA sequences,
       transcript information, and ORF details.

    Args:
        dico_variables (dict): A dictionary containing user-defined parameters
        and paths required for processing.

    Returns:
        bool: True if at least one valid ORF is found and processed successfully,
        False otherwise.

    Notes:
        - If no de novo transcripts or ORFs are found, processing stops early.
        - Outputs are saved in user-defined directories.
        - Various external modules handle specific processing steps.
    """
    validated_step = False
    print(" ")

    ### Step 1 Retrieve all transcripts from data, get orfs, sort orfs by removing the ones that end with the end of transcript and the ones that are reverse
    ## extract links to datasets (genomes and transcriptomes in different format)
    link_to_query_gtf = get_path_gtf(
        dico_variables["path_to_transcriptome_repository"], dico_variables["query"]
    )
    link_to_query_transcriptome = get_path_to_fasta(
        dico_variables["path_to_transcriptome_repository"], dico_variables["query"]
    )
    link_to_query_genome_gff = get_path_gff(
        dico_variables["path_to_genome_repository"], dico_variables["query"]
    )
    link_to_query_genome = get_path_to_fasta(
        dico_variables["path_to_genome_repository"], dico_variables["query"]
    )
    name_output_directory = dico_variables["path_output"]
    name_intermediate_directory = dico_variables["path_output_intermediate"]
    create_new_folder_intermediate(name_output_directory, name_intermediate_directory)

    ## store transcripts
    dict_transcript_fasta_all = store_transcriptome_seq(link_to_query_transcriptome)
    # open gtf file from query transcriptome in readline mode
    gtf_file = openFile(link_to_query_gtf)
    # create dict with genes associated to a list of their spliced variants (transcripts names) only when spliced variants have a TPM value high enough according to user choice
    dict_gene_transcripts_all = associate_gene_and_transcripts(
        gtf_file,
        dico_variables["TPM_threeshold"],
        dict_transcript_fasta_all,
        dico_variables["filter_TE"],
        dico_variables["rm_undir_transc"],
    )
    # create a dictionary with all transcripts propoerties
    dict_transcripts_properties = get_transcripts_properties(gtf_file)
    # get the arbitrary orientation of the "." transcripts in the transcriptome assembly
    get_dot_orientation(
        link_to_query_genome,
        dict_transcripts_properties,
        dict_transcript_fasta_all,
        gtf_file,
        name_intermediate_directory,
    )

    ## Part 2 get transcript overlap to genome elements and withdraw de novo transcripts
    # get dico of all introns and exons sorted.
    dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos = (
        my_bedtool_genes_introns_ordered(link_to_query_genome_gff)
    )
    # creates sublists with intergenic, genic etc transcripts[intergenic_transcripts, intronic_transcripts, antisense_transcript, genic_transcripts]
    sublist_transcripts_overlap = my_bedtool_extract_overlap(
        dict_transcripts_properties,
        dict_gene_transcripts_all,
        dico_sorted_chrom_exon_pos,
        dico_sorted_chrom_intron_pos,
        name_intermediate_directory,
    )
    # creates dictionnaries with denovo transcripts candidates based on the user choice of genomic location (and previous filters)
    (
        dict_gene_transcripts_denovo,
        dict_transcript_fasta_denovo,
        dict_transcrip_status,
        dict_gene_status,
    ) = extract_denovo_transcripts(
        sublist_transcripts_overlap,
        dico_variables["transcript_overlap"],
        dict_transcript_fasta_all,
        dict_gene_transcripts_all,
        name_output_directory,
    )

    # Part3 : retrieve orfs from de novo transcripts and reduce them by properties
    # extracts all ORFs from all transcripts
    dict_all_ORFs = my_get_orfs(name_output_directory)
    if len(dict_all_ORFs) > 0:
        print(
            BRIGHT_BLUE
            + "Total number of correct sense ORFs in denovo transcripts : "
            + RESET
            + str(len(dict_all_ORFs))
        )
        # Creates a dictionary with transcripts as keys and lists of their corresponding ORFs as items
        dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_all_ORFs)
        # creates a list of all genes taken as classes, that contain transcripts and associated ORFs.
        list_gene_object = generate_list_genes_objects(
            dict_gene_transcripts_denovo,
            dict_transcrit_all_orfs,
            dict_all_ORFs,
            dict_transcript_fasta_denovo,
        )
        # goes through each gene object created before and retreive the ORFs we want to keep from all the splice variants.
        dict_transcrit_filtered_orfs, dict_all_ORFs_filtered = sort_orfs_by_properties(
            dico_variables["filter_genic"],
            list_gene_object,
            dico_variables["ORFs_choice"],
            dict_all_ORFs,
            dict_gene_status,
        )
        print(
            BRIGHT_BLUE
            + "Total number of correct sense ORFs in denovo transcripts filtered : "
            + RESET
            + str(len(dict_all_ORFs_filtered))
        )

        # Part 4 : retrieve unspliced seqs of orfs and control steps
        # Creates a dictionary where keys are transcript names, and values are lists containing chromosome, direction of transcription, and a sublist with exon start and stop positions.
        dict_gene_exon_pos = build_dict_unspliced_seqs(gtf_file)
        # reconstructs spliced transcripts based on GTF information, extracts spliced ORFs based on the GTF information
        dict_pos_unspliced_ORFs_plus_stop = rebuild_spliced_orfs(
            link_to_query_genome,
            dict_gene_exon_pos,
            dict_all_ORFs_filtered,
            name_intermediate_directory,
        )
        # add stop to the ORF. This function is due to the fact that get ORF used on the first version removed the stop codon.
        dict_all_ORFs_filtered_with_stop = add_stop_codon_to_spliced_orfs(
            dict_pos_unspliced_ORFs_plus_stop, dict_all_ORFs_filtered
        )
        print(
            BRIGHT_BLUE
            + "Total number of filtered ORFs correctly detected in genomes : "
            + RESET
            + str(len(dict_all_ORFs_filtered_with_stop))
        )

        # Part 5 : build unspliced orfs with intron
        # generates unspliced FASTA sequences with lowered introns based on provided genomic information
        dict_unspliced_orfs_fasta_lowered_introns = (
            build_unspliced_fasta_seq_lowered_intron(
                dict_pos_unspliced_ORFs_plus_stop,
                link_to_query_genome,
                dict_gene_exon_pos,
            )
        )

        # Part 6 : build output
        # generate all output files (fasta + info file) of denovo candidates filtered after step 1.
        create_unspliced_orfs(
            name_output_directory, dict_unspliced_orfs_fasta_lowered_introns
        )
        create_info_file_orfs(
            name_output_directory,
            dict_transcripts_properties,
            dict_gene_status,
            dict_gene_transcripts_denovo,
            dict_transcrit_filtered_orfs,
            dict_pos_unspliced_ORFs_plus_stop,
            sublist_transcripts_overlap,
        )
        create_denovo_nucl_file(name_output_directory, dict_all_ORFs_filtered_with_stop)
        create_denovo_prot_file(name_output_directory)
        create_denovo_transcripts_file(
            name_output_directory, dict_transcript_fasta_denovo
        )
        validated_step = True

    else:
        print(BRIGHT_MAGENTA + "No denovo transcript found ! " + RESET)
    return validated_step


def run_part_1_strat2(dico_variables: dict) -> bool:
    """
    Executes all steps of Part 1 for Strategy 2 in the DESwoMAN pipeline.

    This function processes transcriptomic and genomic data to extract, filter,
    and analyze de novo transcripts and ORFs (Open Reading Frames). It performs
    the following main steps:

    1. **Data Retrieval**: Loads transcript and genome datasets in various formats.
    2. **Transcript Processing**: Extracts and filters transcripts from each transcriptome based on TPM
       threshold and user-defined criteria.
    3. **Genomic Overlap Analysis**: Identifies transcripts overlapping with genomic
       features and filters out unwanted candidates.
    4. **ORF Extraction & Filtering**: Retrieves ORFs from de novo transcripts,
       sorts them based on specific properties, and validates their genomic presence.
    5. **Unspliced ORF Reconstruction**: Aligns ORFs with genome data to ensure
       correctness and adds stop codons if necessary.
    6. **Output Generation**: Creates various output files including FASTA sequences,
       transcript information, and ORF details.

    Args:
        dico_variables (dict): A dictionary containing user-defined parameters
        and paths required for processing.

    Returns:
        bool: True if at least one valid ORF is found and processed successfully,
        False otherwise.

    Notes:
        - If no de novo transcripts or ORFs are found, processing stops early.
        - Outputs are saved in user-defined directories.
        - Various external modules handle specific processing steps.
    """
    validated_step = False
    ## save total final output
    dict_transcripts_properties_total = {}
    dict_gene_status_total = {}
    dict_gene_transcripts_denovo_total = {}
    dict_transcrit_filtered_orfs_total = {}
    dict_pos_unspliced_ORFs_plus_stop_total = {}
    sublist_transcripts_overlap_total = []
    dict_all_ORFs_filtered_with_stop_total = {}
    dict_transcript_fasta_denovo_total = {}
    ## extract links to datasets (genomes and transcriptomes in different format)
    link_to_query_genome_gff = get_path_gff(
        dico_variables["path_to_genome_repository"], dico_variables["query"]
    )
    link_to_query_genome = get_path_to_fasta(
        dico_variables["path_to_genome_repository"], dico_variables["query"]
    )
    list_query_name_transc = get_list_name_query_rep(
        dico_variables["path_to_transcriptome_repository"]
    )
    name_output_directory = dico_variables["path_output"]
    name_intermediate_directory = dico_variables["path_output_intermediate"]
    create_new_folder_intermediate(name_output_directory, name_intermediate_directory)

    # get dico of all introns and exons sorted from the commun genome file
    dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos = (
        my_bedtool_genes_introns_ordered(link_to_query_genome_gff)
    )
    # here goes through all given transcriptomes
    for ID_transcriptome in list_query_name_transc:
        print("")
        print(YELLOW + "-----------" + RESET)
        print(
            BRIGHT_BLUE
            + "Name query : "
            + RESET
            + BRIGHT_CYAN
            + ID_transcriptome
            + RESET
        )
        link_to_query_gtf = get_path_gtf(
            dico_variables["path_to_transcriptome_repository"], ID_transcriptome
        )
        link_to_query_transcriptome = get_path_to_fasta(
            dico_variables["path_to_transcriptome_repository"], ID_transcriptome
        )

        ## store transcripts
        dict_transcript_fasta_all = store_transcriptome_seq(link_to_query_transcriptome)
        gtf_file = openFile(
            link_to_query_gtf
        )  # Note I have also to deal with GFF files
        # create dict with genes associated to a list of their spliced variants (transcripts names) only when spliced variants have a TPM value high enough according to user choice
        dict_gene_transcripts_all = associate_gene_and_transcripts(
            gtf_file,
            dico_variables["TPM_threeshold"],
            dict_transcript_fasta_all,
            dico_variables["filter_TE"],
            dico_variables["rm_undir_transc"],
        )
        # create a dictionary with all transcripts propoerties
        dict_transcripts_properties = get_transcripts_properties(gtf_file)
        # get the arbitrary orientation of the "." transcripts in the transcriptome assembly
        get_dot_orientation(
            link_to_query_genome,
            dict_transcripts_properties,
            dict_transcript_fasta_all,
            gtf_file,
            name_intermediate_directory,
        )

        ## Part 2 get transcript overlap to genome elements and withdraw de novo transcripts
        # creates sublists with intergenic, genic etc transcripts[intergenic_transcripts, intronic_transcripts, antisense_transcript, genic_transcripts]
        sublist_transcripts_overlap = my_bedtool_extract_overlap(
            dict_transcripts_properties,
            dict_gene_transcripts_all,
            dico_sorted_chrom_exon_pos,
            dico_sorted_chrom_intron_pos,
            name_intermediate_directory,
        )
        # creates dictionnaries with denovo transcripts candidates based on the user choice of genomic location (and previous filters)
        (
            dict_gene_transcripts_denovo,
            dict_transcript_fasta_denovo,
            dict_transcrip_status,
            dict_gene_status,
        ) = extract_denovo_transcripts(
            sublist_transcripts_overlap,
            dico_variables["transcript_overlap"],
            dict_transcript_fasta_all,
            dict_gene_transcripts_all,
            name_output_directory,
        )

        # Part3 : retrieve orfs from de novo transcripts and reduce them by properties
        # extracts all ORFs from all transcripts
        dict_all_ORFs = my_get_orfs(name_output_directory)
        if len(dict_all_ORFs) > 0:
            print(
                BRIGHT_BLUE
                + "Total number of correct sense ORFs in denovo transcripts : "
                + RESET
                + str(len(dict_all_ORFs))
            )
            # Creates a dictionary with transcripts as keys and lists of their corresponding ORFs as items
            dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_all_ORFs)
            # creates a list of all genes taken as classes, that contain transcripts and associated ORFs.
            list_gene_object = generate_list_genes_objects(
                dict_gene_transcripts_denovo,
                dict_transcrit_all_orfs,
                dict_all_ORFs,
                dict_transcript_fasta_denovo,
            )
            # goes through each gene object created before and retreive the ORFs we want to keep from all the splice variants.
            dict_transcrit_filtered_orfs, dict_all_ORFs_filtered = (
                sort_orfs_by_properties(
                    dico_variables["filter_genic"],
                    list_gene_object,
                    dico_variables["ORFs_choice"],
                    dict_all_ORFs,
                    dict_gene_status,
                )
            )
            print(
                BRIGHT_BLUE
                + "Total number of correct sense ORFs in denovo transcripts filtered : "
                + RESET
                + str(len(dict_all_ORFs_filtered))
            )

            # Part 4 : retrieve unspliced seqs of orfs and control steps
            # Creates a dictionary where keys are transcript names, and values are lists containing chromosome, direction of transcription, and a sublist with exon start and stop positions.
            dict_gene_exon_pos = build_dict_unspliced_seqs(gtf_file)
            # reconstructs spliced transcripts based on GTF information, extracts spliced ORFs based on the GTF information
            dict_pos_unspliced_ORFs_plus_stop = rebuild_spliced_orfs(
                link_to_query_genome,
                dict_gene_exon_pos,
                dict_all_ORFs_filtered,
                name_intermediate_directory,
            )
            # add stop to the ORF. This function is due to the fact that get ORF used on the first version removed the stop codon.
            dict_all_ORFs_filtered_with_stop = add_stop_codon_to_spliced_orfs(
                dict_pos_unspliced_ORFs_plus_stop, dict_all_ORFs_filtered
            )
            print(
                BRIGHT_BLUE
                + "Total number of filtered ORFs correctly detected in genomes : "
                + RESET
                + str(len(dict_all_ORFs_filtered_with_stop))
            )

            # Part 6 : build output
            # reformat all transcripts, genes and ORF names in all dictionaries so that they contain the transcriptome ID
            implement_dict_transcripts_properties(
                dict_transcripts_properties,
                dict_transcripts_properties_total,
                ID_transcriptome,
            )
            implement_dict_gene_status(
                dict_gene_status, dict_gene_status_total, ID_transcriptome
            )
            implement_dict_gene_transcripts_denovo(
                dict_gene_transcripts_denovo,
                dict_gene_transcripts_denovo_total,
                ID_transcriptome,
            )
            implement_dict_transcrit_filtered_orfs(
                dict_transcrit_filtered_orfs,
                dict_transcrit_filtered_orfs_total,
                ID_transcriptome,
            )
            implement_dict_pos_unspliced_ORFs_plus_stop(
                dict_pos_unspliced_ORFs_plus_stop,
                dict_pos_unspliced_ORFs_plus_stop_total,
                ID_transcriptome,
            )
            sublist_transcripts_overlap_total = implement_sublist_transcripts_overlap(
                sublist_transcripts_overlap,
                sublist_transcripts_overlap_total,
                ID_transcriptome,
            )
            implement_dict_all_ORFs_filtered_with_stop(
                dict_all_ORFs_filtered_with_stop,
                dict_all_ORFs_filtered_with_stop_total,
                ID_transcriptome,
            )
            implement_dict_transcript_fasta_denovo(
                dict_transcript_fasta_denovo,
                dict_transcript_fasta_denovo_total,
                ID_transcriptome,
            )

    if len(dict_gene_transcripts_denovo_total) > 0:
        # generate all output files (fasta + info file) of denovo candidates filtered after step 1.
        create_info_file_orfs(
            name_output_directory,
            dict_transcripts_properties_total,
            dict_gene_status_total,
            dict_gene_transcripts_denovo_total,
            dict_transcrit_filtered_orfs_total,
            dict_pos_unspliced_ORFs_plus_stop_total,
            sublist_transcripts_overlap_total,
        )
        create_denovo_nucl_file(
            name_output_directory, dict_all_ORFs_filtered_with_stop_total
        )
        create_denovo_prot_file(name_output_directory)
        create_denovo_transcripts_file(
            name_output_directory, dict_transcript_fasta_denovo_total
        )
        validated_step = True
    else:
        print(BRIGHT_MAGENTA + "No denovo transcript found ! " + RESET)
    return validated_step
