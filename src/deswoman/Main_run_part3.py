import os
from deswoman.module_generate_outgroup_dict import (
    get_dot_special_case_dic,
    create_outgroup_dict_strategy1,
)
from deswoman.module_blast_and_diamond import (
    perform_blast_to_genome,
    get_transcription_hits,
    perform_blast_prot_transcripts,
)
from deswoman.module_store_all_nc_hits import (
    make_dico_len_unspliced_denovo,
    extract_best_blast_outputs,
    extract_all_blast_outputs,
    store_de_novo_informations,
)
from deswoman.module_gffread_extract_and_rec_blast import (
    extract_all_prots_genome,
    make_reciprocal_prot_blasts,
    make_dico_blast_for_reciprocal_hits,
)
from deswoman.module_sort_gff_gene_order import store_genome_informations, sort_gff_dic
from deswoman.module_synteny import search_all_denovo_hits
from deswoman.module_prepare_global_alignment import (
    build_dico_seq_denovo_filtered_candidates,
    build_dico_seq_NcHomologs,
)
from deswoman.module_global_alignment_properties import (
    get_transcription_indication,
    main_alignment_function,
    create_final_output_file,
)
from deswoman.module_handle_strat3_hits import (
    build_target_and_query_prot_list,
    make_dico_name_denovo_strat2,
    handle_blast_inputs_transc,
    build_final_file_strat2_step3,
)
from deswoman.module_remove_extra_files import remove_intermediates
from deswoman.module_colors import *


def run_part3_strat1(dico_variables: dict) -> None:
    """
    Run Part 3 of Strategy 1: Homology search and synteny analysis for de novo sequences against target genomes.

    This function performs several key tasks in Part 3 of Strategy 1, which includes homology searches,
    reciprocal BLASTs, synteny analysis, and sequence alignments for de novo candidate sequences against
    target species genomes. It processes multiple species, extracting genomic and transcriptomic information,
    conducting BLAST searches, and aligning sequences based on synteny and homology criteria.

    Args:
        dico_variables (dict): A dictionary containing user-defined parameters, including:
            - "path_output": Directory where output files are saved.
            - "path_output_intermediate": Directory for intermediate output files.
            - "path_to_genome_repository": Path to the genome repository (for both query and outgroup genomes).
            - "query": Name of the query species for alignment.
            - "synteny_window": Window size for synteny analysis.
            - "premature_stop": Flag to handle premature stop codons during analysis.
            - "rm_undir_transc": Flag to remove undirected transcripts during processing.
            - "rec_best_hit": Best hit criteria for reciprocal BLAST (optional).
            - "dot_special_case": Special case for dot orientation (optional).

    Returns:
        None: This function modifies the output files directly and does not return a value.

    Workflow:
        1. **Homology Search**: Performs BLAST searches for de novo sequences against the target genome (DNA)
           and optionally transcriptome sequences.
        2. **Reciprocal BLAST**: Performs reciprocal BLAST searches between query and target genomes for protein hits.
        3. **Synteny Analysis**: Identifies syntenic hits between de novo sequences and the target genome.
        4. **Sequence Alignment and mutation detection**: Aligns de novo sequences to target genome and transcriptome sequences and detect and stores mutations results.
        5. **Output Generation**: Creates the final alignment file based on the results of homology and synteny analysis.

    Notes:
        - If no synteny criteria are provided (synteny_window = 0), only the best hits are considered.
        - If transcriptome data is available, it will be incorporated into the analysis.
        - If no syntenic hits are detected, a message is printed, and no final file is created.
    """
    # final dictionary where final output will be stored
    dico_format_all_data = {}
    # withdraw work directories
    name_output_directory = dico_variables["path_output"]
    name_intermediate_directory = dico_variables["path_output_intermediate"]
    # get dico with the sign of the dot
    dico_dot_special_case = get_dot_special_case_dic(name_intermediate_directory)
    # build dico with each target name associated to genome (fasta and gtf), and, if present, transcriptome
    dico_target_species_lines = create_outgroup_dict_strategy1(dico_variables)
    # generate path where to store blast output from step 3
    path_recip_blast = name_intermediate_directory + "/blast_denovo_to_target_genome"
    # create the folder
    if os.path.isdir(path_recip_blast) == False:
        os.system("mkdir " + path_recip_blast)
    # run then for each target species of the dictionnary
    for pop_species_name in dico_target_species_lines:
        print("")
        print(YELLOW + "-----------" + RESET)
        print(
            BRIGHT_BLUE
            + "Name target : "
            + RESET
            + BRIGHT_CYAN
            + pop_species_name
            + RESET
        )
        # create dico that will store best hits
        dico_denovo_best_hit = {}
        # make dico with denovo candidates as key and the len of the nucl seq with intron as item.
        dico_name_size_denovo = make_dico_len_unspliced_denovo(name_output_directory)
        # perform BLAST to genome and return whether the output file is empty or not
        check_empty = perform_blast_to_genome(
            pop_species_name,
            dico_target_species_lines[pop_species_name]["genome_fasta"],
            name_output_directory,
        )
        print(
            BRIGHT_BLUE
            + "Homology search of denovo candidate against target genome performed"
            + RESET
        )
        if check_empty != 0:
            # store best blast hits if no synteny critaria
            if dico_variables["synteny_window"] == 0:
                dico_denovo_best_hit = extract_best_blast_outputs(
                    pop_species_name, name_intermediate_directory
                )
            else:
                # store all blast hits (dico[namedenoo] : [namehit1, namehit2, etc])
                dict_my_blast_selected_outputs = extract_all_blast_outputs(
                    pop_species_name, name_intermediate_directory
                )
                # extracts necessary informations of all denovo candidates in info file
                dico_denovo_informations = store_de_novo_informations(
                    name_output_directory
                )
                # extract outgroup proteins and store them in the "old_prot" directory.
                extract_all_prots_genome(
                    pop_species_name,
                    dico_variables["path_to_genome_repository"],
                    name_intermediate_directory,
                )
                # same for prots of query genome
                extract_all_prots_genome(
                    dico_variables["query"],
                    dico_variables["path_to_genome_repository"],
                    name_intermediate_directory,
                )
                print(BRIGHT_BLUE + "Proteins extracted" + RESET)

                # perform reciprocal blast between query and target annotated proteins
                make_reciprocal_prot_blasts(
                    dico_variables["query"],
                    pop_species_name,
                    name_intermediate_directory,
                )
                # Store hits from query against target and target against query in 2 dicts.
                dico_blast1, dico_blast2 = make_dico_blast_for_reciprocal_hits(
                    dico_variables["rec_best_hit"], name_intermediate_directory
                )
                print(
                    BRIGHT_BLUE
                    + "Homology detection between annotated genes of query and target genomes performed"
                    + RESET
                )

                # extract genes informations and sort gff gene order for query and target
                # get dicos of all mRNAs corresponding to all extracted prots from target. dicos are : dico_gff[chromosome] = [[start,stop,gene_ID], etc].
                dico_gff_focal = store_genome_informations(
                    dico_variables["query"], dico_variables["path_to_genome_repository"]
                )
                # sort dico from the first to last in the genomic positions
                dico_gff_focal_sorted = sort_gff_dic(dico_gff_focal)
                # same two steps for outgroup genome
                dico_gff_outgroup = store_genome_informations(
                    pop_species_name, dico_variables["path_to_genome_repository"]
                )
                dico_gff_outgroup_sorted = sort_gff_dic(dico_gff_outgroup)

                ## extract syntenic hits. Makes a dict with the denovo as key and its best syntenic hit as item.
                dico_denovo_best_hit = search_all_denovo_hits(
                    dico_denovo_informations,
                    dict_my_blast_selected_outputs,
                    dico_gff_focal_sorted,
                    dico_gff_outgroup_sorted,
                    dico_variables["synteny_window"],
                    dico_blast1,
                    dico_blast2,
                )

            if len(dico_denovo_best_hit) > 0:
                ## Build dico of denovo as keys and their sequence (with intron if intron) as item
                dico_denovo_seq_intron = build_dico_seq_denovo_filtered_candidates(
                    dico_denovo_best_hit, name_output_directory
                )
                ## Build dico of denovo homologous seq that was selected in previous steps as keys and their sequence (with intron if intron) as item
                dico_NcHit = build_dico_seq_NcHomologs(
                    dico_denovo_best_hit,
                    dico_target_species_lines[pop_species_name]["genome_fasta"],
                    pop_species_name,
                )

                ## Align and retrieve infos
                # if transcriptome was provided, extracts its coordinates,chrom,orientation
                dico_target_transcript_coordinate = get_transcription_indication(
                    dico_target_species_lines,
                    pop_species_name,
                    dico_dot_special_case,
                    dico_variables["rm_undir_transc"],
                )
                # make a dico with de novo candidate and all target transcripts to whoch they have a nucl BLAST hit with the coverage associated.
                dico_target_transcript_blast_hits = get_transcription_hits(
                    dico_target_species_lines, pop_species_name, name_output_directory
                )
                # dico_format_all_data is implemented for each pop and each de novo, with all mutations observed.
                dico_format_all_data = main_alignment_function(
                    pop_species_name,
                    dico_name_size_denovo,
                    dico_denovo_best_hit,
                    dico_denovo_seq_intron,
                    dico_NcHit,
                    dico_target_transcript_coordinate,
                    dico_target_transcript_blast_hits,
                    dico_variables["premature_stop"],
                    dico_format_all_data,
                )
    if len(dico_format_all_data) > 0:
        # create final file.
        create_final_output_file(
            dico_format_all_data, name_output_directory, name_intermediate_directory
        )
    else:
        print(BRIGHT_MAGENTA + "No syntenic nchits detected ! " + RESET)
    remove_intermediates(name_intermediate_directory)


def run_part3_strat2(dico_variables: dict) -> None:
    """
    Run Part 3 of Strategy 2: Homology search and orthogroup analysis for de novo proteins against target transcriptomes.

    This function performs homology searches between de novo proteins and target transcriptomes,
    and groups proteins into orthogroups based on their sequence and genomic location similarity. It builds the necessary
    input files, performs BLAST searches, and organizes the resulting data into orthogroups for further analysis.

    Args:
        dico_variables (dict): A dictionary containing user-defined parameters, including:
            - "path_output": Directory where output files are saved.
            - "path_output_intermediate": Directory for intermediate output files.
            - "path_to_transcriptome_repository": Path to the transcriptome repository (for target species).
            - "query": Name of the query species for alignment (optional).
            - "synteny_window": Window size for synteny analysis (optional).
            - "premature_stop": Flag to handle premature stop codons during analysis (optional).
            - "rec_best_hit": Best hit criteria for reciprocal BLAST (optional).

    Returns:
        None: This function modifies the output files directly and does not return a value.

    Workflow:
        1. **Protein File Construction**: Constructs query and target protein files containing all de novo proteins.
        2. **BLAST Homology Search**: Performs BLAST searches to compare de novo proteins against the target transcriptome.
        3. **Orthogroup Assignment**: Groups proteins into orthogroups based on BLAST results.
        4. **Output Generation**: Creates the final output files containing the orthogroups and BLAST results.
    """
    name_output_directory = dico_variables["path_output"]
    name_intermediate_directory = dico_variables["path_output_intermediate"]
    dico_orthogroups = {}
    # create a query and target protein file with all de novo prots (same in both)
    build_target_and_query_prot_list(name_intermediate_directory, name_output_directory)
    # perform BLAST prot between query and target proteins (the two files are similar).
    perform_blast_prot_transcripts(name_intermediate_directory, name_output_directory)
    print(BRIGHT_BLUE + "ORFs homology between transcripts performed" + RESET)
    # build a dictionary with the name of the de novo prots.
    dico_name_size_denovo = make_dico_name_denovo_strat2(name_output_directory)
    # here, dico_orthogroups is implemented with the number of the orthogroup as a key and a list with all homologs as a target
    handle_blast_inputs_transc(
        dico_orthogroups,
        dico_name_size_denovo,
        name_intermediate_directory,
        name_output_directory,
    )
    build_final_file_strat2_step3(
        dico_orthogroups, name_intermediate_directory, name_output_directory
    )
    remove_intermediates(name_intermediate_directory)
