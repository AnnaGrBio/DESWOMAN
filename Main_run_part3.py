import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_assess_input import *
from module_generate_outgroup_dict import *
from module_blast_and_diamond import *
from module_store_all_nc_hits import *
from module_gffread_extract_and_rec_blast import *
from module_sort_gff_gene_order import *
from module_synteny import *
from module_prepare_global_alignment import *
from module_global_alignment_properties import *
from module_handle_strat3_hits import *

cwd = os.getcwd()


## Assess and store variables given by the user
Strategy = openFile("Strategy.setup")
dico_variables = withdraw_strategy_data(Strategy)
dico_format_all_data = {}

if dico_variables["strategy"] == "1":
    dico_target_species_lines = create_outgroup_dict_strategy1(dico_variables)
    c = 0
    os.system ("mkdir blast_denovo_to_target_genome")
    for pop_species_name in dico_target_species_lines.keys():
        print (pop_species_name)
        dico_denovo_best_hit = {}
        dico_name_size_denovo = make_dico_len_unspliced_denovo()
        perform_blast_to_genome(pop_species_name, dico_target_species_lines[pop_species_name]["genome_fasta"])
        if dico_variables["synteny_window"] == 0:
            # store best blast hits
            dico_denovo_best_hit = extract_best_blast_outputs(pop_species_name)
        else:
            # store all blast hits
            dict_my_blast_selected_outputs = extract_all_blast_outputs(pop_species_name)
            
            dico_denovo_informations = store_de_novo_informations()
            
            # extract outgroup proteins
            extract_all_prots_genome(pop_species_name, dico_variables["path_to_genome_repository"])
            print ("lala")
            extract_all_prots_genome(dico_variables["query"], dico_variables["path_to_genome_repository"])

            # perform reciprocal blast
            make_reciprocal_prot_blasts(dico_variables["query"], pop_species_name)
            dico_blast1, dico_blast2 = make_dico_blast()

            # extract genes and sort gff gene order for query and target
            dico_gff_focal = store_genome_informations(dico_variables["query"], dico_variables["path_to_genome_repository"])
            dico_gff_focal_sorted = sort_gff_dic(dico_gff_focal)
            dico_gff_outgroup = store_genome_informations(pop_species_name, dico_variables["path_to_genome_repository"])
            dico_gff_outgroup_sorted = sort_gff_dic(dico_gff_outgroup)

            ## extract syntenic hits
            dico_denovo_best_hit = search_all_denovo_hits(dico_denovo_informations,dict_my_blast_selected_outputs,dico_gff_focal_sorted, dico_gff_outgroup_sorted, dico_variables["synteny_window"], dico_blast1, dico_blast2)

            ## Build sequences files (denovo and best hit)
        dico_denovo_seq_intron = build_dico_seq(dico_denovo_best_hit)
        dico_NcHit = build_dico_seq_NcHomologs(dico_denovo_best_hit, dico_target_species_lines[pop_species_name]["genome_fasta"], pop_species_name)

        ## Align and retrieve infos
        dico_target_transcript_coordinate = get_transcription_indication(dico_target_species_lines, pop_species_name)
        dico_target_transcript_blast_hits = get_transcription_hits(dico_target_species_lines, pop_species_name)
        dico_format_all_data = main_alignment_function(pop_species_name, dico_name_size_denovo, dico_denovo_best_hit, dico_denovo_seq_intron, dico_NcHit, dico_target_transcript_coordinate, dico_target_transcript_blast_hits, dico_variables["premature_stop"], dico_format_all_data)
    create_final_big_file(dico_format_all_data)

elif dico_variables["strategy"] == "2":
    list_query_name_transc = get_list_name_query_rep(dico_variables["path_to_transcriptome_repository"])
    dico_orthogroups = {}
    for ID_transcriptome in list_query_name_transc:
        build_target_prot_list(ID_transcriptome)
        perform_blast_prot_transcripts()
        dico_denovo_prot_hits = handle_blast_inputs_transc()
        dico_denovo_prot_hits_reduced = reduce_by_location(dico_denovo_prot_hits)
        fill_dico_orthogroups(dico_denovo_prot_hits_reduced, dico_orthogroups)
        build_final_file_strat2(dico_orthogroups)
        




    

