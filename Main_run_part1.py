
from Bio.SeqRecord import SeqRecord
import os
from module_assess_input import *
from module_transcripts_data_and_threeshold import *
#from module_transcripts_overlap import *
#from module_orfs_old import * ## see
from module_orfs import * ## see
from module_unspliced_orfs import *
from module_lower_intron import *
from module_output_part1 import *
from module_my_bedtool import *
from module_reformat_dico_with_IDs import *

cwd = os.getcwd()




## Assess and store variables given by the user
Strategy = openFile("Strategy.setup")
dico_variables = withdraw_strategy_data(Strategy)



if dico_variables["strategy"] == "1":
# Part 1 : 
# -Retrieve all transcripts from data, get orfs, sort orfs by removing the ones that end with the end of transcript and the ones that are reverse

    
    link_to_query_gtf = get_path_gtf(dico_variables["path_to_transcriptome_repository"], dico_variables["query"])
    link_to_query_transcriptome = get_path_to_fasta(dico_variables["path_to_transcriptome_repository"],dico_variables["query"])
    link_to_query_genome_gff = get_path_gff(dico_variables["path_to_genome_repository"], dico_variables["query"])
    link_to_query_genome = get_path_to_fasta(dico_variables["path_to_genome_repository"],dico_variables["query"])

    dict_transcript_fasta_all = store_transcriptome(link_to_query_transcriptome)
    gtf_file = openFile(link_to_query_gtf) # Note I have also to deal with GFF files
    dict_gene_transcripts_all = associate_gene_and_transcripts(gtf_file,dico_variables["TPM_threeshold"], dict_transcript_fasta_all, dico_variables["filter_TE"])
    dict_transcripts_properties = transcripts_properties(gtf_file)

    ## Part 2 get transcript overlap to genome elements and withdraw de novo transcripts
    dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos = my_bedtool_genes_introns_ordered(link_to_query_genome_gff)
    sublist_transcripts_overlap = my_bedtool_extract_overlap(dict_transcripts_properties,dict_gene_transcripts_all, dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos)
    dict_gene_transcripts_denovo,dict_transcript_fasta_denovo,dict_transcrip_status,dict_gene_status = extract_denovo_transcripts(sublist_transcripts_overlap, dico_variables["transcript_overlap"],dict_transcript_fasta_all,dict_gene_transcripts_all)


    # Part3 : retrieve orfs from de novo transcripts and reduce them by properties
    #dict_all_ORFs = build_orf_file_rename_discard_reverse("denovo_transcript.fa",dico_variables["query"]) # with old
    dict_all_ORFs = my_get_orfs("denovo_transcript.fa") # with new
    dict_all_ORFs_purge1 = remove_orfs_end_transcript(dict_all_ORFs,dict_transcript_fasta_denovo)
    dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_all_ORFs_purge1)
    list_gene_object = generate_list_genes_objects(dict_gene_transcripts_denovo,dict_transcrit_all_orfs,dict_all_ORFs_purge1,dict_transcript_fasta_denovo)

    # option_list can include [["kozac_threeshold",5000], ["kozac_highest"], ["longest"], ["start_first"], ["utr_size,50,50"], ["duplicate_handle"], ["gene_denovo"]]
    dict_transcrit_filtered_orfs,dict_all_ORFs_filtered = sort_orfs_by_properties(dico_variables["filter_genic"], list_gene_object, dico_variables["ORFs_choice"],dict_all_ORFs_purge1, dict_gene_status)                


    # Part 4 : retrieve unspliced seqs of orfs and control steps
    dict_gene_exon_pos = build_dict_unspliced_seqs(gtf_file)
    print ("Total number of ORFs in denovo transcripts : " + str(len(dict_all_ORFs)))
    print ("Total number of correct sense ORFs in denovo transcripts : " + str(len(dict_all_ORFs_purge1)))
    # Here, very carefull. I have to investigate further this step, because according to the assembler that has been used for the transcriptome assembly, the transcript can be extracted directly from the genome , and therefore correspond completely to the genome, or like with trinity correspond to the rna seq data, can can show small differences to the genome. https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly
    # Therefore, i can not reconstruct the orf from the genome directly if it does not correspond, or more precisely, i have to extract the stop codon directly to the rna.
    dict_pos_unspliced_ORFs_plus_stop = rebuild_spliced_orfs(link_to_query_genome, dict_gene_exon_pos, dict_all_ORFs_filtered)
    dict_all_ORFs_filtered_with_stop = add_stop_codon_to_spliced_orfs(dict_pos_unspliced_ORFs_plus_stop, dict_all_ORFs_filtered)
    print ("Total number of correct sense ORFs in denovo transcripts filtered : " + str(len(dict_all_ORFs_filtered_with_stop)))

    # Part 5 : build unspliced orfs with intron
    #dict_intron_exon_pos, dict_transcript_chrom = build_dico_intron_exon_specified(dict_gene_exon_pos)
    print ("lala")
    dict_unspliced_orfs_fasta_lowered_introns = build_unspliced_fasta_seq_lowered_intron(dict_pos_unspliced_ORFs_plus_stop, link_to_query_genome, dict_gene_exon_pos)
    print ("lili")
    # Part 6 : build output 
    create_new_folder()
    create_unspliced_orfs(dict_unspliced_orfs_fasta_lowered_introns)
    create_info_file_orfs(dict_transcripts_properties, dict_gene_status, dict_gene_transcripts_denovo, dict_transcrit_filtered_orfs,dict_pos_unspliced_ORFs_plus_stop,sublist_transcripts_overlap)
    create_denovo_nucl_file(dict_all_ORFs_filtered_with_stop)
    create_denovo_prot_file()
    create_denovo_transcripts_file(dict_transcript_fasta_denovo)
    ## Note need getORF, biopython, KCS-predicted.tsv, bedtoools stored in the directory

if dico_variables["strategy"] == "2":

    ## save total final output 
    dict_transcripts_properties_total = {}
    dict_gene_status_total = {}
    dict_gene_transcripts_denovo_total = {}
    dict_transcrit_filtered_orfs_total = {}
    dict_pos_unspliced_ORFs_plus_stop_total = {}
    sublist_transcripts_overlap_total = []
    dict_all_ORFs_filtered_with_stop_total = {}
    dict_transcript_fasta_denovo_total = {}
    ##
    link_to_query_genome_gff = get_path_gff(dico_variables["path_to_genome_repository"], dico_variables["query"])
    link_to_query_genome = get_path_to_fasta(dico_variables["path_to_genome_repository"],dico_variables["query"])
    list_query_name_transc = get_list_name_query_rep(dico_variables["path_to_transcriptome_repository"])
    dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos = my_bedtool_genes_introns_ordered(link_to_query_genome_gff)
    for ID_transcriptome in list_query_name_transc:
        print (ID_transcriptome)
        link_to_query_gtf = get_path_gtf(dico_variables["path_to_transcriptome_repository"], ID_transcriptome)
        link_to_query_transcriptome = get_path_to_fasta(dico_variables["path_to_transcriptome_repository"],ID_transcriptome)
        dict_transcript_fasta_all = store_transcriptome(link_to_query_transcriptome)
        gtf_file = openFile(link_to_query_gtf) # Note I have also to deal with GFF files
        dict_gene_transcripts_all = associate_gene_and_transcripts(gtf_file,dico_variables["TPM_threeshold"], dict_transcript_fasta_all, dico_variables["filter_TE"])
        dict_transcripts_properties = transcripts_properties(gtf_file)
        ## Part 2 get transcript overlap to genome elements and withdraw de novo transcripts
        #build_bed_files_reference(link_to_query_genome_gff)
        #sublist_transcripts_overlap = get_transcripts_overlap(dict_transcripts_properties,dict_gene_transcripts_all)
        
        sublist_transcripts_overlap = my_bedtool_extract_overlap(dict_transcripts_properties,dict_gene_transcripts_all, dico_sorted_chrom_exon_pos, dico_sorted_chrom_intron_pos)
        dict_gene_transcripts_denovo,dict_transcript_fasta_denovo,dict_transcrip_status,dict_gene_status = extract_denovo_transcripts(sublist_transcripts_overlap, dico_variables["transcript_overlap"],dict_transcript_fasta_all,dict_gene_transcripts_all)
       
        # Part3 : retrieve orfs from de novo transcripts and reduce them by properties
        dict_all_ORFs = my_get_orfs("denovo_transcript.fa") # with new
        dict_all_ORFs_purge1 = remove_orfs_end_transcript(dict_all_ORFs,dict_transcript_fasta_denovo)
        dict_transcrit_all_orfs = build_dict_orfs_per_transcript(dict_all_ORFs_purge1)
        list_gene_object = generate_list_genes_objects(dict_gene_transcripts_denovo,dict_transcrit_all_orfs,dict_all_ORFs_purge1,dict_transcript_fasta_denovo)

        # option_list can include [["kozac_threeshold",5000], ["kozac_highest"], ["longest"], ["start_first"], ["utr_size,50,50"], ["duplicate_handle"], ["gene_denovo"]]
        dict_transcrit_filtered_orfs,dict_all_ORFs_filtered = sort_orfs_by_properties(dico_variables["filter_genic"], list_gene_object, dico_variables["ORFs_choice"],dict_all_ORFs_purge1, dict_gene_status)                


        # Part 4 : retrieve unspliced seqs of orfs and control steps
        dict_gene_exon_pos = build_dict_unspliced_seqs(gtf_file)
        print ("Total number of ORFs in denovo transcripts : " + str(len(dict_all_ORFs)))
        print ("Total number of correct sense ORFs in denovo transcripts : " + str(len(dict_all_ORFs_purge1)))
        # Here, very carefull. I have to investigate further this step, because according to the assembler that has been used for the transcriptome assembly, the transcript can be extracted directly from the genome , and therefore correspond completely to the genome, or like with trinity correspond to the rna seq data, can can show small differences to the genome. https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly
        # Therefore, i can not reconstruct the orf from the genome directly if it does not correspond, or more precisely, i have to extract the stop codon directly to the rna.
        dict_pos_unspliced_ORFs_plus_stop = rebuild_spliced_orfs(link_to_query_genome, dict_gene_exon_pos, dict_all_ORFs_filtered)
        dict_all_ORFs_filtered_with_stop = add_stop_codon_to_spliced_orfs(dict_pos_unspliced_ORFs_plus_stop, dict_all_ORFs_filtered)
        print ("Total number of correct sense ORFs in denovo transcripts filtered : " + str(len(dict_all_ORFs_filtered_with_stop)))

        # Part 6 : build output 

        implement_dict_transcripts_properties(dict_transcripts_properties, dict_transcripts_properties_total, ID_transcriptome)
        implement_dict_gene_status(dict_gene_status, dict_gene_status_total, ID_transcriptome)
        implement_dict_gene_transcripts_denovo(dict_gene_transcripts_denovo, dict_gene_transcripts_denovo_total, ID_transcriptome)
        implement_dict_transcrit_filtered_orfs(dict_transcrit_filtered_orfs, dict_transcrit_filtered_orfs_total, ID_transcriptome)
        implement_dict_pos_unspliced_ORFs_plus_stop(dict_pos_unspliced_ORFs_plus_stop, dict_pos_unspliced_ORFs_plus_stop_total, ID_transcriptome)
        sublist_transcripts_overlap_total = implement_sublist_transcripts_overlap(sublist_transcripts_overlap, sublist_transcripts_overlap_total, ID_transcriptome)
        print (len(sublist_transcripts_overlap_total[0]))
        
        
        implement_dict_all_ORFs_filtered_with_stop(dict_all_ORFs_filtered_with_stop, dict_all_ORFs_filtered_with_stop_total, ID_transcriptome)
        implement_dict_transcript_fasta_denovo(dict_transcript_fasta_denovo, dict_transcript_fasta_denovo_total, ID_transcriptome)

    create_new_folder()
    create_info_file_orfs(dict_transcripts_properties_total, dict_gene_status_total, dict_gene_transcripts_denovo_total, dict_transcrit_filtered_orfs_total,dict_pos_unspliced_ORFs_plus_stop_total,sublist_transcripts_overlap_total)
    create_denovo_nucl_file(dict_all_ORFs_filtered_with_stop_total)
    create_denovo_prot_file()
    create_denovo_transcripts_file(dict_transcript_fasta_denovo_total)


print ("DESMAN Part 1 : DONE!")
