import os
from Bio import SeqIO



def implement_dict_transcripts_properties(dict_transcripts_properties, dict_transcripts_properties_total, ID_transcriptome):
    for transcript_name in dict_transcripts_properties.keys():
        new_transcript_name = ID_transcriptome+"_" + transcript_name
        dict_transcripts_properties_total[new_transcript_name] = dict_transcripts_properties[transcript_name]


def implement_dict_gene_status(dict_gene_status, dict_gene_status_total, ID_transcriptome):
    for gene_name in dict_gene_status.keys():
        new_gene_name = ID_transcriptome + "_" + gene_name
        dict_gene_status_total[new_gene_name] = dict_gene_status[gene_name]


def implement_dict_gene_transcripts_denovo(dict_gene_transcripts_denovo, dict_gene_transcripts_denovo_total, ID_transcriptome):
    for gene_name in dict_gene_transcripts_denovo.keys():
        new_gene_name = ID_transcriptome + "_" + gene_name
        new_list_transcs = []
        for transcript in dict_gene_transcripts_denovo[gene_name]:
            new_transcript_name = ID_transcriptome + "_" + transcript
            new_list_transcs.append(new_transcript_name)
        dict_gene_transcripts_denovo_total[new_gene_name] = new_list_transcs


def implement_dict_transcrit_filtered_orfs(dict_transcrit_filtered_orfs, dict_transcrit_filtered_orfs_total, ID_transcriptome):
    for transcript_name in dict_transcrit_filtered_orfs.keys():
        new_transcript_name = ID_transcriptome + "_" + transcript_name
        list_new_orfs_names = []
        for orf in dict_transcrit_filtered_orfs[transcript_name]:
            new_orf_name = ID_transcriptome + "_" + orf
            list_new_orfs_names.append(new_orf_name)
        dict_transcrit_filtered_orfs_total[new_transcript_name] = list_new_orfs_names


def implement_dict_pos_unspliced_ORFs_plus_stop(dict_pos_unspliced_ORFs_plus_stop, dict_pos_unspliced_ORFs_plus_stop_total, ID_transcriptome):
    for orf_name in dict_pos_unspliced_ORFs_plus_stop.keys():
        new_orf_name = ID_transcriptome + "_" + orf_name
        dict_pos_unspliced_ORFs_plus_stop_total[new_orf_name] = dict_pos_unspliced_ORFs_plus_stop[orf_name]


def implement_sublist_transcripts_overlap(sublist_transcripts_overlap, sublist_transcripts_overlap_total, ID_transcriptome):
    if len(sublist_transcripts_overlap_total) == 0:
        sublist_transcripts_overlap_total = [[], [], [], []]
    for transcript in sublist_transcripts_overlap[0]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[0].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[1]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[1].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[2]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[2].append(new_transcript_name)
    for transcript in sublist_transcripts_overlap[3]:
        new_transcript_name = ID_transcriptome + "_" + transcript
        sublist_transcripts_overlap_total[3].append(new_transcript_name)
    return sublist_transcripts_overlap_total


def implement_dict_all_ORFs_filtered_with_stop(dict_all_ORFs_filtered_with_stop, dict_all_ORFs_filtered_with_stop_total, ID_transcriptome):
    for orf_name in dict_all_ORFs_filtered_with_stop.keys():
        new_orf_name = ID_transcriptome + "_" + orf_name
        dict_all_ORFs_filtered_with_stop_total[new_orf_name] = dict_all_ORFs_filtered_with_stop[orf_name]


def implement_dict_transcript_fasta_denovo(dict_transcript_fasta_denovo, dict_transcript_fasta_denovo_total, ID_transcriptome):
    for transcript_name in dict_transcript_fasta_denovo.keys():
        new_transcript_name = ID_transcriptome + "_" + transcript_name
        dict_transcript_fasta_denovo_total[new_transcript_name] = dict_transcript_fasta_denovo[transcript_name]