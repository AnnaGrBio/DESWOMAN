import os
#os.chdir("/home/agrandch/Desktop/postdoc/Projet_principal/DESMAN/Get_and_sort_transcripts_ORFs")
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def create_new_folder():
    os.system ("mkdir DESMAN_denovo_output")
    os.system("rm denovo_transcript.fa")
    
def create_unspliced_orfs(dict_unspliced_orfs_fasta_lowered_introns):
    my_file = open("denovo_unspliced_lowered_introns.fa", "w")
    for name in dict_unspliced_orfs_fasta_lowered_introns.keys():
        my_file.write(">"+name+"\n")
        my_file.write(dict_unspliced_orfs_fasta_lowered_introns[name]+"\n")
    my_file.close()
    os.system ("mv denovo_unspliced_lowered_introns.fa DESMAN_denovo_output")
    
    
def create_denovo_nucl_file(dict_all_ORFs_filtered_with_stop):
    my_file = open("denovo_nucl.fa", "w")
    for name in dict_all_ORFs_filtered_with_stop.keys():
        my_file.write(">"+name+"\n")
        my_file.write(dict_all_ORFs_filtered_with_stop[name].upper()+"\n")
    my_file.close()

def create_denovo_transcripts_file(dict_transcript_fasta_denovo):
    my_file = open("denovo_transcripts.fa", "w")
    for name in dict_transcript_fasta_denovo.keys():
        my_file.write(">"+name+"\n")
        my_file.write(dict_transcript_fasta_denovo[name].upper()+"\n")
    my_file.close()
    os.system ("mv denovo_transcripts.fa DESMAN_denovo_output")
    
def create_denovo_prot_file():
    my_file = open("denovo_protein.fa", "w")
    for record in SeqIO.parse("denovo_nucl.fa", "fasta"):
        protein_id = record.id
        protein = str(record.seq.translate(to_stop=True))
        my_file.write(">"+str(protein_id)+"\n")
        my_file.write(protein+"\n")
    my_file.close()
    os.system ("mv denovo_nucl.fa DESMAN_denovo_output")
    os.system ("mv denovo_protein.fa DESMAN_denovo_output")
        
    
def create_info_file_orfs(dict_transcripts_properties, dict_gene_status, dict_gene_transcripts_denovo, dict_transcrit_filtered_orfs,dict_pos_unspliced_ORFs_plus_stop,sublist_transcripts_overlap):
    my_file = open("information_file.txt", "w")
    my_file.write("gene_name,chrom,direction,gene_status,nb_denovo_transcripts,transcript_name,start_in_genome,end_in_genome,nb_filtered_orf,orf_name,start_in_genome,end_in_genome,genomic_overlap"+"\n")
    for gene_name in dict_gene_transcripts_denovo.keys():
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
            if transcript_name in dict_transcrit_filtered_orfs.keys():
                list_orf = dict_transcrit_filtered_orfs[transcript_name]
                for orf_name in list_orf:
                    if orf_name in dict_pos_unspliced_ORFs_plus_stop.keys():
                        my_gene_name = gene_name
                        my_gene_status = dict_gene_status[my_gene_name]
                        my_transcript_name = transcript_name
                        my_chrom = dict_transcripts_properties[my_transcript_name][0]
                        my_direction = dict_transcripts_properties[my_transcript_name][1]
                        my_nb_denovo_transcripts = str(len(list_transcripts))
                        my_transcript_start = str(dict_transcripts_properties[my_transcript_name][2])
                        my_transcript_stop = str(dict_transcripts_properties[my_transcript_name][3])
                        my_nb_filtered_orfs = str(len(list_orf))
                        my_orf_name = orf_name
                        my_orf_start_in_genome = str(dict_pos_unspliced_ORFs_plus_stop[my_orf_name][0])
                        my_orf_stop_in_genome = str(dict_pos_unspliced_ORFs_plus_stop[my_orf_name][1])
                        my_file.write(my_gene_name+","+my_chrom+","+my_direction+","+my_gene_status+","+my_nb_denovo_transcripts+","+my_transcript_name+","+my_transcript_start+","+my_transcript_stop+","+my_nb_filtered_orfs+","+my_orf_name+","+my_orf_start_in_genome+","+my_orf_stop_in_genome+","+transcript_overlap+"\n")
    my_file.close()
    os.system ("mv information_file.txt DESMAN_denovo_output")

def create_info_file_transcripts(dict_transcript_fasta_denovo, dict_gene_status, dict_transcripts_properties, sublist_transcripts_overlap):
    my_file = open("information_file.txt", "w")
    my_file.write("gene_name,chrom,direction,gene_status,transcript_name,start_in_genome,end_in_genome,genomic_overlap"+"\n")
    for transcript_name in dict_transcript_fasta_denovo.keys():
        gene_name = transcript_name.split(".")[0] + "." + transcript_name.split(".")[1]
        elts_transcripts = dict_transcripts_properties[transcript_name]
        chrom = elts_transcripts[0]
        direction = elts_transcripts[1]
        start = elts_transcripts[2]
        stop = elts_transcripts[3]
        gene_status = dict_gene_status[gene_name]
        transcript_overlap = ""
        if transcript_name in sublist_transcripts_overlap[0]:
            transcript_overlap = "intergenic_transcripts"
        elif transcript_name in sublist_transcripts_overlap[1]:
            transcript_overlap = "intronic_transcripts"
        elif transcript_name in sublist_transcripts_overlap[2]:
            transcript_overlap = "antisense_transcript"
        else:
            transcript_overlap = "genic_transcripts"
        my_file.write(gene_name+","+chrom+","+direction+","+gene_status+","+transcript_name+","+start+","+stop+","+transcript_overlap+"\n")
    my_file.close()
    os.system ("mv information_file.txt DESMAN_denovo_output")





                    
def renew_gtf_file(): # To BUild!!!!! Note                
    pass





