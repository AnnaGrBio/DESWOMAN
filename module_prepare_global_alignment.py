import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq




def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 

def build_dico_seq(dico_denovo_best_hit):
    dico_denovo_seq_intron = {}
    for seq_record in SeqIO.parse("DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa", "fasta"):
        if str(seq_record.id) in dico_denovo_best_hit.keys():
            dico_denovo_seq_intron[str(seq_record.id)] = str(seq_record.seq)
    return dico_denovo_seq_intron 


def build_dico_seq_NcHomologs(dico_denovo_best_hit, link_to_genome, pop_name):
    dico_NcHit = {}
    dico_genome = {}
    for seq_record in SeqIO.parse(link_to_genome, "fasta"):
        dico_genome[str(seq_record.id)] = str(seq_record.seq)
    for denovo in dico_denovo_best_hit.keys():
        NcHit_name = denovo # + "-" + pop_name + "-Nc"
        hit_coords = dico_denovo_best_hit[denovo]
        chrom_scaf = hit_coords.split("-")[2]
        start = int(hit_coords.split("-")[0])-2
        stop = int(hit_coords.split("-")[1])+2
        seq_hit_extended = dico_genome[chrom_scaf][start:stop]
        #print (seq_hit_extended)
        if hit_coords.split("-")[3] == "r":
            seq_hit_extended = Seq(seq_hit_extended)
            seq_hit_extended_rev_comp = str(seq_hit_extended.reverse_complement())
            dico_NcHit[NcHit_name] = seq_hit_extended_rev_comp
        else:
            dico_NcHit[NcHit_name] = seq_hit_extended
    return dico_NcHit
