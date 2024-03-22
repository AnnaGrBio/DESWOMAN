import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import MultipleSeqAlignment
from Bio.pairwise2 import format_alignment



def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def chech_ali(seq):
    for numbers in range(0,len(seq)):
        if seq[numbers] != "-":
            break
    start = seq[numbers:numbers+3].upper()
    for numbers in range(len(seq)-1,0,-1):
        if seq[numbers] != "-":
            break
    stop = seq[numbers-2:numbers+1].upper()
    if start == "ATG" and stop == "TAG" or start == "ATG" and stop == "TAA" or start == "ATG" and stop == "TGA":
        #print (start)
        return True
        
    else:
        return False
    

def searchATG(alignment):
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
        return "A"
    else:
        return "P"
    

def searchStop(alignment):
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
        #print (stop_nchit)
        return "A"
    else:
        return "P"
    

def searchIndels(alignment):
    denovo = alignment[0]
    nchit = alignment[1]
    nb_gap_denovo = 0
    for numbers_start in range(0,len(denovo)):
        if denovo[numbers_start] != "-":
            break
    for numbers_stop in range(len(denovo)-1,0,-1):
        if denovo[numbers_stop] != "-":
            break
    seq_denovo_cut = denovo[numbers_start:numbers_stop+1]
    seq_nchit_cut = nchit[numbers_start:numbers_stop+1]
    nb_gap_denovo = seq_denovo_cut.count('-')
    nb_gap_nchit = seq_nchit_cut.count('-')
    indels = nb_gap_denovo - nb_gap_nchit
    if indels%3 == 0:
        frameshift = "A"
    else:
        frameshift = "P"
    return indels, frameshift


def searchSubs(alignment):
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


def searchPreStop(alignment, perc_seq_accepted):
    list_stop = ["TAG", "TAA", "TGA"]
    denovo = alignment[0]
    nchit = alignment[1]
    limit_search = int(len(nchit) * perc_seq_accepted / 100)
    nb_gap_nchit = 0
    premature_stop = ""
    pos_premature_stop = 0
    for numbers_start in range(0,len(denovo)):
        if denovo[numbers_start] != "-":
            break
    for numbers_stop in range(len(denovo)-1,0,-1):
        if denovo[numbers_stop] != "-":
            break
    triplet_denovo = ""
    triplet_nchit = ""
    for i in range(numbers_start, limit_search,3):
        triplet_nchit = nchit[i:i+3].upper()
        nb_gap_nchit += triplet_nchit.count('-')
        if triplet_nchit in list_stop :
            triplet_denovo = denovo[i:i+3].upper()
            if triplet_denovo not in list_stop:
                premature_stop = triplet_nchit
                pos_premature_stop = i
                break
    if premature_stop == "":
        prem_stop_val = "A"
    else:
        prem_stop_val = "P"
    return prem_stop_val, pos_premature_stop

        

    
def measure_size_intron(denovo_seq):
    size_intron = 0
    for i in denovo_seq:
        if i == "a" or i == "t" or i == "g" or i == "c":
            size_intron += 1
    return size_intron


def splice_alignment(alignment):
    denovo = alignment[0]
    nchit = alignment[1]
    new_seq_denovo = ""
    new_seq_hit = ""
    list_introns = []
    seq_intron = ""
    for number in range(0,len(denovo)):
        nucl_denovo = denovo[number]
        if nucl_denovo != "a" and nucl_denovo != "t" and nucl_denovo != "g" and nucl_denovo != "c" and nucl_denovo != "n":
            if seq_intron != "":
                list_introns.append(seq_intron)
                seq_intron = ""
            new_seq_denovo += nucl_denovo
            new_seq_hit += nchit[number]
        else:
            seq_intron += nucl_denovo
    if seq_intron != "":
        list_introns.append(seq_intron)
    alignment = [new_seq_denovo,new_seq_hit,alignment[2],alignment[3],alignment[4]]
    return alignment,list_introns


def search_transcript_overlap(dico_target_transcript_coordinate,start_nc,stop_nc,direction_nc):
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
    return value


def validateTranscription(name_denovo, name_nchom, dico_target_transcript_coordinate, dico_target_transcript_blast_hits):
    valueTranscription = "NA"
    if len(dico_target_transcript_coordinate) > 0:
        valueTranscription = "A"
        list_elts_nchom = name_nchom.split("-")
        start = int(list_elts_nchom[0])
        stop = int(list_elts_nchom[1])
        chrom = list_elts_nchom[2]
        direction = list_elts_nchom[3]
        sub_dico_transcripts = dico_target_transcript_coordinate[chrom]
        for transcript in sub_dico_transcripts.keys():
            transc_status = search_transcript_overlap(sub_dico_transcripts[transcript],start,stop,direction)
            if transc_status == "reverse":
                valueTranscription = transc_status
                break
            elif transc_status == "forward":
                if name_denovo in dico_target_transcript_blast_hits.keys():
                    dico_denovo_hit = dico_target_transcript_blast_hits[name_denovo]
                    if transcript in dico_denovo_hit.keys():
                        valueTranscription = dico_denovo_hit[transcript]
                        break
    if valueTranscription == "NA":
        print (valueTranscription)
    return valueTranscription


def searchSpliSite():
    return "NA"


def searchMotifs():
    return "NA"


def main_alignment_function(pop_species_name, dico_name_size_denovo, dico_denovo_best_hit, dic_denovo, dico_NcHit, dico_target_transcript_coordinate, dico_target_transcript_blast_hits, value_prem_stop, dico_format_all_data):
    score = 0
    for denovo_name in dico_name_size_denovo.keys():
        if denovo_name in dico_NcHit.keys() and denovo_name in dic_denovo.keys():
            validated_hit = "P"
            #print(denovo_name)
            #print (dico_denovo_best_hit[denovo_name])
            denovo_seq = dic_denovo[denovo_name]
            nchit_seq = dico_NcHit[denovo_name].upper()
            correct_ali = False
            intron = False
            alignment = False
            if "a" not in denovo_seq and "t" not in denovo_seq and "g" not in denovo_seq and "c" not in denovo_seq:
                alignment = pairwise2.align.globalms(denovo_seq, nchit_seq, 2, -1, -1.5, -.1)
            else:
                size_intron = measure_size_intron(denovo_seq)
                if size_intron < 3000:
                    alignment = pairwise2.align.globalms(denovo_seq, nchit_seq, 2, -1, -1.5, -.1)
                    intron = True
            if alignment != False:
                for ali in alignment:
                    correct_ali = chech_ali(ali[0])
                    if correct_ali == True:
                        alignment = ali
                        break
                if correct_ali == False:
                    alignment = alignment[0]
                    score += 1
                if intron == True:
                    #print(format_alignment(*alignment))
                    alignment,list_introns = splice_alignment(alignment)
                    #if len(list_introns) > 1:
                        #print (len(list_introns))
                        #print (denovo_seq)
                        #print ("**********")
                    #print(format_alignment(*alignment))
                    
                ATG = searchATG(alignment)
                stop = searchStop(alignment)
                indels, frameshift = searchIndels(alignment)
                nb_subs = searchSubs(alignment)
                premature_stop, pos_premature_stop = searchPreStop(alignment,value_prem_stop)
                transcription = validateTranscription(denovo_name, dico_denovo_best_hit[denovo_name], dico_target_transcript_coordinate, dico_target_transcript_blast_hits)
                splicingSite = searchSpliSite()
                motifs = searchMotifs()
                #print(format_alignment(*alignment))
            else:
                ATG = "NA"
                stop = "NA"
                indels, frameshift = "NA","NA"
                nb_subs = "NA"
                premature_stop, pos_premature_stop = "NA","NA"
                transcription = "NA"
                splicingSite = "NA"
                motifs = "NA"
        else:
            validated_hit = "A"
            ATG = "NA"
            stop = "NA"
            indels, frameshift = "NA","NA"
            nb_subs = "NA"
            premature_stop, pos_premature_stop = "NA","NA"
            transcription = "NA"
            splicingSite = "NA"
            motifs = "NA"

        list_properties = [validated_hit, ATG, stop, indels, frameshift, nb_subs, premature_stop, pos_premature_stop, transcription, splicingSite, motifs]
        if denovo_name not in dico_format_all_data.keys():
            dico_format_all_data[denovo_name] = {pop_species_name:list_properties}
        else:
            dico_format_all_data[denovo_name][pop_species_name] = list_properties
    return dico_format_all_data
            
   #     
    #print (score)


def get_transcription_indication(dico_target_species_lines, target_name):
    dico_transcript_position = {}
    if "transcriptome_gtf" in dico_target_species_lines[target_name].keys():
        my_gtf_file = openFile(dico_target_species_lines[target_name]["transcriptome_gtf"])
        for line in my_gtf_file:
            if line[0] != "#":
                elts_line = line.split()
                if elts_line[2] == "transcript":
                    chrom = elts_line[0]
                    start = int(elts_line[3])
                    stop = int(elts_line[4])
                    signe_direction = elts_line[6]
                    if signe_direction == "+":
                        direction = "f"
                    else:
                        direction = "r"
                    transcript_name = elts_line[11].split("\"")[1]
                    if chrom not in dico_transcript_position.keys():
                        dico_transcript_position[chrom] = {transcript_name:[start,stop,direction]}
                    else:
                        dico_transcript_position[chrom][transcript_name] = [start,stop,direction]
    return dico_transcript_position


def create_final_big_file(dico):
    os.system ("rm -r old_prot_blast")
    os.system ("rm -r blast_denovo_to_target_genome")
    my_file = open("Table_output_step3.txt", "w")
    my_file.write("denovo,target_species,validated_hit,start,stop,Indels,frameshift,substitutions,premature_stop,pos_premature_stop,transcription_status,splicing_site,motifs" + "\n")
    #[validated_hit, ATG, stop, indels, frameshift, nb_subs, premature_stop, pos_premature_stop, transcription, splicingSite, motifs]
    for gene in dico.keys():
        for species in dico[gene]:
            list_elts = dico[gene][species]
            my_file.write(gene + "," + species)
            for elt in list_elts:
                my_file.write("," + str(elt))
            my_file.write("\n")
    my_file.close()
    os.system ("mv Table_output_step3.txt DESMAN_denovo_output")
