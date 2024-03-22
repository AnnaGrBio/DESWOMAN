import os
#os.chdir("/home/agrandch/Desktop/postdoc/Projet_principal/DESMAN/Get_and_sort_transcripts_ORFs")
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import MultipleSeqAlignment
from Bio.pairwise2 import format_alignment




def openFile(NameFile):
    
    """
    open file in lecture mode 
    """
    
    F=open(NameFile, "r")
    L=F.readlines()
    return L  

#the paramaters and functions writen here will have to be downloaded directly from the previous code. 
# they are re-writen here for now cause i need them also in this part of the code

link_to_close_genomes = "/home/agrandch/Desktop/postdoc/Projet_principal/DESMAN/Get_and_sort_transcripts_ORFs/outgroup_genomes"

def extract_genomes(link_to_close_genomes):
    # this function has to be absolutely changed as now it just returns the name of one genome outgroup YE but should return all genomes from the database
    return ["YE_finalGenome.masked.fa"]
list_genomes = extract_genomes(link_to_close_genomes)

parameters_blast_genome = {"type" : "nucl", "e_value" : "0.05", "coverage" : 80}
############################################################


#name_output = ID_genome + "_all_bast_output.txt"
#os.system ("mv " + name_output + " blast_denovo_to_target_genome")

### part 1 make dico of correct blast outputs (all of them if no coverage conditions, some of them if coverage conditions)









### synteny part
    









        


            
###
            











# make alignements




#################################




        



main_alignment_function(dico_denovo_seq_intron, dico_NcHit)



## add a chromosome constrint for synteny


