from Bio.SeqRecord import SeqRecord
import os
from module_assess_input import *
from module_blast_and_diamond import *
from module_select_orphans import *

cwd = os.getcwd()
Strategy = openFile("Strategy.setup")
dico_variables = withdraw_strategy_data(Strategy)

#print (dico_variables["link_database_outgroup_prot"])
#print (dico_variables["parameters_database_prot"])
#print (dico_variables["link_database_outgroup_nucl"])
#print (dico_variables["parameters_database_nucl"])



print ("DESMAN : Start Part 2")
list_name_to_remove_nucl = []
list_name_to_remove_prot = []
list_name_to_remove_total = []
if "link_database_outgroup_prot" in dico_variables.keys():
    print ("performing BLASTp")
    perform_blast_prot(dico_variables["link_database_outgroup_prot"],dico_variables["parameters_database_prot"])
    list_name_to_remove_prot = retrieve_name_hit_prot_blast()
    
if "link_database_outgroup_nucl" in dico_variables.keys():
    print ("performing BLASTn")
    perform_blast_nucl(dico_variables["link_database_outgroup_nucl"],dico_variables["parameters_database_nucl"])
    list_name_to_remove_nucl = retrieve_name_hit_nucl_blast(dico_variables["parameters_database_nucl"])
    


list_name_to_remove_total = merge_hit_lists(list_name_to_remove_prot, list_name_to_remove_nucl)

dico_size_denovo_nucl = get_file_size("DESMAN_denovo_output/denovo_nucl.fa")
dico_size_denovo_prot = get_file_size("DESMAN_denovo_output/denovo_protein.fa")

print ("Number of denovo candidate with a hit : " + str(len(list_name_to_remove_total)) + " out of " + str(len(dico_size_denovo_prot)))
reshufe_files_in_denovo(list_name_to_remove_total)
# seems i have to re-do the coverage as the one give by blast is wrong

print ("DESMAN Part 2 : DONE!")
