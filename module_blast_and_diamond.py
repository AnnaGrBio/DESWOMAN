import os


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  
    

def perform_blast_nucl(path_to_db, parameters_db):
    # Constructing the command to create a BLAST database
    command1 = "makeblastdb -in " + path_to_db + " -dbtype nucl"
    
    # Executing the command to create the BLAST database and redirecting output to a file
    os.system(command1 + ">> blast_n_makebd.txt")
    
    # Setting the name for the output file of the BLAST search
    name_output = "denovo_blast_output_nucl.txt"
    
    # Constructing the BLAST command based on parameters provided
    if parameters_db["strand"] == None:  # Checking if the strand parameter is provided
        command2 = "blastn -evalue " + parameters_db["e_value"] + " -query " + "DESMAN_denovo_output/denovo_nucl.fa -db " + path_to_db + " -out " + name_output + " -outfmt \"7 qacc sacc evalue Identities qstart qend sstart send qcovs\""
    else:
        command2 = "blastn -strand " + parameters_db["strand"] + " -evalue " + parameters_db["e_value"] + " -query " + "DESMAN_denovo_output/denovo_nucl.fa -db " + path_to_db + " -out " + name_output + " -outfmt \"7 qacc sacc evalue Identities qstart qend sstart send qcovs\""
    
    # Executing the BLAST search command and redirecting output to a file
    os.system(command2 + ">> blast_n_output.txt")
    
    # Moving the output file to a specific directory
    os.system("mv " + name_output + " DESMAN_denovo_output")


def perform_blast_prot_transcripts():
  command1 = "./diamond makedb --in Intermediate_prot_BLAST/target_prot.fa -d target_prot"
  os.system (command1)
  print ("lala")
  name_output = "denovo_blast_output_prot.txt"
  command2 = "./diamond blastp -d target_prot -q Intermediate_prot_BLAST/query_prot.fa --more-sensitive -o diamond_transc_prot.out"
  os.system (command2)
  os.system ("mv diamond_transc_prot.out Intermediate_prot_BLAST")
  os.system ("rm *.dmnd")


def perform_blast_prot(path_to_db, parameters_db):
    # Extracting the name of the database from the path
    list_path = path_to_db.split("/")
    name_database = list_path[len(list_path) - 1].split(".")[0]
    
    # Constructing the command to create a Diamond database
    command1 = "./diamond makedb --in " + path_to_db + " -d " + name_database
    
    # Executing the command to create the Diamond database
    os.system(command1)
    
    # Setting the name for the output file of the Diamond search
    name_output = "denovo_blast_output_prot.txt"
    
    # Constructing the Diamond command based on parameters provided
    command2 = "./diamond blastp -d " + name_database + " -q " + "DESMAN_denovo_output/denovo_protein.fa " + parameters_db["mode"] + " -o " + name_output
    
    # Executing the Diamond search command
    os.system(command2)
    
    # Moving the output file to a specific directory
    os.system("mv " + name_output + " DESMAN_denovo_output")
    
    # Removing temporary Diamond database files
    os.system("rm *.dmnd")

    

def perform_blast_to_genome(ID_genome, link_to_my_target_genome):
    # Defining parameters for BLAST to genome
    parameters_blast_genome = {"type": "nucl", "e_value": "0.01", "coverage": None}
    
    # Setting the name for the output file of the BLAST search
    name_output = ID_genome + "_all_bast_output.txt"
    
    # Creating a BLAST database from the target genome
    os.system("makeblastdb -in " + link_to_my_target_genome + " -dbtype nucl")
    
    # Constructing the BLAST command based on parameters provided
    command = "blastn -evalue " + parameters_blast_genome["e_value"] + " -query DESMAN_denovo_output/denovo_unspliced_lowered_introns.fa" + " -db " + link_to_my_target_genome + " -out " + name_output + " -outfmt \"7 qacc sacc evalue Identities qstart qend sstart send qcovs\""
    
    # Executing the BLAST search command
    os.system(command)
    
    # Moving the output file to a specific directory
    os.system("mv " + name_output + " blast_denovo_to_target_genome")

    # Removing temporary BLAST database files
    list_path_to_folder = link_to_my_target_genome.split("/")
    path_to_folder = ""
    for i in list_path_to_folder[0:len(list_path_to_folder) - 1]:
        path_to_folder += i + "/"
    os.system("rm " + path_to_folder + "*.ndb")
    os.system("rm " + path_to_folder + "*.nin")
    os.system("rm " + path_to_folder + "*.not")
    os.system("rm " + path_to_folder + "*.nsq")
    os.system("rm " + path_to_folder + "*.ntf")
    os.system("rm " + path_to_folder + "*.nto")
    os.system("rm " + path_to_folder + "*.nhr")
    
    
def get_transcription_hits(dico_target_species_lines, pop_species_name):
    dico_transcript_hits = {}  # Dictionary to store transcript hits
    
    # Checking if transcriptome GTF file exists in the dictionary for the given species
    if "transcriptome_gtf" in dico_target_species_lines[pop_species_name].keys():
        
        # Obtaining the path to the transcriptome FASTA file
        link_transcripts_fasta = dico_target_species_lines[pop_species_name]["transcriptome_fasta"]
        
        # Creating a BLAST database from the transcriptome FASTA file
        command1 = "makeblastdb -in " + link_transcripts_fasta + " -dbtype nucl"
        os.system(command1)
        
        # Performing BLAST search against the transcriptome
        command2 = "blastn -evalue 0.01 -query DESMAN_denovo_output/denovo_nucl.fa -db " + link_transcripts_fasta + " -out transcript_blast_output.txt -outfmt \"7 qacc sacc evalue Identities qstart qend sstart send qcovs\""
        os.system(command2)
        
        # Removing temporary BLAST database files
        list_path_to_folder = link_transcripts_fasta.split("/")
        path_to_folder = ""
        for i in list_path_to_folder[0:len(list_path_to_folder)-1]:
            path_to_folder += i + "/"
        os.system("rm " + path_to_folder + "*.ndb")
        os.system("rm " + path_to_folder + "*.nin")
        os.system("rm " + path_to_folder + "*.not")
        os.system("rm " + path_to_folder + "*.nsq")
        os.system("rm " + path_to_folder + "*.ntf")
        os.system("rm " + path_to_folder + "*.nto")
        os.system("rm " + path_to_folder + "*.nhr")
        
        # Processing BLAST output to extract hits
        my_blast_hits = openFile("transcript_blast_output.txt")
        for line in my_blast_hits:
            if line[0] != "#":
                elts_line = line.split()
                query = elts_line[0]
                target = elts_line[1]
                cov = int(elts_line[7])
                value_cov = ""
                
                # Determining coverage type
                if cov == 100:
                    value_cov = "complete"
                else:
                    value_cov = "partial"
                
                # Storing hits in the dictionary
                if query not in dico_transcript_hits.keys():
                    dico_transcript_hits[query] = {target: value_cov}
                else:
                    if target not in dico_transcript_hits[query].keys():
                        dico_transcript_hits[query][target] = value_cov
        
        # Removing temporary BLAST output file
        os.system("rm transcript_blast_output.txt")
    
    return dico_transcript_hits


