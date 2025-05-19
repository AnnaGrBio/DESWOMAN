import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from module_colors import *


__author__ = "Anna Grandchamp"
__contributor__=""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def search_surrounding_gene(start : int, stop : int, chrom : str, dico_gff_focal_sorted : dict, window : int) -> list:
    """
    Search for genes surrounding a candidate gene within a specified window on the same chromosome.
    Returns a list of surrounding genes based on the given window size.

    The function assumes that the provided `dico_gff_focal_sorted` dictionary contains 
    a sorted list of genes for each chromosome (sorted by gene start position). The window size 
    determines how many surrounding genes are considered to the left and right of the candidate gene.

    Parameters:
    -----------
    start : int
        The start position of the candidate gene.
    
    stop : int
        The stop position of the candidate gene (not used directly in this function, but included for completeness).
    
    chrom : str
        The chromosome (or scaffold) where the candidate gene is located.
    
    dico_gff_focal_sorted : dict
        A dictionary where keys are chromosome names and values are lists of genes, each represented as a 
        list of gene information (with the start position as the first element). The genes within each 
        list are sorted by their start position.
    
    window : int
        The size of the window around the candidate gene. This defines how many genes before and after 
        the candidate gene will be considered as "surrounding" genes.

    Returns:
    --------
    list_surrounding : list
        A list of surrounding genes, including those within the specified window around the candidate gene. 
        The genes are returned as they appear in the `dico_gff_focal_sorted` dictionary.

    Example:
    --------
    start = 500
    stop = 1000
    chrom = "chromosome_1"
    dico_gff_focal_sorted = {"chromosome_1": [[100, 200, "geneA"], [300, 400, "geneB"], [500, 600, "candidate_gene"], [700, 800, "geneC"], [900, 1000, "geneD"]]}
    window = 1
    surrounding_genes = search_surrounding_gene(start, stop, chrom, dico_gff_focal_sorted, window)
    # surrounding_genes will contain "geneB", "candidate_gene", and "geneC" (1 gene before and 1 gene after the candidate gene).
    """
    list_surrounding = []
    if chrom in dico_gff_focal_sorted:
        #the matrix is a list of all genes from a given chromosome/scaffold in an increasing order of position
        matrix = dico_gff_focal_sorted[chrom]
        # the classe variable indicates the candidate was positioned between its neigbooring genes
        if len(matrix) >= window*2: # means if the chromosome contains more genes that 2 times the window
            matrix_start = []
            dico_associated = {}
            # matrix_start implements the start of all genes in the chrom
            for sublist_gene_info in matrix:
                matrix_start.append(int(sublist_gene_info[0])) # the matrix implements the start position of the gene
                dico_associated [int(sublist_gene_info[0])] = sublist_gene_info
            # next the matrix of all gene start position will implement the start position of the candidate denovo gene/ hit
            matrix_start.append(int(start))
            # we then sort the matrix so that our candidate is positioned between the genes in an increase start position
            matrix_start = list(set(matrix_start))
            matrix_start.sort()
            # we then determine at which position was positioned the candidate denovo/hit
            index_of_denovo = matrix_start.index(start)
            # we then get the index of the gene at the last left extramity of the matrix. 
            # ex : if the denovo is at position 10, and the window is 3, then the gene in the window most on the left is at index 7
            start_before = index_of_denovo - window
            # Exception if the denovo is too on the left in the matrix we retrieve 0 or less genes on the left
            if start_before < 0:
                start_before = 0
            for i in range(start_before,index_of_denovo):
                list_surrounding.append(dico_associated[matrix_start[i]])
            # now we repeat the same thing but on the right
            start_after = index_of_denovo + window
            if start_after >= len(matrix_start) - 1:
                start_after = len(matrix_start) - 1
            for i in range(index_of_denovo + 1, start_after + 1):
                list_surrounding.append(dico_associated[matrix_start[i]])
        else:
            list_surrounding = matrix
    else:
        list_surrounding = []
    return list_surrounding


def search_overlapping_gene(start : int, stop : int, chrom : str, dico_gff_focal_sorted : dict) -> list:
    """
    Search for genes that overlap with a given candidate gene or homologous hit on the same chromosome.
    
    The function checks if the start and stop positions of the candidate gene (or homologous hit) 
    overlap with any genes present in the provided chromosome's gene list. Overlapping genes are 
    identified by comparing the position intervals.

    Parameters:
    -----------
    start : int
        The start position of the candidate gene or homologous hit.
    
    stop : int
        The stop position of the candidate gene or homologous hit.
    
    chrom : str
        The chromosome (or scaffold) on which the candidate gene is located.
    
    dico_gff_focal_sorted : dict
        A dictionary where keys are chromosome names (strings) and values are lists of genes 
        on that chromosome. Each gene is represented as a list, where the first element is 
        the start position and the second element is the stop position of the gene.

    Returns:
    --------
    list_overlap : list
        A list of genes that overlap with the provided candidate gene or homologous hit.
        The genes are returned as they appear in the `dico_gff_focal_sorted` dictionary.

    Example:
    --------
    start = 500
    stop = 1000
    chrom = "chromosome_1"
    dico_gff_focal_sorted = {"chromosome_1": [[100, 200, "geneA"], [300, 400, "geneB"], [500, 600, "candidate_gene"], [700, 800, "geneC"], [900, 1000, "geneD"]]}
    overlapping_genes = search_overlapping_gene(start, stop, chrom, dico_gff_focal_sorted)
    # overlapping_genes will contain "geneB", "candidate_gene", and "geneC", as they overlap with the region 500-1000.
    """
    list_overlap = []
    if chrom in dico_gff_focal_sorted:
        matrix = dico_gff_focal_sorted[chrom]
        #if len(matrix) == 1:
            #list_overlap.append(matrix[0])
        #else:
        for sublist_gene_info in matrix:
            start_target = int(sublist_gene_info[0])
            stop_target = int(sublist_gene_info[1])
            if start <= start_target:
                if stop >= start_target:
                    list_overlap.append(sublist_gene_info)
            elif start >= start_target and stop <= stop_target:
                list_overlap.append(sublist_gene_info)
            elif start <= stop_target and stop >= stop_target:
                list_overlap.append(sublist_gene_info)
    else:
        list_overlap = []
    return list_overlap


def determine_neigboor_function(denovo_genomic_infos : list, dico_gff_focal_sorted : dict, window : int) -> list:
    """
    Determines the neighboring genes of a given de novo gene or its homologous hit based on its genomic type 
    (intergenic, genic, antisense, or intronic). It finds surrounding genes or overlapping genes depending on 
    whether the gene is intergenic or within another gene.

    Parameters:
    -----------
    denovo_genomic_infos : list
        A list containing genomic information about the de novo gene or homologous hit. 
        The list must contain the following elements:
        - [0] : Genic type of the de novo gene (e.g., "intergenic_transcripts", "genic", "antisense", "intronic").
        - [1] : Chromosome name (string).
        - [2] : Start position of the gene (integer).
        - [3] : Stop position of the gene (integer).
    
    dico_gff_focal_sorted : dict
        A dictionary where keys are chromosome names (strings) and values are lists of genes in increasing order 
        of position on the chromosome. Each gene is represented as a list, where the first element is the 
        start position, the second element is the stop position, and the third element is the gene name.
    
    window : int
        The size of the window around an intergenic gene to search for neighboring genes. 
        This is used when the gene is classified as intergenic.
    
    Returns:
    --------
    list_neighboor : list
        A list of neighboring genes (either overlapping or surrounding) depending on the genomic type of the de novo gene.
    """
    start = int(denovo_genomic_infos[2])
    stop = int(denovo_genomic_infos[3])
    chrom = denovo_genomic_infos[1]
    if denovo_genomic_infos[0] == "intergenic_transcripts":
        list_neighboor = search_surrounding_gene(start,stop,chrom,dico_gff_focal_sorted,window)
    else:
        list_neighboor = search_overlapping_gene(start,stop,chrom,dico_gff_focal_sorted)
    return list_neighboor


def get_synteny_score(list_neighboor_denovo : list, list_neighboor_hit : list, dico_blastA : dict, dico_blastB : dict) -> int:
    """
    Evaluates if at least one of the genes neighboring the de novo gene and the hit are syntenic.
    Synteny is determined by checking if there is at least one reciprocal blast hit between the genes 
    in the neighborhood of the de novo gene and the homologous hit.

    Parameters:
    -----------
    list_neighboor_denovo : list
        A list of neighboring genes around the de novo gene. Each gene is represented as a list, where 
        the third element is the gene name.

    list_neighboor_hit : list
        A list of neighboring genes around the homologous hit. Each gene is represented as a list, where 
        the third element is the gene name.

    dico_blastA : dict
        A dictionary where keys are gene names (strings) and values are lists of homologous genes 
        from the de novo gene's neighborhood, as determined by blast hits.

    dico_blastB : dict
        A dictionary where keys are gene names (strings) and values are lists of homologous genes 
        from the homologous hit's neighborhood, as determined by blast hits.

    Returns:
    --------
    nb_reciprocal_hit : int
        The number of reciprocal blast hits found between the neighboring genes of the de novo gene and 
        the homologous hit. If at least two reciprocal hits are found, the homologous hit is considered syntenic.
    """
    new_list_neighboor_denovo = []
    #Retrieve the name of each genes from the neigboring genes of denovo
    for list_gene in list_neighboor_denovo:
        new_list_neighboor_denovo.append(list_gene[2])
    new_list_neighboor_hit = []
    #Retrieve the name of each genes from the neigboring genes of homologous hit
    for list_gene in list_neighboor_hit:
        new_list_neighboor_hit.append(list_gene[2])
    nb_reciprocal_hit = 0
    end_function = False
    # For each genes in the denovo neighboor, check if they have reciprocal blast to neigboor of homologous hit
    for gene_neighboor_denovo in new_list_neighboor_denovo:
        if gene_neighboor_denovo in dico_blastA:
            list_homologs_gene_neigboor_denovo = dico_blastA[gene_neighboor_denovo]
            for homologs_gene_neigboor_denovo in list_homologs_gene_neigboor_denovo:
                # keep advancing analysis only if the genes homologs to the genes surounding denovo are among genes surrounding target homolog hit
                if homologs_gene_neigboor_denovo in new_list_neighboor_hit:
                    if homologs_gene_neigboor_denovo in dico_blastB and gene_neighboor_denovo in dico_blastB[homologs_gene_neigboor_denovo]:
                            nb_reciprocal_hit += 1
                            # If there re two recirpocal hots, we can stop the function, the homologous hit is syntenic
                            if nb_reciprocal_hit == 2:
                                end_function = True
                            break
            if end_function == True:
                break
    return nb_reciprocal_hit


def search_all_denovo_hits(dico_denovo_informations : dict, dict_my_blast_selected_outputs : dict, dico_gff_focal_sorted : dict, dico_gff_outgroup_sorted : dict, window : int, dico_blast1 : dict, dico_blast2 : dict) -> dict:
    """
    Detects syntenic homologous hits for de novo genes by evaluating the neighborhoods of both the de novo gene and its potential homologs.
    It uses reciprocal blast hits to determine synteny between de novo genes and their homologs in the target genome.
    
    Parameters:
    -----------
    dico_denovo_informations : dict
        A dictionary where keys are de novo gene names and values are genomic information for each gene, 
        such as its genomic location and the chromosome it belongs to.
        
    dict_my_blast_selected_outputs : dict
        A dictionary where keys are de novo gene names and values are lists of the homologous hits identified by BLAST.

    dico_gff_focal_sorted : dict
        A dictionary where keys are chromosome names and values are lists of genes in the target genome sorted by their start positions. 
        This is used to find neighboring genes for each de novo gene.

    dico_gff_outgroup_sorted : dict
        A dictionary where keys are chromosome names and values are lists of genes in the outgroup genome sorted by their start positions.
        This is used to find neighboring genes for each homologous hit.
        
    window : int
        The size of the genomic window used to define the neighborhood of a gene. It determines how many genes around a given gene are considered 
        to be in its neighborhood.

    dico_blast1 : dict
        A dictionary of blast results for the first set of genes. This is used to check for reciprocal hits from the de novo gene to its homologous neighbors.
        
    dico_blast2 : dict
        A dictionary of blast results for the second set of genes. This is used to check for reciprocal hits from the homologous neighbors back to the de novo gene.

    Returns:
    --------
    dico_denovo_best_hit : dict
        A dictionary where keys are de novo gene names and values are the best syntenic homologous hits found for each de novo gene.
    """
    # the dictionnary will store the best syntenic hit for each denovo genes.
    dico_denovo_best_hit = {}
    nb_syntenic = 0
    nb_total_denovo = 0
    for denovo in dico_denovo_informations:
        # counts total number of denovo examinated
        nb_total_denovo += 1
        denovo_genomic_infos = dico_denovo_informations[denovo]
        # searches syntenic hit only if the denovo has blast hits to the target genome.
        if denovo in dict_my_blast_selected_outputs:
            list_denovo_hits = dict_my_blast_selected_outputs[denovo]
            # extract annotated genes neighboor to the denovo. The number of max neighboor extracted depends on the window given by the user.
            list_neighboor_denovo = determine_neigboor_function(denovo_genomic_infos,dico_gff_focal_sorted,window)
            best_hit = []
            best_score = 0
            for hit in list_denovo_hits:
                list_elts_hit = hit.split("-")
                #name, chrom, start, stop
                hit_genomic_infos = [denovo_genomic_infos[0],list_elts_hit[2],list_elts_hit[0],list_elts_hit[1]]
                # extract annotated genes neighboor to the hit. The number of max neighboor extracted depends on the window given by the user.
                list_neighboor_hit = determine_neigboor_function(hit_genomic_infos,dico_gff_outgroup_sorted,window)
                nb_reciprocal_hit = get_synteny_score(list_neighboor_denovo,list_neighboor_hit, dico_blast1, dico_blast2)
                if nb_reciprocal_hit > 0:
                    if denovo_genomic_infos[0] == "intergenic_transcripts":
                        if nb_reciprocal_hit > 1:
                            dico_denovo_best_hit[denovo] = hit
                            nb_syntenic += 1
                            break
                    else:
                        dico_denovo_best_hit[denovo] = hit
                        nb_syntenic += 1
                        break  
    avg = nb_syntenic/nb_total_denovo
    print (BRIGHT_BLUE + "Syntenic hits extracted" + RESET)
    print (BRIGHT_BLUE + "Number of neORF candidate : " + RESET + str(len(dico_denovo_informations)))
    print (BRIGHT_BLUE + "Number of neORF candidate with syntenic hits : " + RESET + str(len(dico_denovo_best_hit)))
    print (BRIGHT_BLUE + "Percentage of candidate neORF with syntenic hits : " + RESET + str(avg))
    return dico_denovo_best_hit
    
