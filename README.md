# Introduction
> [!IMPORTANT]
**DESMAN** is being benchmarked and corrected. Please wait for the publication and for the correct structuring of the code before using it.  

**DESMAN** is a software that detect neORFs (precursors of _de novo_ genes), based on transcriptome data, and study their mutations within populations and/or species. 



**DESMAN** DESMAN (De novo Emergence, Shared Mutations, And Nucleotides), performs runs with 3 main steps: 
- (1) detect neORF candidates in transcriptomes,
- (2) validate the absence of homology to any known gene
- (3) search for syntenic homologous sequences in outgroup genomes (+ optionaly transcriptomes) and analyzing coding mutations between homologs. DESMAN is available as a user-friendly graphical user interface, offering users a high degree of flexibility with various options. 
Regardless of its flexibility, **DESMAN** can operate under two distinct strategies.

- **Strategy 1** : In Strategy 1, the user possesses: 1 transcriptome assembled with a reference genome and several optional outgroup genomes (with corresponding transcriptomes, if available). The transcriptome has been assembled by mapping stranded RNA-seq data to the reference genome, and the user aims to determine whether this transcriptome contains neORFs. Additionally, the user seeks to detect whether these neORF candidates can be found in the outgroup genomes and to what extent they are conserved. 
For example, let's consider a user who sequenced RNA from *Drosophila melanogaster* and assembled a transcriptome by mapping the RNA-seq data to the *D. melanogaster* reference genome. Using Strategy 1, the user will: 
- Extract all putative neORFs from the transcriptome.
- Select neORFs that show no homology to known proteins in Drosophila and, optionally, in outgroup species.
- Search for syntenic homologous sequences in outgroup genomes (either from *D. melanogaster* or from outgroup species) and analyze the mutations between the neORFs and their syntenic homologous sequences to study de novo emergence. 

- **Strategy 2** : In Strategy 2, the user possesses several transcriptomes assembled with a single reference genome. The user has sequenced multiple RNA-seq datasets (for example, RNA-seq data from different organs or conditions of one species), and the transcriptomes were assembled by mapping the RNA-seq data to the same reference genome. With Strategy 2, the user aims to extract all candidate neORFs from each transcriptome. In the second step, the user seeks to identify which candidates are detected across multiple transcriptomes, and, for example, transcribed under different conditions. 
For example, let's say a user sequenced RNA-seq data from six different *D. melanogaster* samples, each extracted under different conditions, and assembled six different transcriptomes by mapping the RNA to the *D. melanogaster* reference genome. Using Strategy 2, the user will:: 
- Extract all putative neORFs from the transcriptome.
- Select neORFs that show no homology to known proteins in Drosophila and, optionally, in outgroup species.
- Determine which neORFs are expressed across multiple transcriptomes.

![Flowchart](flowchart.png)
