# Introduction
> [!IMPORTANT]
**DESMAN** is being benchmarked and corrected. Please wait for the publication and for the correct structuring of the code before using it.  

**DESMAN** (De novo Emergence, Shared Mutations, And Nucleotides) is a software that detect neORFs (precursors of _de novo_ genes), based on transcriptome data, and study their mutations within populations and/or species. 

# Usage

**DESMAN** (De novo Emergence, Shared Mutations, And Nucleotides), performs runs with 3 main steps: 
- (1) detect neORF candidates in transcriptomes,
- (2) validate the absence of homology to any known gene
- (3) search for syntenic homologous sequences in outgroup genomes (+ optionaly transcriptomes) and analyzing coding mutations between homologs.

DESMAN is available as a user-friendly graphical user interface, offering users a high degree of flexibility with various options. Regardless of its flexibility, **DESMAN** can operate under two distinct strategies.

**$${\color{orange}Strategy 1 : \space}$$** \
The user possesses: 1 transcriptome assembled with a reference genome and several optional outgroup genomes (with corresponding transcriptomes, if available). The transcriptome has been assembled by mapping stranded RNA-seq data to the reference genome, and the user aims to determine whether this transcriptome contains neORFs. Additionally, the user seeks to detect whether these neORF candidates can be found in the outgroup genomes and to what extent they are conserved. 

For example, let's consider a user who sequenced RNA from *Drosophila melanogaster* and assembled a transcriptome by mapping the RNA-seq data to the *D. melanogaster* reference genome. Using Strategy 1, the user will: 
* Extract all putative neORFs from the transcriptome.
* Select neORFs that show no homology to known proteins in Drosophila and, optionally, in outgroup species.
* Search for syntenic homologous sequences in outgroup genomes (either from *D. melanogaster* or from outgroup species) and analyze the mutations between the neORFs and their syntenic homologous sequences to study de novo emergence. 

**$${\color{orange}Strategy 2 : \space}$$** \
The user possesses several transcriptomes assembled with a single reference genome. The user has sequenced multiple RNA-seq datasets (for example, RNA-seq data from different organs or conditions of one species), and the transcriptomes were assembled by mapping the RNA-seq data to the same reference genome. With Strategy 2, the user aims to extract all candidate neORFs from each transcriptome. In the second step, the user seeks to identify which candidates are detected across multiple transcriptomes, and, for example, transcribed under different conditions. 

For example, let's say a user sequenced RNA-seq data from six different *D. melanogaster* samples, each extracted under different conditions, and assembled six different transcriptomes by mapping the RNA to the *D. melanogaster* reference genome. Using Strategy 2, the user will:: 
* Extract all putative neORFs from the transcriptome.
* Select neORFs that show no homology to known proteins in Drosophila and, optionally, in outgroup species.
* Determine which neORFs are expressed across multiple transcriptomes.

> [!NOTE]
A manual is available to access more precisely all steps made by DESMAN

# Flowchart

![Flowchart](flowchart.png)

# Install DESMAN

The is 2 main strategy to run DESMAN : 
* Install manually DESMAN dependencies on your computer before running DESMAN
* Use the container where everything is set-up
> [!WARNING]
> DESMAN runs only under linux or OS distributions

> [!IMPORTANT]
> We encourage users to work with the container (**Option 2**)

## Option 1. Install manually

DESMAN requires the installation of the following softwares
* BLAST (v2.12 or +) 
* DIAMOND (v2.0.14.152 or +) 
* GffRead (v0.12.7 or +) 


Moreover, DESMAN is developped in python (v3.0 or +). It therefore requires to have python installed, with the following packages:

* Biopython (vs 1.83) 
* customtkinter (5.2.2) 

### 1. Install [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (linux)

```console
user@comp:~/directory$ rpm -ivh ncbi-blast-2.2.18-1.x86\_64.rpm
```

### 2. Install [DIAMOND](https://github.com/bbuchfink/diamond/wiki) (linux)

```console
user@comp:~/directory$ wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
user@comp:~/directory$ tar xzf diamond-linux64.tar.gz
```


### 3. Install [GffRead](https://github.com/gpertea/gffread) (linux)

```console
user@comp:~/directory$ cd /some/build/dir
user@comp:~/directory$ git clone https://github.com/gpertea/
user@comp:~/directory$ gffread
user@comp:~/directory$ cd gffread
user@comp:~/directory$ make release
```

### 4. Install [biopython](https://biopython.org/wiki/Download) (linux)

```console
user@comp:~/directory$ pip install biopython
```

### 5. Install [customTkinter](https://pypi.org/project/customtkinter/0.3/) (linux)

```console
user@comp:~/directory$ sudo apt-get install python3-tk
user@comp:~/directory$ pip install customtkinter==0.3
```

## Option 2. Use the container


# Run DESMAN

DESMAN runs with python3. To run DESMAN, the user must call main_program.py, followed by the chosen Strategy. The strategies must be either **Strategy1** or **Strategy2**

```console
user@comp:~/directory$ python3 DESMAN.py Strategy1
```

# DESMAN options

## Mandatory parameters

## Optional parameters

| option      | choices | default     | description     |
| :---        |    :----:   |  :---: | :---  |
| Minimum TPM threshold      | float       | 0.5   |  The threshold for transcript expression. Any transcripts with expression levels lower than this value will not be considered.  |
| Filter transcripts overlapping with TEs   | boolean        | False      |This option filters the transcripts that overlap with TEs. This option only works if the genomes have been masked for TEs during annotation |
| Remove transcripts with unknown (".") orientation   | boolean        | False      |This option removes all transcripts whose orientation is unknown ("." in the GTF file) |
| Genomic position of candidate neORFs (based on transcript position)   | Intergenic  Intronic  Antisense  Genic     | Intergenic      |This option allows to select the genomic location of the transcripts. the options are cumulative for this setting |
| Remove genic splice variants   | boolean   | False      |DESMAN will filter the ORFs located on *de novo* transcripts whose one or more spliced variants overlap with a gene. |
| ORFs to keep in candidate transcripts   | All ORFs, ORF start first, ORF higest kozac, ORFs longest   | longest      |Select specific ORFs in transcripts |
| Set up a minimum size for 5'UTR and 3'UTR  |int, int   | 0,0      |With this option, you set a minimum size for the UTRs in the transcripts. By default, DESMAN does not impose any size filter for ORFs. Therefore, the selected ORF(s) in the transcript can start at the beginning of the transcript or have a stop codon at the end of the transcript |
| Set up a minimum size for 5'UTR and 3'UTR  |int, int   | 0,0      |With this option, you set a minimum size for the UTRs in the transcripts. By default, DESMAN does not impose any size filter for ORFs. Therefore, the selected ORF(s) in the transcript can start at the beginning of the transcript or have a stop codon at the end of the transcript |
