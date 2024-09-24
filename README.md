# DME: The Drug Metabolising Enzyme software

# Introduction

The microbiome can break down or modify xenobiotics (i.e., externally administered drugs), regulate 
host gene expression, and modulate xenobiotic absorption. The drug metabolism is not limited to orally 
administered drugs, the microbiome also converts drug metabolites destined for excretion via the gut, 
including drug conjugates from the liver

This application is to predict drug-metabolizing enzymes from protein or nucleotide data based on 
homology models. The application uses reference data from the 
[Human Microbiome Drug Metabolism (HMDM) Database](https://hmdm.mcmaster.ca).

# Overview of DME


![dme overview](https://github.com/raphenya/dme/blob/main/docs/images/DME.png?raw=true)


If DNA sequences are submitted, DME first predicts complete open reading frames (ORFs) using 
[Prodigal](https://github.com/hyattpd/Prodigal) (ignoring those less than 30 bp) and analyzes 
the predicted protein sequences. This includes a secondary correction by DME if Prodigal undercalls 
the correct start codon to ensure complete genes are predicted. However, if Prodigal fails to predict 
an ORF for a gene, this will produce a false negative result.

If protein sequences are submitted, DME skips ORF prediction and uses the protein sequences directly.

The DME currently supports HMDM's [protein homolog models](https://hmdm.mcmaster.ca/cvterms/54) (use of BLASTP or [DIAMOND](https://ab.inf.uni-tuebingen.de/software/diamond) bitscore cut-offs to detect functional homologs of enzymes).


|    Example                                               | Enzyme                                            |
| ---------------------------------------------------------| --------------------------------------------------|
|    Protein Homolog Model                                 | [TDC](https://hmdm.mcmaster.ca/cvterms/24)        |


The DME analyzes genome or proteome sequences under three paradigms: **Perfect**, **Strict**, and **Loose** (a.k.a. Discovery).

# Installation

The tool requires Python >= 3.11 and conda >= 4.12.0. The latest release can be installed directly from pip or this repository.

```
pip install dme
```

or 

```
python3 -m pip install /path/to/dme-1.0.0.tar.gz
```


# Dependencies

- Python 3.6
- NCBI BLAST 2.9.0
- zlib
- Prodigal 2.6.3
- DIAMOND 0.8.36
- Biopython 1.78
- filetype 1.0.0+
- pytest 3.0.0+ 
- pandas 0.15.0+
- Matplotlib 2.1.2+
- seaborn 0.8.1+ 
- pyfaidx 0.5.4.1+
- pyahocorasick 1.1.7+
- OligoArrayAux 3.8
- samtools 1.9
- bamtools 2.5.1
- bedtools 2.27.1
- Jellyfish 2.2.10
- Bowtie2 2.3.4.3
- BWA 0.7.17 (r1188)
- KMA 1.3.4


# DME Usage

```
usage: dme <command> [<args>]
            commands are:
               ---------------------------------------------------------------------------------------
               Database
               ---------------------------------------------------------------------------------------
               auto_load Automatically loads HMDM database, annotations and k-mer database
               load      Loads HMDM database, annotations and k-mer database
               clean     Removes BLAST databases and temporary files
               database  Information on installed hmdm database
               galaxy    Galaxy project wrapper

               ---------------------------------------------------------------------------------------
               BLAST
               ---------------------------------------------------------------------------------------
               blast     Runs NCBI BLAST algorithm on HMDM data

               ---------------------------------------------------------------------------------------
               Genomic
               ---------------------------------------------------------------------------------------

               main     Runs dme application
               tab      Creates a Tab-delimited from dme results
               parser   Creates categorical JSON files DME wheel visualization
               heatmap  Heatmap for multiple analysis

               ---------------------------------------------------------------------------------------
               Annotations
               ---------------------------------------------------------------------------------------
               hmdm_annotation       Create fasta files with annotations from hmdm.json
               wildhmdm_annotation   Create fasta files with annotations from variants
               baits_annotation      Create fasta files with annotations from baits (experimental)
               remove_duplicates     Removes duplicate sequences (experimental)

               

Drug Metabolising Enzyme (DME) - 1.0.0

positional arguments:
  {main,tab,parser,load,auto_load,clean,galaxy,database,bwt,tm,hmdm_annotation,wildhmdm_annotation,baits_annotation,remove_duplicates,heatmap,kmer_build,kmer_query}
                        Subcommand to run

options:
  -h, --help            show this help message and exit

Use the Drug Metabolising Enzyme (DME) to predict drug-metabolizing enzymes from protein or nucleotide data based on homology models. Check https://hmdm.mcmaster.ca/download for software and
data updates. Receive email notification of monthly HMDM updates via the HMDM Mailing List (https://mailman.mcmaster.ca/mailman/listinfo/hmdm-l)
```