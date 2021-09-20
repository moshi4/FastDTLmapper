# fastDTLmapper: Fast Genome-wide DTL event mapper  

![Python3](https://img.shields.io/badge/Language-Python_3-steelblue)
![OS](https://img.shields.io/badge/OS-Linux-steelblue)
![License](https://img.shields.io/badge/License-GPL3.0-steelblue)

## Overview

**fastDTLmapper** can fastly map genome-wide DTL(Duplication-Transfer-Loss) event  
from "Rooted species tree" and "Target species fasta files"

## Install

fastDTLmapper is written in **Python3(>=3.7)** and can be run on **Linux**.  

Install by using pip command below.  

    pip install git+https://github.com/moshi4/fastDTLmapper.git

### Dependencies

All dependencies packaged in **bin** directory.

- [OrthoFinder](https://github.com/davidemms/OrthoFinder)  
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
- [trimal](http://trimal.cgenomics.org/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [OptRoot](https://compbio.engr.uconn.edu/software/RANGER-DTL/)
- [Mowgli](http://www.atgc-montpellier.fr/Mowgli/)
- [parallel](https://www.gnu.org/software/parallel/)

## Analysis Pipeline

1. Grouping ortholog sequences using OrthoFinder
2. Align each OG(Ortholog Group) sequences using mafft
3. Trim each OG alignment using trimal
4. Reconstruct each OG gene tree using FastTree
5. Rooting each OG gene tree using OptRoot
6. 'Species tree' & 'each OG gene tree' DTL reconciliation using Mowgli
7. Aggregate and map DTL reconciliation result

## Command Usage

### Run command

    fastDTLmapper.py -i [fasta|genbank directory] -t [newick rooted species file] -o [output directory]

### Options

    -h, --help           show this help message and exit
    -i , --indir         Input Fasta(*.fa|*.faa|*.fasta), Genbank(*.gb|*.gbk|*.genbank) directory
    -t , --tree_file     Input rooted species time(ultrametric) tree file (Newick format)
    -o , --outdir        Output directory
    -p , --process_num   Number of processor
    --dup_cost           Duplication event cost (Default: 2)
    --los_cost           Loss event cost (Default: 1)
    --trn_cost           Transfer event cost (Default: 3)
    --rseed              Number of random seed (Default: 0)

### Example

    fastDTLmapper.py -i [fasta|genbank directory] -t [newick rooted species file] -o [output directory]

## Output detail

| Directory name          | Contents                                                |
| ----------------------- | ------------------------------------------------------- |
| 00_user_data            | Formatted user fasta and tree file                      |
| 01_orthofinder          | OrthoFinder raw output result                           |
| 02_dtl_reconciliation   | align, trim, gene tree, rooted gene tree, dtl result    |
| 03_aggregate_dtl_result | Aggregate and map dtl result                            |
