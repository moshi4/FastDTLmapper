# FastDTLmapper: Fast genome-wide DTL event mapper  

![Python3](https://img.shields.io/badge/Language-Python_3-steelblue)
![OS](https://img.shields.io/badge/OS-Linux-steelblue)
![License](https://img.shields.io/badge/License-GPL3.0-steelblue)

## Table of contents

- [Overview](#overview)
- [Install](#install)
  - [Dependencies](#dependencies)
- [Analysis Pipeline](#analysis-pipeline)
- [Command Usage](#command-usage)
  - [Run command](#run-command)
  - [Options](#options)
  - [Example](#example)
- [Output detail](#output-detail)

## Overview

**FastDTLmapper** can fastly map genome-wide DTL(Duplication-Transfer-Loss) event  
from "Rooted species tree" and "Target species fasta files"

## Install

FastDTLmapper is implemented in **Python3(>=3.7)** and runs on **Linux**.  

Install by using pip command below.  

    pip install git+https://github.com/moshi4/FastDTLmapper.git

### Dependencies

All of the following dependencies are packaged in the **bin** directory.  

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
6. 'Species tree' & 'each OG rooted gene tree' DTL reconciliation using Mowgli
7. Aggregate and map DTL reconciliation result

## Command Usage

### Run command

    FastDTLmapper.py -i [fasta|genbank directory] -t [newick rooted species file] -o [output directory]

### Options

    -h, --help           show this help message and exit
    -i , --indir         Input Fasta(*.fa|*.faa|*.fasta), Genbank(*.gb|*.gbk|*.genbank) directory
    -t , --treefile      Input rooted species time(ultrametric) tree file (Newick format)
    -o , --outdir        Output directory
    -p , --process_num   Number of processor (Default: max processor - 1)
    --dup_cost           Duplication event cost (Default: 2)
    --los_cost           Loss event cost (Default: 1)
    --trn_cost           Transfer event cost (Default: 3)
    --rseed              Number of random seed (Default: 0)

### Example

    FastDTLmapper.py -i [fasta|genbank directory] -t [newick rooted species file] -o [output directory]

## Output detail

| Directory               | Contents                                                    |
| ----------------------- | ----------------------------------------------------------- |
| 00_user_data            | Formatted user fasta and tree files                         |
| 01_orthofinder          | OrthoFinder raw output result                               |
| 02_dtl_reconciliation   | Each OG(Ortholog Group) DTL reconciliation results<br>OGXXXXXXX/ --- One OG output directory example<br>├ \*.fa --- OG fasta file<br>├ \*_aln.fa --- OG alignment fasta file<br>├ \*_aln_trim.fa --- Trimmed OG alignment fasta file<br>├ \*_unrooted_tree.nwk --- OG unrooted gene tree file<br>├ \*_rooted_tree.nwk --- OG rooted gene tree file<br>├ *_gain_loss_map.nwk --- Gain-loss map tree file<br>├ *_dtl_map.nwk --- DTL map tree file<br>└ dtl_reconciliation/ --- DTL reconciliation output directory|
| 03_aggregate_map_dtl    | Genome-wide DTL reconciliation aggregated and mapped result |
