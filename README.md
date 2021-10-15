# FastDTLmapper: Fast genome-wide DTL event mapper  

![Python3](https://img.shields.io/badge/Language-Python_3.7_|_3.8_|_3.9-steelblue)
![OS](https://img.shields.io/badge/OS-Linux-steelblue)
![License](https://img.shields.io/badge/License-GPL3.0-steelblue)  

## Table of contents

- [Overview](#overview)
- [Install](#install)
- [Method Summary](#method-summary)
- [Command Usage](#command-usage)
- [Output Summary](#output-summary)

## Overview

**FastDTLmapper** can fastly map genome-wide DTL(Duplication-Transfer-Loss) event  
from "Rooted species tree" and "Target species fasta files"

## Install

FastDTLmapper is implemented in **Python3(>=3.7)** and runs on **Linux**.  

~~Install with pip:~~

~~pip install fastdtlmapper~~

### Dependencies

All of the following dependencies are packaged in **src/fastdtlmapper/bin** directory.  

- [OrthoFinder](https://github.com/davidemms/OrthoFinder)  
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
- [trimal](http://trimal.cgenomics.org/)
- [IQ-TREE](http://www.iqtree.org/)
- [AnGST](https://github.com/almlab/angst)
- [parallel](https://www.gnu.org/software/parallel/)

## Method Summary

1. Grouping ortholog sequences using OrthoFinder
2. Align each OG(Ortholog Group) sequences using mafft
3. Trim each OG alignment using trimal
4. Reconstruct each OG gene tree using IQ-TREE
5. **Species tree** & **each OG gene tree** DTL reconciliation using AnGST
6. Aggregate and map genome-wide DTL reconciliation result

## Command Usage

### Run Command

    FastDTLmapper -i [fasta|genbank directory] -t [species tree file] -o [output directory]

### Options

    -h, --help           show this help message and exit
    -i , --indir         Input Fasta(*.fa|*.faa|*.fasta), Genbank(*.gb|*.gbk|*.genbank) directory
    -t , --tree          Input rooted species tree file (Newick format)
    -o , --outdir        Output directory
    -p , --process_num   Number of processor (Default: MaxProcessor - 1)
    --dup_cost           Duplication event cost (Default: 2)
    --los_cost           Loss event cost (Default: 1)
    --trn_cost           Transfer event cost (Default: 3)
    --inflation          MCL inflation parameter (Default: 3.0)
    --timetree           Use species tree as timetree (Default: off)
    --rseed              Number of random seed (Default: 0)

**--timetree** enable AnGST timetree option  
Following explanation of timetree is taken from [AnGST manual](<https://github.com/almlab/angst/blob/master/doc/manual.pdf>).
> If the branch lengths on the provided species tree represent times,  
> AnGST can restrict the set of possible inferred gene transfers to  
> only those between contemporaneous lineages  

### Example Command

    FastDTLmapper -i ./example/fasta/ -t ./example/species_tree.nwk -o ./fastdtlmapper_result

## Output Summary

### Output Top Directory Contents

| Top directory           | Contents                                                     |
| ----------------------- | ------------------------------------------------------------ |
| 00_user_data            | Formatted user input fasta and tree files                    |
| 01_orthofinder          | OrthoFinder raw output results                               |
| 02_dtl_reconciliation   | Each OG(Ortholog Group) DTL reconciliation result            |
| 03_aggregate_map_dtl    | Genome-wide DTL reconciliation aggregated and mapped results |
| log                     | Config log and command log files                                 |

### Output Directory Structure & Contents

    .
    ├── 00_user_data/  -- User input data
    │   ├── fasta/     -- Formatted fasta files
    │   └── tree/      -- Formatted newick species tree files
    ├── 01_orthofinder/  -- OrthoFinder raw output results
    ├── 02_dtl_reconciliation/  -- Each OG(Ortholog Group) DTL reconciliation result
    │   ├── OG0000000/
    │   │   ├── OG0000000.fa                 -- OG fasta file
    │   │   ├── OG0000000_aln.fa             -- OG alignment fasta file
    │   │   ├── OG0000000_aln_trim.fa        -- Trimmed OG alignement fasta file
    │   │   ├── OG0000000_dtl_map.nwk        -- OG DTL event mapped tree file
    │   │   ├── OG0000000_gain_loss_map.nwk  -- OG Gain-Loss event mapped tree file
    │   │   ├── angst/                       -- AnGST DTL reconciliation result
    │   │   └── iqtree/                      -- IQ-TREE gene tree reconstruction result
    │   ├── OG0000001/
    │   . 
    │   . 
    │   └── OGXXXXXXX/
    ├── 03_aggregate_map_result/  -- Genome-wide DTL reconciliation aggregated and mapped results
    │   ├── all_dtl_map.nwk              -- Genome-wide DTL event mapped tree file
    │   ├── all_gain_loss_map.nwk        -- Genome-wide Gain-Loss event mapped tree file
    │   ├── all_og_node_event.tsv        -- All OG DTL event record file
    │   ├── all_transfer_gene_count.tsv  -- All transfer(HGT) gene count file
    │   └── all_transfer_gene_list.tsv   -- All transfer(HGT) gene list file
    └── log/
        ├── parallel_cmds/ -- Parallel run command log results
        └── run_config.log -- Program run config log file
