# (Under Construction) FastDTLmapper: Fast genome-wide DTL event mapper  

![Python3](https://img.shields.io/badge/Language-Python_3.7_|_3.8_|_3.9-steelblue)
![OS](https://img.shields.io/badge/OS-Linux-steelblue)
![License](https://img.shields.io/badge/License-GPL3.0-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/fastdtlmapper.svg)](https://pypi.python.org/pypi/fastdtlmapper)
[![Downloads](https://static.pepy.tech/personalized-badge/fastdtlmapper?period=month&units=international_system&left_color=grey&right_color=orange&left_text=Downloads)](https://pepy.tech/project/fastdtlmapper)  
![CI workflow](https://github.com/moshi4/FastDTLmapper/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/moshi4/FastDTLmapper/branch/main/graph/badge.svg?token=ZJ8D747JUY)](https://codecov.io/gh/moshi4/FastDTLmapper)

## Table of contents

- [Overview](#overview)
- [Install](#install)
- [Analysis Pipeline](#analysis-pipeline)
- [Command Usage](#command-usage)
- [Output Contents](#output-contents)
- [Further Analysis](#further-analysis)

## Overview

Gene gain/loss is considered to be one of the most important evolutionary processes
driving adaptive evolution, but it remains largely unexplored.
Therefore, to investigate the relationship between gene gain/loss and adaptive evolution
in the evolutionary process of organisms, I developed a software pipeline **FastDTLmapper**
which automatically estimates and maps genome-wide gene gain/loss.  
FastDTLmapper takes two inputs, 1. *Species tree (Newick format)* 2. *Genomic Protein CDSs (Fasta|Genbank format)*,
and performs genome-wide mapping of DTL(Duplication-Transfer-Loss) events by
DTL reconciliation of species tree and gene trees.  
Additionally, FastDTLmapper can perform
[Plot Gain/Loss Map Figure](#plot-gainloss-map-figure) and
[Functional Analysis (GOEA)](#functional-analysis-goea)
using packaged subtools.  

![demo_all_gain_loss_map.png](https://github.com/moshi4/FastDTLmapper/wiki/images/demo_all_gain_loss_map.png)  
**Fig. Genome-wide gain/loss map result example (all_gain_loss_map.nwk)**  
Each node gain/loss data is mapped in following format (*NodeID | GeneCount [gain=GainCount los=LossCount]*)  
Map data is embeded in newick format bootstrap value field and user can visualize using [SeaView](http://doua.prabi.fr/software/seaview).  

## Install

FastDTLmapper is implemented in **Python3(>=3.7)** and runs on **Linux** (Tested on Ubuntu20.04).  

Install PyPI stable version with pip:

    pip install fastdtlmapper

Install latest development version with pip:

    pip install git+git://github.com/moshi4/FastDTLmapper.git

### Dependencies

Python package dependencies list here (Auto installed with pip).

Well known python package `numpy`, `pandas`, `scipy` and

- [BioPython](https://github.com/biopython/biopython)  
  Utility tools for computational molecular biology
- [GOAtools](https://github.com/tanghaibao/goatools)  
  GOEA(GO Enrichment Analysis) tool
- [ETE3](http://etetoolkit.org/)  
  Tree analysis and visualization tool

Following dependencies are packaged in **[src/fastdtlmapper/bin](https://github.com/moshi4/FastDTLmapper/tree/main/src/fastdtlmapper/bin)** directory.  

- [OrthoFinder](https://github.com/davidemms/OrthoFinder)  
  Orthology inference tool
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
  Sequences alignment tool
- [trimal](http://trimal.cgenomics.org/)  
  Alignment sequences trim tool
- [IQ-TREE](http://www.iqtree.org/)  
  Phylogenetic tree reconstruction tool
- [AnGST](https://github.com/almlab/angst)  
  DTL reconciliation tool (Requires Python 2.X to run)
- [parallel](https://www.gnu.org/software/parallel/)  
  Job parallelization tool (Requires Perl to run)

## Analysis Pipeline

This is summary description. See [wiki]() for details.

1. Grouping ortholog sequences using OrthoFinder
2. Align each OG(Ortholog Group) sequences using mafft
3. Trim each OG alignment using trimal
4. Reconstruct each OG gene tree using IQ-TREE
5. DTL reconciliation of species tree & each OG gene tree using AnGST
6. Aggregate and map genome-wide DTL reconciliation result

## Command Usage

### Basic Run Command

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

#### Timetree Option

If user set this option, input species tree must be ultrametric tree.  
--timetree enable AnGST timetree option below (See [AnGST manual](<https://github.com/almlab/angst/blob/master/doc/manual.pdf>) for details).  
> If the branch lengths on the provided species tree represent times,
> AnGST can restrict the set of possible inferred gene transfers to
> only those between contemporaneous lineages  

#### Input Limitation

- fasta or genbank files (--indir option)  
  :warning: Following characters cannot be included in file name '_', '-', '|', '.'
- species tree file (--tree option)  
  :warning: Species name in species tree must match fasta or genbank file name

### Example Command

    FastDTLmapper -i ./example/fasta/ -t ./example/species_tree.nwk -o ./fastdtlmapper_result

## Output Contents

### Output Top Directory

| Top directory           | Contents                                                     |
| ----------------------- | ------------------------------------------------------------ |
| 00_user_data            | Formatted user input fasta and tree files                    |
| 01_orthofinder          | OrthoFinder raw output results                               |
| 02_dtl_reconciliation   | Each OG(Ortholog Group) DTL reconciliation result            |
| 03_aggregate_map_result | Genome-wide DTL reconciliation aggregated and mapped results |
| log                     | Config log and command log files                                 |

### Output Directory Structure & Files

    .
    ├── 00_user_data/  -- User input data
    │   ├── fasta/     -- Formatted fasta files
    │   └── tree/      -- Formatted newick species tree files
    │
    ├── 01_orthofinder/  -- OrthoFinder raw output results
    │
    ├── 02_dtl_reconciliation/  -- Each OG(Ortholog Group) DTL reconciliation result
    │   ├── OG0000000/
    │   │   ├── OG0000000.fa                 -- OG fasta file
    │   │   ├── OG0000000_aln.fa             -- OG alignment fasta file
    │   │   ├── OG0000000_aln_trim.fa        -- Trimmed OG alignement fasta file
    │   │   ├── OG0000000_dtl_map.nwk        -- OG DTL event mapped tree file
    │   │   ├── OG0000000_gain_loss_map.nwk  -- OG Gain-Loss event mapped tree file
    │   │   ├── angst/                       -- AnGST DTL reconciliation result
    │   │   └── iqtree/                      -- IQ-TREE gene tree reconstruction result
    │   │
    │   ├── OG0000001/
    │   . 
    │   . 
    │   └── OGXXXXXXX/
    │
    ├── 03_aggregate_map_result/  -- Genome-wide DTL reconciliation aggregated and mapped results
    │   ├── all_dtl_map.nwk              -- Genome-wide DTL event mapped tree file
    │   ├── all_gain_loss_map.nwk        -- Genome-wide Gain-Loss event mapped tree file
    │   ├── all_og_node_event.tsv        -- All OG DTL event record file
    │   ├── all_transfer_gene_count.tsv  -- All transfer gene count file
    │   └── all_transfer_gene_list.tsv   -- All transfer gene list file
    │
    └── log/
        ├── parallel_cmds/ -- Parallel run command log results
        └── run_config.log -- Program run config log file

## Further Analysis

### Plot Gain/Loss Map Figure

FastDTLmapper subtool `plot_gain_loss_map` supports for plotting
publication-ready gain/loss map figure as shown below.
Users can plot easily and can output in any format they want by
changing plotting parameter.
See [wiki]() for details.

![demo_plot_all_gain_loss_map.png](https://github.com/moshi4/FastDTLmapper/wiki/images/demo_plot_all_gain_loss_map.png)  
**Fig. Gain/Loss map plot result example**

### Functional Analysis (GOEA)

FastDTLmapper subtool `FastDTLgoea` supports for performing
GOEA(GO Enrichment Analysis) in each node gain/loss genes.
Each node gain/loss gene's significant GOterms are
listed and plotted as shown below.
This GOEA functional analysis is useful for getting glasp of genome-wide
functional trends in genome-wide gain/loss genes. See [wiki]() for details.

![demo_plot_goea.png](https://github.com/moshi4/FastDTLmapper/wiki/images/demo_plot_goea.png)  
**Fig. GOEA plot result example**  
In this example, gain gene's significant over representation
GOterms in N023 node is plotted with color.  
CC indicates GO category of 'Cell Components'. 3 GO category BP,MF,CC exists.
