#!/usr/bin/env python3
import datetime
import shutil
import subprocess as sp
from typing import Dict, List

from lib.config import Config, get_config
from lib.mowgli import MowgliEventMap, NodeEvent
from lib.util import UtilSeq, UtilTree


def main(config: Config):
    """Main function"""
    # TODO: Add check_user_data(config) function
    # 1. fasta count = leaf species count
    # 2. fasta file name = leaf species name
    # 3. Dependent program check

    # 00. Prepare analysis data
    format_user_fasta(config)
    format_user_tree(config)
    # 01. Grouping ortholog sequences using OrthoFinder
    orthofinder_run(config)
    # 02. Align each OG(Ortholog Group) sequences using mafft
    mafft_run(config)
    # 03. Trim each OG alignment using trimal
    trimal_run(config)
    # 04. Reconstruct each OG gene tree using FastTree
    fasttree_run(config)
    # 05. Rooting each OG gene tree using OptRoot
    optroot_run(config)
    # 06. 'Species tree' & 'each OG gene tree' DTL reconciliation using Mowgli
    mowgli_run(config)
    # 07. Aggregate and map DTL reconciliation result
    aggregate_dtl_result(config)


def format_user_tree(config: Config) -> None:
    """Format user newick tree file (Make ultrametric & Add node id)"""
    # Copy raw user input tree file
    shutil.copy(config.tree_file, config.user_tree_file)
    # Make ultrametric tree
    make_ultrametric_cmd = config.make_ultrametric_cmd(
        config.user_tree_file, config.ultrametric_tree_file
    )
    sp.run(make_ultrametric_cmd, shell=True)
    # Add internal node id (e.g. N001,N002,...N0XX)
    UtilTree.add_internal_node_id(config.ultrametric_tree_file, config.nodeid_tree_file)


def format_user_fasta(config: Config) -> None:
    """Format fasta file from fasta or genbank file"""
    infiles = [f for f in config.indir.glob("*") if f.is_file()]
    for infile in infiles:
        fasta_outfile = config.fasta_dir / infile.with_suffix(".fa").name
        filesymbol = infile.stem.replace("|", "_")
        id_prefix = f"{filesymbol}_"
        if infile.suffix in (".fa", ".faa", ".fasta"):
            UtilSeq.add_serial_id(infile, fasta_outfile, id_prefix)
        elif infile.suffix in (".gb", ".gbk", ".genbank"):
            UtilSeq.gbk2cds_fasta(infile, fasta_outfile, "protein", id_prefix)


def orthofinder_run(config: Config) -> None:
    """Run OrthoFinder"""
    print("# 01. Grouping ortholog sequences using OrthoFinder")
    # Get OrthoFinder output directory name
    weekday = datetime.date.today().strftime("%b%d")
    ortho_result_dir = config.fasta_dir / "OrthoFinder" / ("Results_" + weekday)
    # Run OrthoFinder
    orthofinder_cmd = config.orthofinder_cmd(config.fasta_dir)
    sp.run(orthofinder_cmd, shell=True)
    # Move output result & remove unnecessary directory
    shutil.rmtree(config.ortho_dir)
    config.ortho_dir.mkdir(exist_ok=True)
    for data in ortho_result_dir.glob("*"):
        shutil.move(str(data), config.ortho_dir)
    shutil.rmtree(ortho_result_dir.parent)


def mafft_run(config: Config) -> None:
    """Run mafft"""
    print("# 02. Align each OG(Ortholog Group) sequences using mafft")
    mafft_cmd_list = []
    for fasta_file in sorted(config.ortho_group_fasta_dir.glob("*.fa")):
        seq_count = UtilSeq.count_fasta_seq(fasta_file)
        aln_file = config.ortho_aln_dir / fasta_file.name
        if seq_count >= 2:
            mafft_cmd = config.mafft_cmd(fasta_file, aln_file)
            mafft_cmd_list.append(mafft_cmd)
        else:
            shutil.copy(fasta_file, aln_file)

    mafft_parallel_cmd = config.parallel_cmd(mafft_cmd_list)
    sp.run(mafft_parallel_cmd, shell=True)


def trimal_run(config: Config) -> None:
    """Run trimal"""
    print("# 03. Trim each OG alignment using trimal")
    trimal_cmd_list = []
    for aln_file in sorted(config.ortho_aln_dir.glob("*.fa")):
        seq_count = UtilSeq.count_fasta_seq(aln_file)
        aln_trim_file = config.ortho_aln_trim_dir / aln_file.name
        shutil.copy(aln_file, aln_trim_file)
        if seq_count >= 2:
            trimal_cmd = config.trimal_cmd(aln_file, aln_trim_file)
            trimal_cmd_list.append(trimal_cmd)

    trimal_parallel_cmd = config.parallel_cmd(trimal_cmd_list)
    sp.run(trimal_parallel_cmd, shell=True)


def fasttree_run(config: Config) -> None:
    """Run FastTree"""
    print("# 04. Reconstruct each OG gene tree using FastTree")
    fasttree_cmd_list = []
    for aln_trim_file in sorted(config.ortho_aln_trim_dir.glob("*.fa")):
        seq_count = UtilSeq.count_fasta_seq(aln_trim_file)
        gene_tree_file = config.gene_tree_dir / (aln_trim_file.stem + ".nwk")
        if seq_count >= 3:
            fasttree_cmd = config.fasttree_cmd(aln_trim_file, gene_tree_file)
            fasttree_cmd_list.append(fasttree_cmd)

    fasttree_parallel_cmd = config.parallel_cmd(fasttree_cmd_list)
    sp.run(fasttree_parallel_cmd, shell=True)


def optroot_run(config: Config) -> None:
    """Run OptRoot"""
    print("# 05. Rooting each OG gene tree using OptRoot")
    optroot_cmd_list = []
    optroot_tmp_infile_list, optroot_outfile_list = [], []
    for gene_tree_file in sorted(config.gene_tree_dir.glob("*.nwk")):
        UtilTree.convert_bootstrap_float_to_int(gene_tree_file, gene_tree_file)
        optroot_outfile = config.rooted_gene_tree_dir / (gene_tree_file.stem + ".txt")
        optroot_tmp_infile = config.rooted_gene_tree_dir / ("tmp" + gene_tree_file.name)
        concat_cmd = (
            f"cat {config.ultrametric_tree_file} {gene_tree_file} "
            + f"> {optroot_tmp_infile}"
        )
        sp.run(concat_cmd, shell=True)
        optroot_cmd = config.optroot_cmd(optroot_tmp_infile, optroot_outfile)

        optroot_outfile_list.append(optroot_outfile)
        optroot_tmp_infile_list.append(optroot_tmp_infile)
        optroot_cmd_list.append(optroot_cmd)

    optroot_parallel_cmd = config.parallel_cmd(optroot_cmd_list)
    sp.run(optroot_parallel_cmd, shell=True)
    # Delete tmp files
    [tmp_file.unlink() for tmp_file in optroot_tmp_infile_list]

    # Extract rooted gene tree content from optroot result
    for optroot_outfile in optroot_outfile_list:
        rooted_gene_tree_file = optroot_outfile.with_suffix(".nwk")
        extract_cmd = f"cat {optroot_outfile} | grep ';' > {rooted_gene_tree_file}"
        sp.run(extract_cmd, shell=True)
    # Delete unnecessary optroot result files
    [out_file.unlink() for out_file in optroot_outfile_list]


def mowgli_run(config: Config) -> None:
    """Run Mowgli"""
    print("# 06. 'Species tree' & 'each OG gene tree' DTL reconciliation using Mowgli")
    mowgli_cmd_list = []
    for rooted_gene_tree_file in sorted(config.rooted_gene_tree_dir.glob("*.nwk")):
        group_id = rooted_gene_tree_file.stem
        outdir = config.dtl_dir / group_id
        outdir.mkdir(exist_ok=True)
        mowgli_cmd = config.mowgli_cmd(
            config.ultrametric_tree_file, rooted_gene_tree_file, outdir
        )
        mowgli_cmd_list.append(mowgli_cmd)

    mowgli_parallel_cmd = config.parallel_cmd(mowgli_cmd_list)
    sp.run(mowgli_parallel_cmd, shell=True)


def aggregate_dtl_result(config: Config) -> Dict[str, List[NodeEvent]]:
    """Aggregate Mowgli DTL event results"""
    print("# 07. Aggregate and map DTL reconciliation result")
    group_id2all_node_event: Dict[str, List[NodeEvent]] = {}
    for dtl_result_dir in sorted(config.dtl_dir.glob("*")):
        # TODO: SeqNum <= 2 DTL reconciliation and make node event
        mowgli_xml_file = dtl_result_dir / "RecGeneTree.xml"
        mowgli_species_tree_file = dtl_result_dir / "outputSpeciesTree.mpr"
        event_map = MowgliEventMap(mowgli_species_tree_file, mowgli_xml_file)

        group_id = dtl_result_dir.name
        dtl_group_dir = config.dtl_event_map_dir / group_id
        dtl_group_dir.mkdir(exist_ok=True)
        gain_loss_nwk_file = dtl_group_dir / (group_id + "_gain_loss.nwk")
        dtl_nwk_file = dtl_group_dir / (group_id + "_dtl.nwk")
        event_map.write_tree(gain_loss_nwk_file, "gain-loss")
        event_map.write_tree(dtl_nwk_file, "dtl")

        all_node_event = event_map.get_all_node_event()
        group_id2all_node_event[group_id] = all_node_event

    output_content = ""
    for group_id, all_node_event in group_id2all_node_event.items():
        for node_event in all_node_event:
            output_content += f"{group_id}\t{node_event.as_tsv_format}\n"

    with open(config.outdir / "all_group_node_event.tsv", "w") as f:
        f.write(
            "OG_ID\tNODE_ID\tGENE_NUM\tGAIN_NUM\tBRN_NUM\tDUP_NUM\tTRN_NUM"
            + "\tLOS_NUM\tTRN_DETAIL\n"
        )
        f.write(output_content)


if __name__ == "__main__":
    config = get_config()
    main(config)
