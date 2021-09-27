#!/usr/bin/env python3
import datetime
import shutil
import subprocess as sp
from typing import Dict, List

from lib.config import Config, get_config
from lib.input_check import InputCheck
from lib.mowgli import MowgliEventMap, NodeEvent
from lib.reconcilation import Reconciliation as Rec
from lib.util import UtilSeq, UtilTree


def main(config: Config):
    """Main function"""
    # Check user input
    InputCheck(config.indir, config.tree_file).run()
    # # 00. Prepare analysis data
    format_user_tree(config)
    format_user_fasta(config)
    # # 01. Grouping ortholog sequences using OrthoFinder
    orthofinder_run(config)
    # # 02. Align each OG(Ortholog Group) sequences using mafft
    mafft_run(config)
    # # 03. Trim each OG alignment using trimal
    trimal_run(config)
    # # 04. Reconstruct each OG gene tree using FastTree
    fasttree_run(config)
    # # 05. Rooting each OG gene tree using OptRoot
    optroot_run(config)
    # # 06. 'Species tree' & 'each OG gene tree' DTL reconciliation using Mowgli
    mowgli_run(config)
    # 07. Aggregate and map DTL reconciliation result
    group_id2all_node_event = aggregate_dtl_results(config)
    output_aggregate_map_results(config, group_id2all_node_event)


def format_user_tree(config: Config) -> None:
    """Format user newick tree file (Make ultrametric & Add node id)"""
    # Copy raw user input tree file
    shutil.copy(config.tree_file, config.user_tree_file)
    # Make ultrametric tree
    make_ultrametric_cmd = config.make_ultrametric_cmd(
        config.user_tree_file, config.um_tree_file
    )
    sp.run(make_ultrametric_cmd, shell=True)
    # Add internal node id (e.g. N001,N002,...N0XX)
    UtilTree.add_internal_node_id(config.um_tree_file, config.um_nodeid_tree_file)


def format_user_fasta(config: Config) -> None:
    """Format fasta file from fasta or genbank file"""
    infiles = [f for f in config.indir.glob("*") if f.is_file()]
    for infile in infiles:
        fasta_outfile = config.user_fasta_dir / infile.with_suffix(".fa").name
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
    ortho_result_dir = config.ortho_tmpwork_dir / ("Results_" + weekday)
    # Run OrthoFinder
    orthofinder_cmd = config.orthofinder_cmd(config.user_fasta_dir)
    sp.run(orthofinder_cmd, shell=True)
    # Move output result & remove unnecessary directory
    shutil.rmtree(config.ortho_dir)
    config.ortho_dir.mkdir(exist_ok=True)
    for data in ortho_result_dir.glob("*"):
        shutil.move(str(data), config.ortho_dir)
    shutil.rmtree(ortho_result_dir.parent)

    # Copy OrthologGroup fasta to DTL reconciliation analysis directory
    for fasta_file in sorted(config.ortho_group_fasta_dir.glob("*.fa")):
        group_id = fasta_file.stem
        dtl_rec_group_dir = config.dtl_rec_dir / group_id
        dtl_rec_group_dir.mkdir(exist_ok=True)
        shutil.copy(fasta_file, dtl_rec_group_dir)


def mafft_run(config: Config) -> None:
    """Run mafft"""
    print("# 02. Align each OG(Ortholog Group) sequences using mafft")
    mafft_cmd_list = []
    for dtl_rec_group_dir in sorted(config.dtl_rec_dir.glob("*")):
        group_id = dtl_rec_group_dir.name
        fasta_file = dtl_rec_group_dir / (group_id + ".fa")
        aln_file = fasta_file.with_name(group_id + "_aln.fa")
        if UtilSeq.count_fasta_seq(fasta_file) >= 3:
            mafft_cmd = config.mafft_cmd(fasta_file, aln_file)
            mafft_cmd_list.append(mafft_cmd)

    mafft_parallel_cmd = config.parallel_cmd(mafft_cmd_list, "mafft")
    sp.run(mafft_parallel_cmd, shell=True)


def trimal_run(config: Config) -> None:
    """Run trimal"""
    print("# 03. Trim each OG alignment using trimal")
    trimal_cmd_list = []
    for aln_file in sorted(config.dtl_rec_dir.glob("**/*_aln.fa")):
        group_id = aln_file.parent.name
        aln_trim_file = aln_file.parent / (group_id + "_aln_trim.fa")
        trimal_cmd = config.trimal_cmd(aln_file, aln_trim_file)
        trimal_cmd_list.append(trimal_cmd)
        shutil.copy(aln_file, aln_trim_file)

    trimal_parallel_cmd = config.parallel_cmd(trimal_cmd_list, "trimal")
    sp.run(trimal_parallel_cmd, shell=True)

    # If trimal removed sequence, use raw mafft alignment in next step
    for aln_file in sorted(config.dtl_rec_dir.glob("**/*_aln.fa")):
        group_id = aln_file.parent.name
        aln_trim_file = aln_file.parent / (group_id + "_aln_trim.fa")
        aln_seq_num = UtilSeq.count_fasta_seq(aln_file)
        aln_trim_seq_num = UtilSeq.count_fasta_seq(aln_trim_file)
        if aln_seq_num > aln_trim_seq_num:
            shutil.copy(aln_file, aln_trim_file)


def fasttree_run(config: Config) -> None:
    """Run FastTree"""
    print("# 04. Reconstruct each OG gene tree using FastTree")
    fasttree_cmd_list = []
    for aln_trim_file in sorted(config.dtl_rec_dir.glob("**/*_aln_trim.fa")):
        group_id = aln_trim_file.parent.name
        gene_tree_file = aln_trim_file.parent / (group_id + "_unrooted_tree.nwk")
        fasttree_cmd = config.fasttree_cmd(aln_trim_file, gene_tree_file)
        fasttree_cmd_list.append(fasttree_cmd)

    fasttree_parallel_cmd = config.parallel_cmd(fasttree_cmd_list, "FastTree")
    sp.run(fasttree_parallel_cmd, shell=True)


def optroot_run(config: Config) -> None:
    """Run OptRoot"""
    print("# 05. Rooting each OG gene tree using OptRoot")
    optroot_cmd_list = []
    optroot_tmp_infile_list, optroot_outfile_list = [], []
    for gene_tree_file in sorted(config.dtl_rec_dir.glob("**/*_unrooted_tree.nwk")):
        group_id = gene_tree_file.parent.name
        # OptRoot cannot handle float bootstrap value
        UtilTree.convert_bootstrap_float_to_int(gene_tree_file, gene_tree_file)
        # Make OptRoot input tmpfile (species tree & gene tree mixed newick file)
        optroot_tmp_infile = gene_tree_file.parent / ("tmp" + gene_tree_file.name)
        concat_cmd = (
            f"cat {config.um_tree_file} {gene_tree_file} > {optroot_tmp_infile}"
        )
        sp.run(concat_cmd, shell=True)

        optroot_outfile = gene_tree_file.parent / (gene_tree_file.stem + ".txt")
        optroot_cmd = config.optroot_cmd(optroot_tmp_infile, optroot_outfile)

        optroot_outfile_list.append(optroot_outfile)
        optroot_tmp_infile_list.append(optroot_tmp_infile)
        optroot_cmd_list.append(optroot_cmd)

    optroot_parallel_cmd = config.parallel_cmd(optroot_cmd_list, "OptRoot")
    sp.run(optroot_parallel_cmd, shell=True)
    # Delete tmp files
    [tmp_file.unlink() for tmp_file in optroot_tmp_infile_list]

    # Extract rooted gene tree content from optroot result
    for optroot_outfile in optroot_outfile_list:
        group_id = optroot_outfile.parent.stem
        rooted_gene_tree_file = optroot_outfile.with_name(f"{group_id}_rooted_tree.nwk")
        extract_cmd = f"cat {optroot_outfile} | grep ';' > {rooted_gene_tree_file}"
        sp.run(extract_cmd, shell=True)
    # Delete unnecessary optroot result files
    [out_file.unlink() for out_file in optroot_outfile_list]


def mowgli_run(config: Config) -> None:
    """Run Mowgli"""
    print("# 06. 'Species tree' & 'each OG gene tree' DTL reconciliation using Mowgli")
    mowgli_cmd_list = []
    for rooted_tree_file in sorted(config.dtl_rec_dir.glob("**/*_rooted_tree.nwk")):
        outdir = rooted_tree_file.parent / "dtl_reconciliation"
        outdir.mkdir(exist_ok=True)
        mowgli_cmd = config.mowgli_cmd(config.um_tree_file, rooted_tree_file, outdir)
        mowgli_cmd_list.append(mowgli_cmd)

    mowgli_parallel_cmd = config.parallel_cmd(mowgli_cmd_list, "Mowgli")
    sp.run(mowgli_parallel_cmd, shell=True)


def aggregate_dtl_results(config: Config) -> Dict[str, List[NodeEvent]]:
    """Aggregate DTL reconciliation result"""
    print("# 07. Aggregate and map DTL reconciliation result")
    group_id2node_event_list: Dict[str, List[NodeEvent]] = {}
    for dtl_rec_group_dir in sorted(config.dtl_rec_dir.glob("*")):
        group_id = dtl_rec_group_dir.name
        group_fasta_file = dtl_rec_group_dir / (group_id + ".fa")
        seq_count = UtilSeq.count_fasta_seq(group_fasta_file)
        if seq_count >= 3:
            # Get DTL reconciliation event map
            dtl_rec_group_result_dir = dtl_rec_group_dir / "dtl_reconciliation"
            mowgli_xml_file = dtl_rec_group_result_dir / "RecGeneTree.xml"
            mowgli_tree_file = dtl_rec_group_result_dir / "outputSpeciesTree.mpr"
            event_map = MowgliEventMap(mowgli_tree_file, mowgli_xml_file)

            # Map each group DTL event to newick tree
            gain_loss_map_file = dtl_rec_group_dir / (group_id + "_gain_loss_map.nwk")
            dtl_map_file = dtl_rec_group_dir / (group_id + "_dtl_map.nwk")
            event_map.write_tree(gain_loss_map_file, "gain-loss")
            event_map.write_tree(dtl_map_file, "dtl")

            # Get each group all node DTL event
            node_event_list = event_map.get_all_node_event()

        elif seq_count == 2:
            node_event_list = Rec.two_species(
                config.um_nodeid_tree_file,
                group_fasta_file,
                config.los_cost,
                config.trn_cost,
            )
        elif seq_count == 1:
            node_event_list = Rec.one_species(
                config.um_nodeid_tree_file, group_fasta_file
            )

        group_id2node_event_list[group_id] = node_event_list

        # DTL result consistency check
        species_name2seq_num = UtilSeq.get_species_name2seq_num(group_fasta_file)
        for species_name, seq_num in species_name2seq_num.items():
            for node_event in node_event_list:
                # Species 'seq_num' must be equal node event 'gene_num'
                if (
                    species_name == node_event.node_id
                    and seq_num != node_event.gene_num
                ):
                    print(
                        "# DTL ConsistencyCheckError:\n"
                        + f"{group_id} {species_name} seq_num & gene_num not equal!!\n"
                        + f"seq_num = {seq_num}, gene_num = {node_event.gene_num}",
                    )

    return group_id2node_event_list


def output_aggregate_map_results(
    config: Config, group_id2node_event_list: Dict[str, List[NodeEvent]]
):
    """Output aggregated results and mapped dtl reconciliation results"""
    # Output results in TSV format
    output_content = ""
    for group_id, node_event_list in group_id2node_event_list.items():
        for node_event in node_event_list:
            output_content += f"{group_id}\t{node_event.as_tsv_format}\n"
    with open(config.all_node_event_file, "w") as f:
        header = (
            "OG_ID\tNODE_ID\tGENE_NUM\tGAIN_NUM\tBRN_NUM\tDUP_NUM\t"
            + "TRN_NUM\tLOS_NUM\tTRN_DETAIL"
        )
        f.write(header + "\n")
        f.write(output_content)

    # Output results in newick mapped format
    node_id2total_node_event: Dict[str, NodeEvent] = {}
    for node_event_list in group_id2node_event_list.values():
        # Count up the number of events for each node id
        for node_event in node_event_list:
            node_id = node_event.node_id
            if node_id not in node_id2total_node_event.keys():
                node_id2total_node_event[node_id] = NodeEvent(node_id, gene_num=0)
            node_id2total_node_event[node_id] += node_event

    UtilTree.map_node_event(
        config.um_nodeid_tree_file,
        node_id2total_node_event,
        config.all_gain_loss_map_file,
        "gain-loss",
    )
    UtilTree.map_node_event(
        config.um_nodeid_tree_file,
        node_id2total_node_event,
        config.all_dtl_map_file,
        "dtl",
    )


if __name__ == "__main__":
    config = get_config()
    main(config)
    config.output_run_config_log()
