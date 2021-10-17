#!/usr/bin/env python3
import datetime
import shutil
import subprocess as sp
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from fastdtlmapper.angst import AngstEventMap, AngstTransferGene, NodeEvent
from fastdtlmapper.args import Args, get_args
from fastdtlmapper.cmd import Cmd
from fastdtlmapper.input_check import InputCheck
from fastdtlmapper.out_path import OutPath
from fastdtlmapper.reconcilation import Reconciliation as Rec
from fastdtlmapper.setup_binpath import SetupBinpath
from fastdtlmapper.util import UtilFasta, UtilGenbank, UtilTree


def main(args: Optional[Args] = None):
    """Main function"""
    # Get configs (arguments, output_path, command)
    if args is None:
        args = get_args()
    outpath = OutPath(args.outdir)
    cmd = Cmd(args)

    # Setup binary path
    root_binpath = Path(__file__).parent.parent / "bin"
    SetupBinpath(root_binpath)

    # Check user input
    InputCheck(args.indir, args.tree_file).run()
    # # 00. Prepare analysis data
    format_user_tree(args.tree_file, outpath, cmd)
    format_user_fasta(args.indir, outpath.user_fasta_dir)
    # 01. Grouping ortholog sequences using OrthoFinder
    orthofinder_run(outpath, cmd)
    # 02. Align each OG(Ortholog Group) sequences using mafft
    mafft_run(outpath, cmd)
    # 03. Trim each OG alignment using trimal
    trimal_run(outpath, cmd)
    # 04. Reconstruct each OG gene tree using iqtree
    iqtree_run(outpath, cmd)
    # 05. DTL reconciliation using AnGST
    angst_run(outpath, cmd)
    # 06. Aggregate and map DTL reconciliation result
    group_id2all_node_event = aggregate_dtl_results(args, outpath)
    output_aggregate_map_results(outpath, group_id2all_node_event)
    output_aggregate_transfer_results(outpath)

    args.write_log(outpath.run_config_log_file)


def format_user_tree(tree_file: Path, outpath: OutPath, cmd: Cmd) -> None:
    """Format user newick tree file (Make ultrametric & Add node id)"""
    # Copy raw user input tree file
    shutil.copy(tree_file, outpath.user_tree_file)
    # Make ultrametric tree
    make_ultrametric_cmd = cmd.get_make_ultrametric_cmd(
        outpath.user_tree_file, outpath.um_tree_file
    )
    sp.run(make_ultrametric_cmd, shell=True)
    # Add internal node id (e.g. N001,N002,...N0XX)
    UtilTree.add_internal_node_id(outpath.um_tree_file, outpath.um_nodeid_tree_file)


def format_user_fasta(fasta_dir: Path, format_fasta_dir: Path) -> None:
    """Format fasta file from fasta or genbank file"""
    infiles = [f for f in fasta_dir.glob("*") if f.is_file()]
    for infile in infiles:
        fasta_outfile = format_fasta_dir / infile.with_suffix(".fa").name
        filesymbol = infile.stem.replace("|", "_")
        id_prefix = f"{filesymbol}_GENE"
        if infile.suffix in (".fa", ".faa", ".fasta"):
            UtilFasta(infile).add_serial_id(fasta_outfile, id_prefix)
        elif infile.suffix in (".gb", ".gbk", ".genbank"):
            UtilGenbank(infile).convert_cds_fasta(fasta_outfile, "protein", id_prefix)


def orthofinder_run(outpath: OutPath, cmd: Cmd) -> None:
    """Run OrthoFinder"""
    print("\n# 01. Grouping ortholog sequences using OrthoFinder")
    # Get OrthoFinder output directory name
    weekday = datetime.date.today().strftime("%b%d")
    ortho_result_dir = outpath.ortho_tmpwork_dir / ("Results_" + weekday)
    # Run OrthoFinder
    orthofinder_cmd = cmd.get_orthofinder_cmd(outpath.user_fasta_dir)
    sp.run(orthofinder_cmd, shell=True)
    # Move output result & remove unnecessary directory
    shutil.rmtree(outpath.ortho_dir)
    outpath.ortho_dir.mkdir(exist_ok=True)
    for data in ortho_result_dir.glob("*"):
        shutil.move(str(data), outpath.ortho_dir)
    shutil.rmtree(ortho_result_dir.parent)

    # Copy OrthologGroup fasta to DTL reconciliation analysis directory
    for fasta_file in sorted(outpath.ortho_group_fasta_dir.glob("*.fa")):
        group_id = fasta_file.stem
        dtl_rec_group_dir = outpath.dtl_rec_dir / group_id
        dtl_rec_group_dir.mkdir(exist_ok=True)
        shutil.copy(fasta_file, dtl_rec_group_dir)


def mafft_run(outpath: OutPath, cmd: Cmd) -> None:
    """Run mafft"""
    print("\n# 02. Align each OG(Ortholog Group) sequences using mafft")
    mafft_cmd_list = []
    for dtl_rec_group_dir in sorted(outpath.dtl_rec_dir.glob("*")):
        group_id = dtl_rec_group_dir.name
        fasta_file = dtl_rec_group_dir / (group_id + ".fa")
        aln_file = fasta_file.with_name(group_id + "_aln.fa")
        if UtilFasta(fasta_file).seq_num >= 3:
            mafft_cmd = cmd.get_mafft_cmd(fasta_file, aln_file)
            mafft_cmd_list.append(mafft_cmd)

    cmd.run_parallel_cmd(
        mafft_cmd_list, outpath.tmp_parallel_cmds_file, outpath.mafft_output_log_file
    )


def trimal_run(outpath: OutPath, cmd: Cmd) -> None:
    """Run trimal"""
    print("\n# 03. Trim each OG alignment using trimal")
    trimal_cmd_list = []
    for aln_file in sorted(outpath.dtl_rec_dir.glob("**/*_aln.fa")):
        group_id = aln_file.parent.name
        aln_trim_file = aln_file.parent / (group_id + "_aln_trim.fa")
        trimal_cmd = cmd.get_trimal_cmd(aln_file, aln_trim_file)
        trimal_cmd_list.append(trimal_cmd)
        shutil.copy(aln_file, aln_trim_file)

    cmd.run_parallel_cmd(
        trimal_cmd_list, outpath.tmp_parallel_cmds_file, outpath.trimal_output_log_file
    )

    # If trimal removed sequence, use raw mafft alignment in next step
    for aln_file in sorted(outpath.dtl_rec_dir.glob("**/*_aln.fa")):
        group_id = aln_file.parent.name
        aln_trim_file = aln_file.parent / (group_id + "_aln_trim.fa")
        if UtilFasta(aln_file).seq_num > UtilFasta(aln_trim_file).seq_num:
            shutil.copy(aln_file, aln_trim_file)


def iqtree_run(outpath: OutPath, cmd: Cmd) -> None:
    """Run iqtree"""
    print("\n# 04. Reconstruct each OG gene tree using iqtree")
    iqtree_cmd_list = []
    for aln_trim_file in sorted(outpath.dtl_rec_dir.glob("**/*_aln_trim.fa")):
        group_id = aln_trim_file.parent.name
        iqtree_outdir = aln_trim_file.parent / "iqtree"
        iqtree_outdir.mkdir(exist_ok=True)
        iqtree_prefix = iqtree_outdir / group_id
        seq_num = UtilFasta(aln_trim_file).seq_num
        uniqseq_num = UtilFasta(aln_trim_file).uniq_seq_num
        if seq_num >= 4 and uniqseq_num >= 4:
            iqtree_cmd = cmd.get_iqtree_cmd(aln_trim_file, iqtree_prefix)
            iqtree_cmd_list.append(iqtree_cmd)
        elif seq_num >= 4 and uniqseq_num < 4:
            # iqtree bootstrap option cannot handle less than 4 identical sequences
            iqtree_cmd = cmd.get_iqtree_cmd(aln_trim_file, iqtree_prefix, boot=False)
            sp.run(iqtree_cmd, shell=True)
            shutil.copy(
                str(iqtree_prefix) + ".treefile", str(iqtree_prefix) + ".ufboot"
            )
        else:
            # iqtree cannot handle 3 genes tree
            gene_tree_file = str(iqtree_prefix) + ".ufboot"
            UtilTree.make_3genes_tree(aln_trim_file, gene_tree_file)

    cmd.run_parallel_cmd(
        iqtree_cmd_list, outpath.tmp_parallel_cmds_file, outpath.iqtree_output_log_file
    )


def angst_run(outpath: OutPath, cmd: Cmd) -> None:
    """Run AnGST"""
    print("\n# 05. DTL reconciliation using AnGST")
    angst_cmd_list = []
    for boot_tree_file in sorted(outpath.dtl_rec_dir.glob("**/*.ufboot")):
        outdir = boot_tree_file.parent.parent / "angst"
        outdir.mkdir(exist_ok=True)
        angst_cmd = cmd.get_angst_cmd(outpath.um_tree_file, boot_tree_file, outdir)
        angst_cmd_list.append(angst_cmd)

    cmd.run_parallel_cmd(
        angst_cmd_list, outpath.tmp_parallel_cmds_file, outpath.angst_output_log_file
    )


def aggregate_dtl_results(args: Args, outpath: OutPath) -> Dict[str, List[NodeEvent]]:
    """Aggregate DTL reconciliation result"""
    print("\n# 06. Aggregate and map DTL reconciliation result")
    group_id2node_event_list: Dict[str, List[NodeEvent]] = {}
    for dtl_rec_group_dir in sorted(outpath.dtl_rec_dir.glob("*")):
        group_id = dtl_rec_group_dir.name
        group_fasta_file = dtl_rec_group_dir / (group_id + ".fa")
        seq_count = UtilFasta(group_fasta_file).seq_num
        if seq_count >= 3:
            # Get DTL reconciliation event map
            angst_result_dir = dtl_rec_group_dir / "angst"
            event_map = AngstEventMap(outpath.um_nodeid_tree_file, angst_result_dir)

            # Map each group DTL event to newick tree
            gain_loss_map_file = dtl_rec_group_dir / (group_id + "_gain_loss_map.nwk")
            dtl_map_file = dtl_rec_group_dir / (group_id + "_dtl_map.nwk")
            event_map.write_tree(gain_loss_map_file, "gain-loss")
            event_map.write_tree(dtl_map_file, "dtl")

            # Get each group all node DTL event
            node_event_list = list(event_map.nodeid2node_event.values())

        elif seq_count == 2:
            node_event_list = Rec.two_species(
                outpath.um_nodeid_tree_file,
                group_fasta_file,
                args.los_cost,
                args.trn_cost,
            )
        elif seq_count == 1:
            node_event_list = Rec.one_species(
                outpath.um_nodeid_tree_file, group_fasta_file
            )
        else:
            raise ValueError("Unexpected fasta file detected!!")

        group_id2node_event_list[group_id] = node_event_list

        # DTL result consistency check
        species_name2seq_num = UtilFasta(group_fasta_file).species_name2seq_num
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
    outpath: OutPath, group_id2node_event_list: Dict[str, List[NodeEvent]]
):
    """Output aggregated results and mapped dtl reconciliation results"""
    # Output results in TSV format
    output_content = ""
    for group_id, node_event_list in group_id2node_event_list.items():
        for node_event in node_event_list:
            output_content += f"{group_id}\t{node_event.as_tsv_format}\n"
    with open(outpath.all_node_event_file, "w") as f:
        header = (
            "OG_ID\tNODE_ID\tGENE_NUM\tGAIN_NUM\tBRN_NUM\tDUP_NUM\t"
            + "TRN_NUM\tLOS_NUM\tTRN_DETAIL"
        )
        f.write(header + "\n")
        f.write(output_content)

    # Output results in newick mapped format
    nodeid2total_node_event: Dict[str, NodeEvent] = {}
    for node_event_list in group_id2node_event_list.values():
        # Count up the number of events for each node id
        for node_event in node_event_list:
            node_id = node_event.node_id
            if node_id not in nodeid2total_node_event.keys():
                nodeid2total_node_event[node_id] = NodeEvent(node_id, gene_num=0)
            nodeid2total_node_event[node_id] += node_event

    UtilTree.map_node_event(
        outpath.um_nodeid_tree_file,
        nodeid2total_node_event,
        outpath.all_gain_loss_map_file,
        "gain-loss",
    )
    UtilTree.map_node_event(
        outpath.um_nodeid_tree_file,
        nodeid2total_node_event,
        outpath.all_dtl_map_file,
        "dtl",
    )


def output_aggregate_transfer_results(outpath: OutPath) -> None:
    """Output aggregated transfer results"""
    trn_gene_list: List[AngstTransferGene] = []
    trn_fromto2all_gene_id_list = defaultdict(list)
    for angst_result_dir in sorted(outpath.dtl_rec_dir.glob("**/angst/")):
        og_id = angst_result_dir.parent.name
        trn_gene = AngstTransferGene(
            outpath.um_nodeid_tree_file, angst_result_dir, og_id
        )
        trn_gene_list.append(trn_gene)
        for trn_fromto, gene_id_list in trn_gene.trn_fromto2gene_id_list.items():
            trn_fromto2all_gene_id_list[trn_fromto].extend(gene_id_list)

    # Write all transfer derived gene list per line
    all_trn_genes_info = ""
    for trn_gene in trn_gene_list:
        for trn_fromto, gene_id_list in trn_gene.trn_fromto2gene_id_list.items():
            for gene_id in gene_id_list:
                all_trn_genes_info += f"{trn_gene.group_id}\t{gene_id}\t{trn_fromto}\n"
    with open(outpath.all_trn_genes_file, "w") as f:
        header = "OG_ID\tTRANSFER_DERIVED_GENE_ID\tTRANSFER_PATH\n"
        f.write(header)
        f.write(all_trn_genes_info)

    # Write transfer derived gene count summary file
    # Number of genes in descending order
    all_trn_count_info = ""
    for trn_fromto, all_gene_id_list in sorted(
        trn_fromto2all_gene_id_list.items(), key=lambda x: len(x[1]), reverse=True
    ):
        all_gene_num = len(all_gene_id_list)
        all_trn_count_info += (
            f"{trn_fromto}\t{all_gene_num}\t{'|'.join(all_gene_id_list)}\n"
        )
    with open(outpath.all_trn_count_file, "w") as f:
        header = (
            "TRANSFER_PATH\tTRANSFER_DERIVED_GENE_NUM\tTRANSFER_DERIVED_GENE_ID_LIST\n"
        )
        f.write(header)
        f.write(all_trn_count_info)


if __name__ == "__main__":
    main()
