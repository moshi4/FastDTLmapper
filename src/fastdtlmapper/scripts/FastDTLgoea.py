#!/usr/bin/env python3
import csv
import os
import shutil
import subprocess as sp
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd
from fastdtlmapper.goea import GOEA, OgGoAssociation
from fastdtlmapper.goea.args import Args, get_args
from fastdtlmapper.out_path import OutPath


def main(args: Optional[Args] = None):
    """Run GOEA(GO Enrichment Analysis) for FastDTLmapper result"""
    # Get arguments
    if args is None:
        args = get_args()
    outpath = OutPath(args.indir, goea_mode=True)

    # GO annotation using interproscan
    for fasta_file in outpath.user_fasta_dir.glob("*.fa"):
        species_name = fasta_file.stem
        annotation_outfile = outpath.go_annotation_dir / (
            species_name + "_annotation.tsv"
        )
        if not annotation_outfile.exists():
            run_interproscan(
                fasta_file, annotation_outfile, outpath.go_annotation_workdir
            )
    shutil.rmtree(outpath.go_annotation_workdir)

    # Make OG & GO association file
    annotation_file_list = list(outpath.go_annotation_dir.glob("*_annotation.tsv"))
    og_go_association = OgGoAssociation(outpath.og_gene_file, annotation_file_list)
    og_go_association.write_og2go_association(outpath.og2go_association_file)

    # Make each node gain/loss GOEA resource files
    make_node_goea_resource(
        outpath.all_node_event_file, outpath.go_enrichment_dir, "gain"
    )
    make_node_goea_resource(
        outpath.all_node_event_file, outpath.go_enrichment_dir, "loss"
    )

    # Run GOEA for each node gain/loss genes
    GOEA.download_obo(outpath.obo_file)
    run_goatools_goea(
        outpath.go_enrichment_dir,
        outpath.og2go_association_file,
        outpath.obo_file,
        args.plot_pvalue_thr,
        args.plot_max_num,
        args.plot_format,
        args.plot_color,
        args.use_adjusted_pvalue,
    )

    # Generate result summary report
    significant_go_df = pd.DataFrame()
    significant_go_count_info = ""
    for goea_result_file in sorted(outpath.go_enrichment_dir.glob("**/*.tsv")):
        filename = goea_result_file.with_suffix("").name
        node_id, gain_or_loss, go_category = filename.split("_")
        # Extract only significant data
        over_df, under_df = GOEA.extract_significant_goea_result(
            goea_result_file, args.plot_pvalue_thr, args.use_adjusted_pvalue
        )
        # Format dataframe
        over_df = GOEA.format_significant_goea_dataframe(
            over_df, node_id, gain_or_loss, go_category
        )
        under_df = GOEA.format_significant_goea_dataframe(
            under_df, node_id, gain_or_loss, go_category
        )
        # Get significant GO count stats
        over_go_num, under_go_num = len(over_df), len(under_df)
        over_go_list, under_go_list = "|".join(over_df["GO"]), "|".join(under_df["GO"])
        significant_go_count_info += f"{node_id}\t{gain_or_loss}\t{go_category}\t"
        significant_go_count_info += f"over\t{over_go_num}\t{over_go_list}\n"
        significant_go_count_info += f"{node_id}\t{gain_or_loss}\t{go_category}\t"
        significant_go_count_info += f"under\t{under_go_num}\t{under_go_list}\n"
        # Concat all significant GO dataframe
        if significant_go_df.empty:
            significant_go_df = pd.concat([over_df, under_df])
        else:
            significant_go_df = pd.concat([significant_go_df, over_df, under_df])

    # Write all significant GO dataframe
    significant_go_df.to_csv(outpath.significant_go_list_file, sep="\t", index=False)

    # Write significat GO count stats
    with open(outpath.significant_go_count_file, "w") as f:
        header = (
            "NODE_ID\tGAIN/LOSS\tGO_CATEGORY\tOVER/UNDER\t"
            + "SIGNIFICANT_GO_COUNT\tSIGNIFICANT_GO_LIST\n"
        )
        f.write(header)
        f.write(significant_go_count_info)

    # Copy all plot file in one directory
    for plot_file in outpath.go_enrichment_dir.glob(f"**/*.{args.plot_format}"):
        shutil.copy(plot_file, outpath.result_summary_plot_dir)


def run_interproscan(
    fasta_file: Path, annotation_outfile: Path, go_annotation_workdir: Path
) -> None:
    """Run InterProScan

    Args:
        fasta_file (Path): Annotation target fasta file path
        annotation_outfile (Path): Annotation result output file path
        go_annotation_workdir (Path): InterProScan working directory
    """
    ipr_bin = "interproscan.sh"
    cmd = (
        f"{ipr_bin} -i {fasta_file} -o {annotation_outfile} -f tsv "
        + f"-T {go_annotation_workdir} --goterms"
    )
    if not shutil.which(ipr_bin):
        print(f"'{ipr_bin}' not found. Please confirm InterProScan installation!!")
        exit(1)
    sp.run(cmd, shell=True)


def make_node_goea_resource(
    all_og_node_event_file: Path,
    go_enrichment_dir: Path,
    goea_type: str,
) -> None:
    """Make each node GOEA resource file

    Args:
        all_og_node_event_file (Path): All OG's node event file path
        go_enrichment_dir (Path): GOEA output directory path
        goea_type (str): "gain" or "loss"
    """
    if goea_type not in ("gain", "loss"):
        raise ValueError("goea_type must choice from 'gain' or 'loss'")

    # Make node id & node event dict
    node_id2node_event_list = defaultdict(list)
    with open(all_og_node_event_file) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            node_id = row[1]
            # Skip root node (GOEA impossible)
            if node_id == "N001":
                continue
            node_id2node_event_list[node_id].append(row)

    # Make target & all gene list files for goea
    for node_id, node_event_list in node_id2node_event_list.items():
        go_enrichment_node_dir = go_enrichment_dir / node_id
        go_enrichment_node_type_dir = go_enrichment_node_dir / goea_type
        os.makedirs(go_enrichment_node_type_dir, exist_ok=True)
        target_gene_info, all_gene_info = "", ""
        for row in node_event_list:
            og_id = row[0]
            gene_num, gain_num, los_num = int(row[2]), int(row[3]), int(row[7])

            if goea_type == "gain":
                for target_gene_cnt in range(gain_num):
                    target_gene_info += f"{og_id}_{target_gene_cnt:05d}\n"
                for all_gene_cnt in range(gene_num):
                    all_gene_info += f"{og_id}_{all_gene_cnt:05d}\n"
            elif goea_type == "loss":
                for target_gene_cnt in range(los_num):
                    target_gene_info += f"{og_id}_{target_gene_cnt:05d}\n"
                for all_gene_cnt in range(gene_num + los_num - gain_num):
                    all_gene_info += f"{og_id}_{all_gene_cnt:05d}\n"

        with open(go_enrichment_node_type_dir / "goea_target_gene_list.txt", "w") as f:
            f.write(target_gene_info)
        with open(go_enrichment_node_type_dir / "goea_all_gene_list.txt", "w") as f:
            f.write(all_gene_info)


def run_goatools_goea(
    go_enrichment_dir: Path,
    association_file: Path,
    obo_file: Path,
    pvalue_thr: float = 0.05,
    plot_max_num: int = 10,
    plot_format: str = "png",
    plot_color: str = "",
    use_adjusted_pvalue: bool = False,
) -> None:
    """Run goatools GOEA

    Args:
        go_enrichment_dir (Path): GOEA output directory path
        association_file (Path): Gene ID & GO association file path
        obo_file (Path): GO OBO file path
        pvalue_thr (float): Plot GOterm pvalue threshold
        plot_max_num (int): Max GOterm plot number
        plot_format (str): Plot file format
        plot_color (str): Plot specified color
        use_adjusted_pvalue (bool): Use adjusted pvalue or not
    """
    for goea_node_dir in sorted(go_enrichment_dir.glob("*/*/")):
        # Get GOEA target & all list files
        target_gene_list_file = goea_node_dir / "goea_target_gene_list.txt"
        all_gene_list_file = goea_node_dir / "goea_all_gene_list.txt"
        if os.stat(target_gene_list_file).st_size == 0:
            continue
        # Run GOEA
        goea = GOEA(
            target_gene_list_file,
            all_gene_list_file,
            association_file,
            obo_file,
            pvalue_thr,
            plot_max_num,
            plot_format,
            plot_color,
            use_adjusted_pvalue,
        )
        node_id = goea_node_dir.parent.name
        gain_or_loss = goea_node_dir.name
        output_prefix = goea_node_dir / f"{node_id}_{gain_or_loss}"
        goea_result_file_list = goea.run(output_prefix)
        # Plot GOEA significant GOterms
        for goea_result_file in goea_result_file_list:
            goea.plot(goea_result_file, goea_result_file.with_suffix(""))


if __name__ == "__main__":
    main()
