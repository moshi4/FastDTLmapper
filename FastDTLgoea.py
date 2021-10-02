#!/usr/bin/env python3
import argparse
import csv
import os
import re
import shutil
import subprocess as sp
from collections import defaultdict
from pathlib import Path
from typing import Tuple

import pandas as pd

from lib.goea import GOEA, OgGoAssociation


def main(
    indir: Path,
    plot_pvalue_thr: float,
    plot_max_num: int,
    plot_format: str,
    plot_color: str,
    use_adjusted_pvalue: bool,
):
    """Run GOEA(GO Enrichment Analysis) for FastDTLmapper result

    Args:
        indir (Path): FastDTLmapper result directory
        plot_pvalue_thr (float): Plot GOterm pvalue threshold
        plot_max_num (int): Plot GOterm max number
        plot_format (str): Plot file format
        plot_color (str): Plot specified hexcolor
        use_adjusted_pvalue (bool): Use adjusted pvalue or not
    """
    # Output directory
    outdir = indir / "04_functional_analysis"
    go_annotation_dir = outdir / "go_annotation"
    go_annotation_workdir = go_annotation_dir / "work"
    go_enrichment_dir = outdir / "go_enrichment"
    result_summary_dir = outdir / "result_summary"
    result_summary_plot_dir = result_summary_dir / "significant_go_plot"
    shutil.rmtree(go_enrichment_dir, ignore_errors=True)
    shutil.rmtree(result_summary_dir, ignore_errors=True)
    os.makedirs(go_annotation_dir, exist_ok=True)
    os.makedirs(go_annotation_workdir, exist_ok=True)
    os.makedirs(go_enrichment_dir, exist_ok=True)
    os.makedirs(result_summary_plot_dir, exist_ok=True)

    # GO annotation using interproscan
    fasta_dir = indir / "00_user_data" / "fasta"
    if not fasta_dir.exists():
        err_msg = f"Input fasta directory '{fasta_dir}' not found!!\n"
        err_msg += "Please specify FastDTLmapper result directory as input directory."
        raise ValueError(err_msg)

    for fasta_file in fasta_dir.glob("*.fa"):
        species_name = fasta_file.stem
        annotation_outfile = go_annotation_dir / (species_name + "_annotation.tsv")
        if not annotation_outfile.exists():
            run_interproscan(fasta_file, annotation_outfile, go_annotation_workdir)
    shutil.rmtree(go_annotation_workdir)

    # Make OG & GO association file
    og_gene_file = indir / "01_orthofinder" / "Orthogroups" / "Orthogroups.txt"
    annotation_file_list = go_annotation_dir.glob("*_annotation.tsv")
    og2go_association_file = go_enrichment_dir / "og2go_association.txt"
    og_go_association = OgGoAssociation(og_gene_file, annotation_file_list)
    og_go_association.write_og2go_association(og2go_association_file)

    # Make each node gain/loss GOEA resource files
    all_og_node_event_file = indir / "03_aggregate_map_result" / "all_og_node_event.tsv"
    make_node_goea_resource(all_og_node_event_file, go_enrichment_dir, "gain")
    make_node_goea_resource(all_og_node_event_file, go_enrichment_dir, "loss")

    # Run GOEA for each node gain/loss genes
    obo_file = outdir / "go-basic.obo"
    GOEA.download_obo(obo_file)
    run_goatools_goea(
        go_enrichment_dir,
        og2go_association_file,
        obo_file,
        plot_pvalue_thr,
        plot_max_num,
        plot_format,
        plot_color,
        use_adjusted_pvalue,
    )

    # Generate result summary report
    significant_go_df = pd.DataFrame()
    significant_go_count_info = ""
    for goea_result_file in sorted(go_enrichment_dir.glob("**/*.tsv")):
        filename = goea_result_file.with_suffix("").name
        node_id, gain_or_loss, go_category = filename.split("_")
        # Extract only significant data
        over_df, under_df = extract_significant_goea_result(
            goea_result_file, plot_pvalue_thr, use_adjusted_pvalue
        )
        # Format dataframe
        over_df = format_significant_go_dataframe(
            over_df, node_id, gain_or_loss, go_category
        )
        under_df = format_significant_go_dataframe(
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
    significant_go_list_file = result_summary_dir / "significant_go_list.tsv"
    significant_go_df.to_csv(significant_go_list_file, sep="\t", index=False)

    # Write significat GO count stats
    significant_go_count_file = result_summary_dir / "significant_go_count.tsv"
    with open(significant_go_count_file, "w") as f:
        header = (
            "NODE_ID\tGAIN/LOSS\tGO_CATEGORY\tOVER/UNDER\t"
            + "SIGNIFICANT_GO_COUNT\tSIGNIFICANT_GO_LIST\n"
        )
        f.write(header)
        f.write(significant_go_count_info)

    # Copy all plot file in one directory
    for plot_file in go_enrichment_dir.glob(f"**/*.{plot_format}"):
        shutil.copy(plot_file, result_summary_plot_dir)


def get_args() -> argparse.Namespace:
    """Get argument values

    Returns:
        argparse.Namespace: argument values
    """
    parser = argparse.ArgumentParser(
        description="Gain/Loss genes GOEA(GO Enrichment Analysis) tool"
    )
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="FastDTLmapper result directory",
        metavar="IN",
    )
    default_plot_pvalue_thr = 0.05
    parser.add_argument(
        "--plot_pvalue_thr",
        type=float,
        default=default_plot_pvalue_thr,
        help=f"Plot GOterm pvalue threshold (Default: {default_plot_pvalue_thr})",
        metavar="",
    )
    default_max_plot_num = 10
    parser.add_argument(
        "--plot_max_num",
        type=int,
        default=default_max_plot_num,
        help=f"Plot GOterm max number (Default: {default_max_plot_num})",
        metavar="",
    )
    default_plot_format = "png"
    available_format = ["png", "svg", "jpg", "pdf"]
    parser.add_argument(
        "--plot_format",
        type=str,
        default=default_plot_format,
        choices=available_format,
        help=f"Plot file format [{'|'.join(available_format)}] "
        + f"(Default: '{default_plot_format}')",
        metavar="",
    )
    parser.add_argument(
        "--plot_color",
        type=str,
        default="",
        help="Plot specified hexcolor [e.g. '1affdb'] "
        + "(Default: yellow to red gradient color)",
        metavar="",
    )
    parser.add_argument(
        "--adjusted_pvalue",
        help="Use BH adjusted pvalue for plot threshold",
        action="store_true",
    )

    args = parser.parse_args()

    # Plot hexcolor check
    def is_valid_hexcolor(hexcolor: str) -> bool:
        return re.search(r"^(?:[0-9a-fA-F]{3}){1,2}$", hexcolor)

    if args.plot_color and not is_valid_hexcolor(args.plot_color):
        parser.error(
            f"argument --plot_color: invalid hexcolor code '{args.plot_color}'."
        )
    else:
        args.plot_color = "#" + args.plot_color if args.plot_color else ""

    return args


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
    if not shutil.which(ipr_bin):
        err = f"'{ipr_bin}' not found. Please confirm InterProScan installation!!"
        raise RuntimeError(err)

    cmd = (
        f"{ipr_bin} -i {fasta_file} -o {annotation_outfile} -f tsv "
        + f"-T {go_annotation_workdir} --goterms"
    )
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


def extract_significant_goea_result(
    goea_result_file: Path,
    pvalue_thr: float,
    use_adjusted_pvalue: bool,
    min_depth: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract over and under significant goea result

    Args:
        goea_result_file (Path): GOEA result file
        pvalue_thr (float): Pvalue threshold for extract
        use_adjusted_pvalue (bool): Use BH adjusted pvalue or not
        min_depth (int, optional): Minimum depth for extract. Defaults to 2.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: over/under extract dataframes
    """
    df = pd.read_table(goea_result_file)
    pvalue_column_name = "p_fdr_bh" if use_adjusted_pvalue else "p_uncorrected"
    over_df = df[
        (df["enrichment"] == "e")
        & (df[pvalue_column_name] < pvalue_thr)
        & (df["depth"] >= min_depth)
    ]
    under_df = df[
        (df["enrichment"] == "p")
        & (df[pvalue_column_name] < pvalue_thr)
        & (df["depth"] >= min_depth)
    ]
    return (over_df, under_df)


def format_significant_go_dataframe(
    goea_result_df: pd.DataFrame,
    node_id: str,
    gain_or_loss: str,
    go_category: str,
) -> pd.DataFrame:
    """Format significant GOterm dataframe for output

    Args:
        goea_result_df (pd.DataFrame): GOEA result dataframe
        node_id (str): Node id
        gain_or_loss (str): "gain" or "loss"
        go_category (str): "BP" or "MF" or "CC"

    Returns:
        pd.DataFrame: Formatted GOEA result dataframe
    """
    # Rename columns
    rename_list = [
        "GO",
        "GO_CATEGORY",
        "OVER/UNDER",
        "GO_NAME",
        "RATIO_IN_STUDY",
        "RATIO_IN_POP",
        "PVALUE",
        "DEPTH",
        "STUDY_COUNT",
        "BH_ADJUSTED_PVALUE",
        "STUDY_ITEMS",
    ]
    rename_dict = {b: a for b, a in zip(goea_result_df.columns, rename_list)}
    goea_result_df = goea_result_df.rename(columns=rename_dict)
    # Add columns
    goea_result_df["NODE_ID"] = node_id
    goea_result_df["GAIN/LOSS"] = gain_or_loss
    goea_result_df["GO_CATEGORY"] = go_category
    # Replace OVER/UNDER value ("e" -> "over", "p" -> "under")
    goea_result_df = goea_result_df.replace({"OVER/UNDER": {"e": "over", "p": "under"}})
    # Return reorder columns dataframe
    return goea_result_df[
        ["NODE_ID", "GAIN/LOSS", "GO_CATEGORY", "OVER/UNDER", "GO", *rename_list[3:]]
    ]


if __name__ == "__main__":
    args = get_args()
    main(
        args.indir,
        args.plot_pvalue_thr,
        args.plot_max_num,
        args.plot_format,
        args.plot_color,
        args.adjusted_pvalue,
    )
