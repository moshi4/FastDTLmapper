#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import subprocess as sp
from collections import defaultdict
from pathlib import Path

from lib.goea import GOEA, OgGoAssociation


def main(
    indir: Path,
    plot_pvalue_thr: float,
    plot_max_num: int,
    plot_format: str,
    use_adjusted_pvalue: bool,
):
    """Run GOEA(GO Enrichment Analysis) for FastDTLmapper result

    Args:
        indir (Path): FastDTLmapper result directory
        plot_pvalue_thr (float): Plot GOterm pvalue threshold
        plot_max_num (int): Plot GOterm max number
        plot_format (str): Plot file format
        use_adjusted_pvalue (bool): Use adjusted pvalue or not
    """
    # Output directory
    outdir = indir / "04_go_enrichment"
    go_annotation_dir = outdir / "go_annotation"
    go_annotation_workdir = go_annotation_dir / "work"
    go_enrichment_dir = outdir / "go_enrichment"
    os.makedirs(go_annotation_dir, exist_ok=True)
    os.makedirs(go_annotation_workdir, exist_ok=True)
    os.makedirs(go_enrichment_dir, exist_ok=True)

    # GO annotation using interproscan
    fasta_dir = indir / "00_user_data" / "fasta"
    if not fasta_dir.exists():
        raise ValueError(f"Input directory '{fasta_dir}' not found!!")

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
    obo_file = go_enrichment_dir / "go-basic.obo"
    GOEA.download_obo(obo_file)
    run_goatools_goea(
        go_enrichment_dir,
        og2go_association_file,
        obo_file,
        plot_pvalue_thr,
        plot_max_num,
        plot_format,
        use_adjusted_pvalue,
    )


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
        metavar="",
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
        + f"(Default: {default_plot_format})",
        metavar="",
    )
    parser.add_argument(
        "--adjusted_pvalue",
        help="Use BH adjusted pvalue for plot threshold",
        action="store_true",
    )

    return parser.parse_args()


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
            use_adjusted_pvalue,
        )
        output_prefix = goea_node_dir / f"{goea_node_dir.name}_goea"
        goea_result_file_list = goea.run(output_prefix)
        # Plot GOEA significant GOterms
        for goea_result_file in goea_result_file_list:
            goea.plot(goea_result_file, goea_result_file.with_suffix(""))


if __name__ == "__main__":
    args = get_args()
    main(
        args.indir,
        args.plot_pvalue_thr,
        args.plot_max_num,
        args.plot_format,
        args.adjusted_pvalue,
    )
