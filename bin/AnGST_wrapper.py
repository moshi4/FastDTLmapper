#!/usr/bin/env python3
import argparse
import re
import shutil
import subprocess as sp
from pathlib import Path


def main(
    species_tree_file: Path,
    gene_tree_file: Path,
    outdir: Path,
    dup_cost: int,
    los_cost: int,
    trn_cost: int,
    timetree: bool,
):
    """AnGST wrapper run function

    Args:
        species_tree_file (Path): Species tree file path
        gene_tree_file (Path): Gene tree file path
        outdir (Path): Output directory
        dup_cost (int): Duplication event cost
        los_cost (int): Loss event cost
        trn_cost (int): Transfer event cost
        timetree (bool): Timetree or not flag
    """
    # Remove previous directory
    shutil.rmtree(outdir)

    # Make work directory & files
    work_dir = outdir.parent / "work"
    work_species_tree_file = work_dir / "species_tree.nwk"
    work_gene_tree_file = work_dir / "gene_tree.nwk"
    work_penalty_score_file = work_dir / "penalty_score.txt"
    work_input_info_file = work_dir / "input_info.txt"
    work_log_file = work_dir / "log.txt"
    work_dir.mkdir(exist_ok=True)
    shutil.copy(species_tree_file, work_species_tree_file)
    sp.run(f"cat {gene_tree_file} | head -n 100 > {work_gene_tree_file}", shell=True)

    # Fix species tree file for AnGST run
    # AnGST requires branch length in all branch
    with open(work_species_tree_file) as f:
        replace_tree_info = "(" + f.read().replace(";", ":0.01);")
    with open(work_species_tree_file, "w") as f:
        f.write(replace_tree_info)

    # Make penalty score file
    with open(work_penalty_score_file, "w") as f:
        f.write(f"hgt: {trn_cost}\ndup: {dup_cost}\nlos: {los_cost}\nspc: 0.0\n")

    # Make AnGST input info file
    input_info = ""
    input_info += f"species={work_species_tree_file}\n"
    input_info += f"gene={work_gene_tree_file}\n"
    input_info += f"penalties={work_penalty_score_file}\n"
    input_info += f"output={outdir}\n"
    input_info += "ultrametric=True\n" if timetree else "ultrametric=False\n"
    with open(work_input_info_file, "w") as f:
        f.write(input_info)

    # Run AnGST
    angst_run_cmd = f"AnGST.py {work_input_info_file} > {work_log_file}"
    sp.run(angst_run_cmd, shell=True)

    # Fix AnGST unreadable newick format chimeric gene tree
    angst_nwk_file = outdir / "AnGST.newick"
    with open(angst_nwk_file) as f:
        nwk_text = f.readline().replace(");", "")[1:]
        # match = re.search(r"^\((.+\))[^\)]+\)", nwk_text)
        match = re.search(r"^\(.+\)[^\)]+\)", nwk_text)
        fix_nwk_info = match.group() + ";"
    with open(angst_nwk_file, "w") as f:
        f.write(fix_nwk_info)

    # Remove work directory
    shutil.rmtree(work_dir)


def get_args() -> argparse.Namespace:
    """Get argument values

    Returns:
        argparse.Namespace: Argument values
    """
    parser = argparse.ArgumentParser(description="AnGST run wrapper tool")

    parser.add_argument(
        "-s",
        "--species_tree",
        required=True,
        type=Path,
        help="Input species tree file",
        metavar="IN",
    )
    parser.add_argument(
        "-g",
        "--gene_tree",
        required=True,
        type=Path,
        help="Input gene tree file",
        metavar="TREE",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="OUT",
    )
    default_dup_cost = 2
    parser.add_argument(
        "--dup_cost",
        type=int,
        help=f"Duplication event cost (Default: {default_dup_cost})",
        default=default_dup_cost,
        metavar="",
    )
    default_los_cost = 1
    parser.add_argument(
        "--los_cost",
        type=int,
        help=f"Loss event cost (Default: {default_los_cost})",
        default=default_los_cost,
        metavar="",
    )
    default_trn_cost = 3
    parser.add_argument(
        "--trn_cost",
        type=int,
        help=f"Transfer event cost (Default: {default_trn_cost})",
        default=default_trn_cost,
        metavar="",
    )
    parser.add_argument(
        "--timetree",
        help="Use species tree as timetree",
        action="store_true",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    main(
        args.species_tree,
        args.gene_tree,
        args.outdir,
        args.dup_cost,
        args.los_cost,
        args.trn_cost,
        args.timetree,
    )
