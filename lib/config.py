import argparse
import os
import random
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List


@dataclass
class Config:
    """Config Class"""

    indir: Path
    tree_file: Path
    outdir: Path
    process_num: int
    dup_cost: int
    los_cost: int
    trn_cost: int
    rseed: int

    def __post_init__(self):
        # TODO: Fix directory structure
        # 00_user_data -> dir:fasta, tree
        # 01_orthofinder (Raw result)
        # 02_dtl_reconciliation -> dir:OGname
        #    -> file: OGfasta, aln, trim,  tree, rooted_tree, dir: mowgli
        # 03_aggregate_dtl_result -> file: all_dtl_map, all_gl_map, all_node_event_tsv
        # log -> dir: dependent programe names -> file: command.txt, output.txt
        self.fasta_dir = self.outdir / "00_fasta"
        self.ortho_dir = self.outdir / "01_orthofinder"
        self.ortho_aln_dir = self.outdir / "02_ortho_align"
        self.ortho_aln_trim_dir = self.outdir / "03_ortho_align_trim"
        self.gene_tree_dir = self.outdir / "04_gene_tree"
        self.rooted_gene_tree_dir = self.outdir / "05_rooted_gene_tree"
        self.dtl_dir = self.outdir / "06_dtl_reconciliation"
        self.dtl_event_map_dir = self.outdir / "07_dtl_event_map"

        self.ortho_tmpwork_dir = self.fasta_dir / "OrthoFinder"

        self.user_tree_file = self.outdir / "user_tree.nwk"
        self.ultrametric_tree_file = self.outdir / "user_ultrametric_tree.nwk"
        self.nodeid_tree_file = self.outdir / "user_ultrametric_nodeid_tree.nwk"

        self.ortho_group_fasta_dir = self.ortho_dir / "Orthogroup_Sequences"

        self.parallel_cmds_tmp_file = self.outdir / "parallel_cmds_tmp.txt"

        self._makedirs()
        self._add_bin_path()

    def _makedirs(self) -> None:
        """Create config output directories"""
        # Check & Delete OrthoFinder previous work dir
        if self.ortho_tmpwork_dir.exists():
            shutil.rmtree(self.ortho_tmpwork_dir)

        # Make configured output directory
        for target_dir in (
            self.fasta_dir,
            self.ortho_dir,
            self.ortho_aln_dir,
            self.ortho_aln_trim_dir,
            self.gene_tree_dir,
            self.rooted_gene_tree_dir,
            self.dtl_dir,
            self.dtl_event_map_dir,
        ):
            os.makedirs(target_dir, exist_ok=True)

    def _add_bin_path(self) -> None:
        """Add bin programs path"""
        # Bin path list
        bin_path = Path(__file__).parent.parent / "bin"
        mafft_path = bin_path / "mafft"
        ortho_finder_path = bin_path / "OrthoFinder"
        ortho_finder_tool_path = ortho_finder_path / "tools"
        # Add bin path
        env_path = os.environ["PATH"]
        env_path = f"{bin_path}:{env_path}"
        env_path = f"{mafft_path}:{env_path}"
        env_path = f"{ortho_finder_path}:{env_path}"
        env_path = f"{ortho_finder_tool_path}:{env_path}"
        # Set fixed path
        os.environ["PATH"] = env_path

    def make_ultrametric_cmd(self, nwk_tree_infile: str, nwk_tree_outfile: str) -> str:
        """Get make_ultrametric run command"""
        return (
            f"make_ultrametric.py {nwk_tree_infile} {nwk_tree_outfile} "
            + "> /dev/null 2>&1"
        )

    def orthofinder_cmd(self, fasta_indir: str) -> str:
        """Get OrthoFinder run command"""
        return f"orthofinder.py -og -f {fasta_indir} -t {self.process_num}"

    def mafft_cmd(self, fasta_infile: str, aln_outfile: str) -> str:
        """Get mafft run command"""
        return f"mafft --auto --quiet {fasta_infile} > {aln_outfile}"

    def trimal_cmd(self, aln_infile: str, aln_trim_outfile: str) -> str:
        """Get trimal run command"""
        return (
            f"trimal -in {aln_infile} -out {aln_trim_outfile} -automated1 "
            + "> /dev/null 2>&1"
        )

    def fasttree_cmd(self, aln_infile: str, tree_outfile: str) -> str:
        """Get FastTree run command"""
        return f"FastTree -quiet -seed {self.rseed} {aln_infile} > {tree_outfile}"

    def optroot_cmd(self, gene_tree_infile: str, rooted_gene_tree_outfile: str) -> str:
        """Get OptRoot run command"""
        return (
            f"OptRoot -q -i {gene_tree_infile} -o {rooted_gene_tree_outfile} "
            + f"-D {self.dup_cost} -L {self.los_cost} -T {self.trn_cost} "
            + f"-r --seed {self.rseed} --type 2 > /dev/null 2>&1"
        )

    def mowgli_cmd(
        self, species_tree_infile: str, gene_tree_infile, outdir: str
    ) -> str:
        """Get Mowgli run command

        Note:
            Mowgli cannot handle node id tagged tree (e.g. N001)
        """
        return (
            f"Mowgli -s {species_tree_infile} -g {gene_tree_infile} -o {outdir} "
            + f"-d={self.dup_cost} -l={self.los_cost} -t={self.trn_cost} "
            + "-n 1 -T 80"
        )

    def parallel_cmd(self, cmd_list: List[str]) -> None:
        """Get parallel run command from command list"""
        # TODO: Logging parallel command list
        # Write parallel cmd list file
        random.shuffle(cmd_list)
        with open(self.parallel_cmds_tmp_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        return (
            f"parallel --no-notice --bar --eta -a {self.parallel_cmds_tmp_file} "
            + f"-j {self.process_num} "
            + f"&& rm {self.parallel_cmds_tmp_file}"
        )


def get_config() -> Config:
    """Get config from arguments

    Returns:
        Config: Config Class
    """
    parser = argparse.ArgumentParser(
        description="Fast Genome-wide DTL event mapping tool"
    )

    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=Path,
        help="Input Fasta(*.fa|*.faa|*.fasta), Genbank(*.gb|*.gbk|*.genbank) directory",
        metavar="",
    )
    parser.add_argument(
        "-t",
        "--tree_file",
        required=True,
        type=Path,
        help="Input rooted species time(ultrametric) tree file (Newick format)",
        metavar="",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory",
        metavar="",
    )
    default_processor_num = os.cpu_count() - 1
    parser.add_argument(
        "-p",
        "--process_num",
        type=int,
        help=f"Number of processor (Default: {default_processor_num})",
        default=default_processor_num,
        metavar="",
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
    default_rseed = 0
    parser.add_argument(
        "--rseed",
        type=int,
        help=f"Number of random seed (Default: {default_rseed})",
        default=default_rseed,
        metavar="",
    )

    args = parser.parse_args()

    return Config(
        args.indir,
        args.tree_file,
        args.outdir,
        args.process_num,
        args.dup_cost,
        args.los_cost,
        args.trn_cost,
        args.rseed,
    )
